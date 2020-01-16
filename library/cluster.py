"""
Produce a hierarchical tree (e.g. UPGMA tree)
from a distance matrix of genomes.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""
import numpy
import subprocess
import os
import pickle
from multiprocessing import Manager, Pool, Process

from scipy.cluster.hierarchy import linkage, fcluster

from library.sketch import MASH_LOC
from library.genome_selection.strategy import filter_genomes

SKETCH_END = ".msh"

def dist(sketch_1, sketch_2, sketch_prog=MASH_LOC):
    cp = subprocess.run([
        sketch_prog, "dist", sketch_1, sketch_2
    ], capture_output=True)
    return (cp.stdout.decode().strip().split(), 
            cp.returncode)
    
    
def cluster(taxid, sketch_locs):
    """
    
    Parameters
    ----------
    sketch_locs: dict
        A mapping from accessions to sketch locations.
        
    """
    n = len(sketch_locs)
    genome_pos_map = {}
    pos_genome_map = {}
    for i, genome in enumerate(sketch_locs.keys()):
        genome_pos_map[genome] = i
        pos_genome_map[i] = genome
    distances = numpy.zeros((int(n * (n-1) / 2)))
    # distances = numpy.zeros((n, n))
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            out, ret = dist(sketch_locs[pos_genome_map[i]], 
                            sketch_locs[pos_genome_map[j]])
            if not ret:
                return None
            distances[c] = float(out[2])
            # distances[j][i] = float(out[2])
            c += 1
    return taxid, linkage(distances, method='average'), genome_pos_map


def get_sketch_locs(root_dir, taxid, index):
    genome_list = list(index["taxids"][taxid].keys())
    genome_list, _ = filter(genome_list, index)
    sketch_locs = {}
    for acc in genome_list:
        loc = os.path.join(root_dir, index["genomes"][acc]['location'])
        files = (file for file in os.listdir(loc) if os.path.isfile(os.path.join(loc, file)))
        for f in files:
            if f.endswith(SKETCH_END):
                sketch_locs[acc] = os.path.join(loc, f)
                break
    return sketch_locs


def cluster_all(root_dir, root_taxid, tree, index, output_dir, processes=1):
    def traverse(taxid):
        children = tree[taxid]
        if not children:
            return [taxid]
        taxid_list = []
        for child in children:
            taxid_list += traverse(child)
        return taxid_list
    species_list = traverse(root_taxid)
    pool = Pool(processes=processes)
    pd_list = []
    for taxid in species_list:
        sketch_locs = get_sketch_locs(root_dir, taxid, index)
        pd_list.append(pool.apply_async(cluster, args=(taxid, sketch_locs)))
    for pd in pd_list:
        taxid, Z, genome_pos_map = pd.get()
        try:
            os.makedirs(os.path.join(output_dir, str(taxid)))
        except FileExistsError:
            pass
        numpy.save(os.path.join(output_dir, str(taxid), str(taxid) + "_Z"), Z)
        with open(os.path.join(output_dir, str(taxid), str(taxid) + "_map.pck"), "wb") as genome_file: 
            pickle.dump(genome_pos_map, genome_file)
    return len(pd_list)
    