"""
Save a distance matrix of genomes for each species.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""
import os
import pickle
import subprocess
import sys
from multiprocessing import Pool

import numpy

from library.genome_selection.strategy import filter_genomes
from library.sketch import MASH_LOC

# The ending of sketch files.
SKETCH_END = ".msh"


def dist(sketch_1, sketch_2, sketch_prog=MASH_LOC):
    """
    Compute a distance between sketched genomes.

    Parameters
    ----------
    sketch_1: str
        The location of a sketch file.
    sketch_2: str
        The location of another sketch file.
    sketch_prog: str
        The location of a program for calculating distances.
        (e.g. mash)

    Returns
    -------
    (list, int)
        A list of [sketch_1, sketch_2, distance, p-value, jaccard_ratio]
        An int representing the return code of sketch_prog

    """
    cp = subprocess.run([sketch_prog, "dist", sketch_1, sketch_2],
                        capture_output=True)
    return (cp.stdout.decode().strip().split(), cp.returncode)


def distance_matrix(taxid, sketch_locs, square=True):
    """
    Create a distance matrix of genomes and a mapping of positions to genomes.

    Parameters
    ----------
    taxid: int
        The taxid corresponding to genomes to calculate distance
    sketch_locs: dict
        A mapping from accessions to sketch locations.
    square: boolean
        If True, return a symmetric square matrix,
        if False, return a flat distance matrix.

    Returns
    -------
    (int, int, arr, dict)
        int A return code
        int The taxid
        arr The distance matrix
        dict A mapping from positions to genomes.

    """
    n = len(sketch_locs)
    if n == 1:
        # Is it better to return an empty matrix?
        return (1, None, None, None)
#         return (0, taxid, 
#                 numpy.zeros((1,1)), 
#                 {0: list(sketch_locs.keys())[0]})
    elif n == 0:
        return (1, None, None, None)
    genome_pos_map = {}
    pos_genome_map = {}
    for i, genome in enumerate(sketch_locs.keys()):
        genome_pos_map[genome] = i
        pos_genome_map[i] = genome
    if not square:
        distances = numpy.zeros((int(n * (n - 1) / 2)))
    else:
        distances = numpy.zeros((n, n))
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            out, ret = dist(sketch_locs[pos_genome_map[i]],
                            sketch_locs[pos_genome_map[j]])
            if ret:
                return ret, out
            if not square:
                distances[c] = float(out[2])
            else:
                distances[j][i] = float(out[2])
                distances[i][j] = float(out[2])
            c += 1
    return 0, taxid, distances, pos_genome_map


def get_sketch_locs(root_dir, taxid, index):
    """
    Get a mapping of accessions to sketch file locations.

    Parameters
    ----------
    root_dir: str
        The location of the root directory that contains genomes.
    taxid: int
        The taxonomic id that determines which genomes to find.
        All (unfiltered) genomes in the subtree rooted at taxid
        will be returned.
    index: dict
        A Radogest index.

    Returns
    -------
    dict
        A mapping from accessions to sketch file locations.

    """
    genome_list = list(index["taxids"][taxid].keys())
    genome_list, _ = filter_genomes(genome_list, index)
    sketch_locs = {}
    for acc in genome_list:
        loc = os.path.join(root_dir, index["genomes"][acc]['location'][1:])
        files = (file for file in os.listdir(loc)
                 if os.path.isfile(os.path.join(loc, file)))
        for f in files:
            if f.endswith(SKETCH_END):
                sketch_locs[acc] = os.path.join(loc, f)
                break
    return sketch_locs


def dist_all(root_dir, root_taxid, tree, index, output_dir, processes=1):
    """
    Compute distance matrices for all species under
    the given taxonomic id and save them to disk.

    Parameters
    ----------
    root_dir: str
        The root directory where genomes are contained.
    root_taxid: int
        The taxonomic id that determines the subtree under which
        all distance matrices will be calculated.
    tree: dict
        The Radogest tree data structure.
    index: dict
        The Radogest index data structure.
    output_dir: str
        The location of a directory to store the distance matrices.
    processes: int
        The number of parallel processes to use.

    Returns
    -------
    int, int
        The number of attempted distance matrices computed followed
        by the number of distance matrices that failed to compute.

    """
    def traverse(taxid):
        """A depth first search traversal of the taxonoic tree
        to find all species taxonomic ids."""
        children = tree[taxid]
        if not children:
            return [taxid]
        taxid_list = []
        for child in children:
            taxid_list += traverse(child)
        return taxid_list

    print("Root directory {}".format(root_dir), file=sys.stderr)
    root_taxid = int(root_taxid)
    species_list = traverse(root_taxid)
    print("Found {} species: {}".format(len(species_list), species_list),
          file=sys.stderr)
    pool = Pool(processes=processes)
    pd_list = []
    failed = 0
    for taxid in species_list:
        sketch_locs = get_sketch_locs(root_dir, taxid, index)
        if len(sketch_locs) > 1:
            pd_list.append(
                pool.apply_async(distance_matrix, args=(taxid, sketch_locs)))
    for pd in pd_list:
        out = pd.get()
        if not out[0]:
            _, taxid, X, pos_genome_map = out
            try:
                os.makedirs(os.path.join(output_dir, str(taxid)))
            except FileExistsError:
                pass
            numpy.save(os.path.join(output_dir, str(taxid),
                                    str(taxid) + "_X"), X)
            with open(
                    os.path.join(output_dir, str(taxid),
                                 str(taxid) + "_map.pck"),
                    "wb") as genome_file:
                pickle.dump(pos_genome_map, genome_file)
        else:
            failed += 1
    return len(pd_list), failed
