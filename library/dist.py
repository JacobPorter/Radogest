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

from library.genome_selection.strategy import filter_genomes, select_genomes
from library.sketch import MASH_LOC

# The ending of sketch files.
SKETCH_END = ".msh"

# The maximum distance matrices to store in memory at once.
# If this number is high, then the process may use a lot more memory.
MAX_JOBS = 64


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
    (str, int)
        A string representing the distance between sketch_1 and sketch_2

        An int representing the return code of sketch_prog

    """
    cp = subprocess.run([sketch_prog, "dist", sketch_1, sketch_2],
                        capture_output=True)
    return (cp.stdout.decode().strip().split()[2], cp.returncode)


def distance_cell(cell_id, loc_i, loc_j):
    """
    Calculate the pairwise distance between genomes.

    Parameters
    ----------
    cell_id: tuple, str, etc.
        Some identifier for the pairwise distance.
    loc_i: str
        The location of a genome.
    loc_j: str
        The location of a genome.

    Returns
    -------
    ret: int
        The return code for calculating the distance.
    cell_id: tuple, str, etc.
        Some identifier for the pairwise distance.
        The same as the input.
    dist: float
        The distance between genomes.

    """
    mash_dist, ret = dist(loc_i, loc_j)
    if ret:
        return ret, cell_id, 0.0
    return ret, cell_id, float(mash_dist)


def init_array_mapping(genome_list):
    """
    Gets the genome to pos mapping and the initial distance matrix.

    Parameters
    ---------
    genome_list: list
        A list of genome accessions.

    Returns
    -------
    n: int
        The number of genomes to get sketches for.
    distances: array
        A zero array to store distance values.
    genome_pos_map: dict
        A mapping from genomes to positions in the distances matrix.
    pos_genome_map: dict
        A mapping from positions in the distance matrix to genomes.

    """
    n = len(genome_list)
    if n == 1:
        return None
    if n == 0:
        return None
    genome_pos_map = {}
    pos_genome_map = {}
    for i, genome in enumerate(genome_list):
        genome_pos_map[genome] = i
        pos_genome_map[i] = genome
    distances = numpy.zeros((n, n))
    return n, distances, genome_pos_map, pos_genome_map


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
        return (1, None, None, None)
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
                distances[c] = float(out)
            else:
                distances[j][i] = float(out)
                distances[i][j] = float(out)
            c += 1
    return 0, taxid, distances, pos_genome_map


def get_sketch_locs(root_dir,
                    taxid,
                    index,
                    down_select="random",
                    select_amount=None):
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
    if not (down_select == "none" or select_amount <= 0
            or select_amount is None):
        genome_list = select_genomes(genome_list,
                                     index,
                                     down_select=down_select,
                                     select_amount=select_amount)
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


def dist_all(root_dir,
             root_taxid,
             tree,
             index,
             output_dir,
             down_select="random",
             select_amount=None,
             processes=1):
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
    down_select: str
        The type of down selection to use.  random, sort, none
    select_amount: int
        The number of genomes to select when down_select is not 'none'.
    processes: int
        The number of parallel processes to use.


    Returns
    -------
    int, int
        The number of attempted distance matrices computed followed
        by the number of distance matrices that failed to compute.

    """
    def traverse(taxid):
        """A depth first search traversal of the taxonomic tree
        to find all species taxonomic ids."""
        children = tree[taxid]
        if not children:
            return [taxid]
        taxid_list = []
        for child in children:
            taxid_list += traverse(child)
        return taxid_list

    def assign_d(cell_id, d_list, taxid_d_map, dist):
        """Assign a value to a cell in a matrix."""
        taxid, i, j = cell_id
        d = d_list[taxid_d_map[taxid]]
        d[i][j] = dist
        d[j][i] = dist

    def write_dist_map(taxid_d_map, d_list, map_list):
        """Write a distance matrix."""
        written = 0
        for taxid in taxid_d_map:
            X = d_list[taxid_d_map[taxid]]
            mapping = map_list[taxid_d_map[taxid]]
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
                pickle.dump(mapping, genome_file)
            written += 1
        return written

    print("Root directory {}".format(root_dir), file=sys.stderr)
    root_taxid = int(root_taxid)
    species_list = traverse(root_taxid)
    print("Found {} species: {}".format(len(species_list), species_list),
          file=sys.stderr)
    pool = Pool(processes=processes)
    pd_list = []
    total_written = 0
    failed = 0
    d_list = []
    map_list = []
    jobs = 0
    taxid_d_map = {}
    for taxid in species_list:
        sketch_locs = get_sketch_locs(root_dir, taxid, index, down_select,
                                      select_amount)
        init = init_array_mapping(list(sketch_locs.keys()))
        if not init:
            failed += 1
            continue
        n, d, _, pos_genome_map = init
        # d_list.append((d, 0, (n * (n -1 ) / 2)))
        d_list.append(d)
        map_list.append(pos_genome_map)
        taxid_d_map[taxid] = len(d_list) - 1
        for i in range(n):
            loc_i = sketch_locs[pos_genome_map[i]]
            for j in range(i + 1, n):
                loc_j = sketch_locs[pos_genome_map[j]]
                pd_list.append(
                    pool.apply_async(distance_cell,
                                     args=((taxid, i, j), loc_i, loc_j)))
        jobs += 1
        if jobs == MAX_JOBS:
            for pd in pd_list:
                _, cell_id, mash_dist = pd.get()
                assign_d(cell_id, d_list, taxid_d_map, mash_dist)
            total_written += write_dist_map(taxid_d_map, d_list, map_list)
            jobs = 0
            pd_list = []
            d_list = []
            map_list = []
            taxid_d_map = {}
    for pd in pd_list:
        _, cell_id, mash_dist = pd.get()
        assign_d(cell_id, d_list, taxid_d_map, mash_dist)
    total_written += write_dist_map(taxid_d_map, d_list, map_list)
    return total_written, failed


# def dist_all(root_dir, root_taxid, tree, index, output_dir,
#              down_select="random", select_amount=None, processes=1):
#     """
#     Compute distance matrices for all species under
#     the given taxonomic id and save them to disk.
#
#     Parameters
#     ----------
#     root_dir: str
#         The root directory where genomes are contained.
#     root_taxid: int
#         The taxonomic id that determines the subtree under which
#         all distance matrices will be calculated.
#     tree: dict
#         The Radogest tree data structure.
#     index: dict
#         The Radogest index data structure.
#     output_dir: str
#         The location of a directory to store the distance matrices.
#     processes: int
#         The number of parallel processes to use.
#
#     Returns
#     -------
#     int, int
#         The number of attempted distance matrices computed followed
#         by the number of distance matrices that failed to compute.
#
#     """
#     def traverse(taxid):
#         """A depth first search traversal of the taxonoic tree
#         to find all species taxonomic ids."""
#         children = tree[taxid]
#         if not children:
#             return [taxid]
#         taxid_list = []
#         for child in children:
#             taxid_list += traverse(child)
#         return taxid_list
#
#     print("Root directory {}".format(root_dir), file=sys.stderr)
#     root_taxid = int(root_taxid)
#     species_list = traverse(root_taxid)
#     print("Found {} species: {}".format(len(species_list), species_list),
#           file=sys.stderr)
#     pool = Pool(processes=processes)
#     pd_list = []
#     failed = 0
#     for taxid in species_list:
#         sketch_locs = get_sketch_locs(root_dir, taxid, index,
#                                       down_select, select_amount)
#         if len(sketch_locs) > 1:
#             pd_list.append(
#                 pool.apply_async(distance_matrix, args=(taxid, sketch_locs)))
#     for pd in pd_list:
#         out = pd.get()
#         if not out[0]:
#             _, taxid, X, pos_genome_map = out
#             try:
#                 os.makedirs(os.path.join(output_dir, str(taxid)))
#             except FileExistsError:
#                 pass
#             numpy.save(os.path.join(output_dir, str(taxid),
#                                     str(taxid) + "_X"), X)
#             with open(
#                     os.path.join(output_dir, str(taxid),
#                                  str(taxid) + "_map.pck"),
#                     "wb") as genome_file:
#                 pickle.dump(pos_genome_map, genome_file)
#         else:
#             failed += 1
#     return len(pd_list) - failed, failed
