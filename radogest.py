#!/usr/bin/env python
"""
Radogest: random genome sampler for trees.

Commands:
download
make_index
update_index
down_select
make_tree
split_fasta
sample
parallel_sample

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import datetime
import sys
import pickle
import os
import json
import subprocess
import random
import string
import sqlite3
import errno
import shutil

from random import shuffle
from multiprocessing import Pool
from collections import defaultdict

from library.ncbi_genome_download.ncbi_genome_download.core import argument_parser as ncbi_argument_parser
from library.ncbi_genome_download.ncbi_genome_download import __version__ as ncbi_version
from library.ncbi_genome_download.ncbi_genome_download.__main__ import run_ncbi

from tqdm import tqdm
from ete3 import NCBITaxa

# from library.genome_selection.strategy import ProportionalRandom
# from library.genome_selection.strategy import QualitySortingLeaf
# from library.genome_selection.strategy import QualitySortingTree

from SeqIterator.SeqIterator import SeqWriter
from SeqIterator.SeqIterator import SeqReader

ncbi = NCBITaxa()

FASTA_ENDINGS = ['fasta', 'fa', 'fna', 'faa']

# get_paths
"""Stores the paths to tools and the data."""

SAMTOOLS = "/home/jsporter/Applications/samtools-1.8/"
UCSC = "/home/jsporter/Applications/UCSC/"
BEDTOOLS = "/home/jsporter/Applications/bedtools2/bin/"
GENOMES = "/groups/fungcat/datasets/current/fasta/Genomes/"
GENOMESAA = "/groups/fungcat/datasets/current/fasta/GenomesAA/"

# If the sample falls below NUMBER_CUTOFF, get NUMBER_SAMPLE samples instead.
NUMBER_CUTOFF = 150
NUMBER_SAMPLE = 300
# The actual number of random samples taken is multipled by this multiplier.
# This is done so that samples with N's in them can be excluded.
SAMPLE_MULTIPLIER = 1.2

def get_paths():
    """Give the locations of the paths."""
    return {'SAMTOOLS': SAMTOOLS,
            'UCSC': UCSC,
            'BEDTOOLS': BEDTOOLS,
            'GENOMES': GENOMESAA}

# Tools and genomes locations.
path_dict = get_paths()

GENOMES_LOCATION = path_dict['GENOMES']
BEDTOOLS_DIR = path_dict['BEDTOOLS']
UCSC_DIR = path_dict['UCSC']
# 
# SAMTOOLS = path_dict['SAMTOOLS']
# UCSC = path_dict['UCSC']


def add_index(index, genome_info, verbose=True):
    """
    Add the contents of JSON genome information to an index.

    Parameters
    ----------
    index: dictionary
        a dictionary representing an index of taxonomic ids to genomes
    genome_info: iterable of dictionaries
        each record represents a genome
    verbose: boolean
        if set to True, prints out periodic counts of records processed.

    Returns
    -------
    int
        A count of the records processed.

    """
    ncbi = NCBITaxa()
    count = 0
    for genome_record in tqdm(genome_info):
        group = genome_record[1]
        genome_record = genome_record[0]
        accession = genome_record["assembly_accession"]
        index['genomes'][accession][
            'location'] = '/{section}/{group}/{accession}/'.format(
                section=genome_record['section'],
                group=group, accession=accession)
        for key in genome_record:
            index['genomes'][accession][key] = genome_record[key]
        for taxid in ncbi.get_lineage(int(genome_record["taxid"])):
            if accession not in index['taxids'][taxid]:
                index['taxids'][taxid][accession] = True
        # if verbose and count % COUNT_THRESHOLD == 0:
        #    sys.stderr.write("Processed {} records.\n".format(count))
        count += 1
    return count

"""
Run faidx from samtools on all genomes below the given directory.
Decompress gzip comressed fasta files.  Create an index
of genomic and taxonomic information.
This program requires SAMTOOLS and UCSC tools.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def process_fai(fai_file_fd):
    """
    Use an fai file to update the index.

    Parameters
    ----------
    fai_file_fd: file object
        An open file object for the fai file.

    Returns
    -------
    contigs, genome_length
        The number of contigs followed by the number of nucleotide bases.

    """
    contigs = 0
    genome_length = 0
    for line in fai_file_fd:
        contigs += 1
        genome_length += int(line.split('\t')[1])
    return contigs, genome_length


def name_ends(name, endings, addition=""):
    """
    Determine if a string ends with something.

    Parameters
    ----------
    name: str
        The string to determine if it ends with something.
    endings: iterable
        A list of endings to check.
    addition: str
        A constant string to search for that is added to all endings.

    Returns
    -------
    True if name ends with something in endings plus the addition.

    """
    for end in endings:
        if name.endswith(end + addition):
            return True
    return False


def make_fai_individual(root_directory, path, files, twoBitToFa=False,
                        leave_compressed=False, verbose=0, index_only=False):
    """
    Make the fai files and decompress the fasta files.

    Parameters
    ----------
    root_directory: str
        A string giving the location of the root directory where genomes
        are found.
    path: str
        The path to a directory containing files
    files: iterable
        An iterable of files under the path.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.

    Returns
    -------
    count, contig, genome_length
        A dictionary that gives which types of files were processed.
        The number of contigs in the fasta file that was processed.
        The number of bases in the fasta file.

    """
    count = defaultdict(int)
    contigs = None
    genome_length = None
    fasta_exists = False
    for f in files:
        if (name_ends(f, FASTA_ENDINGS, addition=".gz") or
                name_ends(f, FASTA_ENDINGS, addition="")):
            fasta_exists = True
            break
    if not fasta_exists:
        return None, None, None
    # if twoBitToFa:
    #     skip = False
    #     for name in files:
    #         if (name_ends(name, FASTA_ENDINGS)):
    #             skip = True
    #             break
    #     if skip:
    #         for name in files:
    #             if name.endswith('.2bit'):
    #                 os.remove(os.path.join(path, name))
    #         return count, contigs, genome_length
    if not index_only:
        for name in files:
            if (name.endswith('fai')):
                os.remove(os.path.join(path, name))
    for name in files:
        if index_only:
            if name.endswith('.fai'):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                fai_fd = open(os.path.join(path, name))
                contigs, genome_length = process_fai(fai_fd)
                count["fai"] += 1
            else:
                continue
        else:
            if (not leave_compressed and
                    name_ends(name, FASTA_ENDINGS, addition='.gz')):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                subprocess.run(["gunzip", os.path.join(path, name)])
                name = name[0:len(name) - 3]
                count["gzip"] += 1
            if name.endswith('.2bit'):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                new_name = name[0:len(name)-5]
                fa_file = os.path.join(path, new_name)
                twobit_process = subprocess.run([UCSC + "twoBitToFa",
                                                 os.path.join(path, name),
                                                 fa_file])
                if twobit_process.returncode == 0:
                    os.remove(os.path.join(path, name))
                name = new_name
                count["2bit"] += 1
            if (name_ends(name, FASTA_ENDINGS) and not twoBitToFa):
                fa_file = os.path.join(path, name)
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                subprocess.run([SAMTOOLS + "samtools", "faidx", fa_file])
                count["fai"] += 1
                contigs, genome_length = process_fai(
                    open(os.path.join(path, name + '.fai')))
    return count, contigs, genome_length


def make_fai_multi(index, root_directory, twoBitToFa=False,
                   leave_compressed=False, verbose=0,
                   index_only=False, workers=16):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.
    workers: int
        The number of workers to use for multiprocessing.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    Examples
    --------
        make_fai({}, '/home/jsporter/Genomes/')

    """
    count = {"gzip": 0, "2bit": 0, "fai": 0}
    counter = 0
    pool = Pool(processes=workers)
    pd_list = []
    for path, _, files in os.walk(root_directory):
        # print(path, files)
        pd = pool.apply_async(make_fai_individual, (root_directory,
                                                    path,
                                                    files,
                                                    twoBitToFa,
                                                    leave_compressed,
                                                    verbose,
                                                    index_only))
        pd_list.append(pd)
    pool.close()
    pool.join()
    for pd in pd_list:
        count_local, contigs, genome_length = pd.get()
        accession = os.path.basename(path)
        if contigs is not None and genome_length is not None:
            index['genomes'][accession]['contig_count'] = contigs
            index['genomes'][accession]['base_length'] = genome_length
        else:
            print(path, files)
        for key in count_local:
            count[key] += count_local[key]
        counter += 1
        if counter >= 5000 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
    sys.stderr.flush()
    return count, index


def remove_accessions(index, accessions_to_remove, verbose=0):
    """
    Remove genomes from the index.

    Parameters
    ----------
    index: dict
        The genomes index.
    accessions_to_remove: dict
        A dict keyed on accessions to remove from the index.

    Returns
    -------
    dict
        The genomes index is returned.

    """
    remove_count = 0
    print("Removing empty genome folders from the index.", file=sys.stderr)
    count = 0
    for taxid in index['taxids']:
        my_accessions_to_remove = []
        for accession in index['taxids'][taxid]:
            if accession in accessions_to_remove:
                my_accessions_to_remove.append(accession)
        for accession in my_accessions_to_remove:
            del index['taxids'][taxid][accession]
            remove_count += 1
        count += 1
        if verbose >= 1 and count % 5000 == 0:
            print("Processed: {}".format(count), file=sys.stderr)
            sys.stderr.flush()
    return index, remove_count


def make_fai(index, root_directory, twoBitToFa=False,
             leave_compressed=False, verbose=0,
             index_only=False):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    output_fd: writable (such as a file)
        A writable object such as an open file to periodically save the index.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    Examples
    --------
        make_fai({}, '/home/jsporter/Genomes/')

    """
    count = {"gzip": 0, "2bit": 0, "fai": 0}
    counter = 0
    accessions_to_remove = {}
    for path, _, files in os.walk(root_directory):
        counter += 1
        if counter % 5000 == 0 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
        count_local, contigs, genome_length = make_fai_individual(
            root_directory, path, files, twoBitToFa,
            leave_compressed, verbose, index_only)
        accession = os.path.basename(path)
        if count_local is None:
            accessions_to_remove[accession] = True
        else:
            if contigs is not None and genome_length is not None:
                index['genomes'][accession]['contig_count'] = contigs
                index['genomes'][accession]['base_length'] = genome_length
            for key in count_local:
                count[key] += count_local[key]
    index, remove_count = remove_accessions(index, accessions_to_remove,
                                            verbose=verbose)
    count['remove_count'] = remove_count
    sys.stderr.flush()
    return count, index

"""
Construct a taxonomic rank tree hierarchy from the genomes index.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def count_levels(tree, roots):
    """
    Count the number of taxid that are not leaves under the taxid under roots.

    Parameters
    ----------
    tree: dict
        The tree hierarchy.
    roots: list
        A list of taxonomic ids to search under.

    Returns
    -------
    dict
        A dictionary indexed on a taxonomic id root from roots.
        The first number is the number of taxonomic ids under root.
        The second number is the number of taxonomic ids under root that
        have only a single child.

    """
    counter_dict = {root: [0, 0] for root in roots}
    taxid_list = []
    for taxid in tree:
        lineage = ncbi.get_lineage(taxid)
        for root in roots:
            if root in lineage:
                taxid_list.append(taxid)
                counter_dict[root][0] = counter_dict[root][0] + 1
                if len(tree[taxid]) == 1:
                    counter_dict[root][1] = counter_dict[root][1] + 1
                break
    return counter_dict, taxid_list


def remove_only_children(tree, roots):
    """
    Remove nodes with only one child.  Reconnect parent nodes to their
    new children.

    Parameters
    ----------
    tree: dict
        The taxonomic tree.
    roots: iterable
        The root nodes in the tree to traverse down.

    Returns
    -------
    int
        The number of nodes removed.

    """
    nodes_deleted = 0
    for root in roots:
        deletion_list = []
        remapping = defaultdict(list)
        remove_only_children_traversal(tree, root, deletion_list, remapping)
        for taxid in remapping:
            for old_child, new_child in remapping[taxid]:
                tree[taxid].remove(old_child)
                tree[taxid].append(new_child)
        for taxid in deletion_list:
            del tree[taxid]
        nodes_deleted += len(deletion_list)
    return nodes_deleted


def remove_only_children_traversal(tree, node, deletion_list, remapping):
    """
    Traverse the taxonomic tree (depth-first search) to identify nodes
    with only one child.  Gives a remapping from parents to children.

    Parameters
    ----------
    tree: dict
        The taxonomic tree.
    node: int
        The rooth node to start searching.
    deletion_list: list
        A list to store taxonomic ids to delete.
    remapping: defaultdict(list)
        A mapping between parents and their [(old_child, new_child), ..., ]

    Returns
    -------
    (int, int)
        A taxonomic node and its number of children.

    """
    children = tree[node]
    if len(children) == 0:
        return (node, 0)
    for child in children:
        ret_node, amount = remove_only_children_traversal(tree, child,
                                                          deletion_list,
                                                          remapping)
        if ret_node != child and len(children) != 1:
            remapping[node].append((child, ret_node))
    if len(children) == 1:
        deletion_list.append(node)
        return ret_node, amount
    else:
        return (node, len(children))


def check_leaves_species(taxid, tree):
    """
    Validate that the leaves of the tree are all species.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id.
    tree: dict
        A dictionary representing the taxonomic tree.

    Returns
    -------
    list
        A list of taxonomic ids that are leaves and that are not species.

    """
    if not tree[taxid] and 'species' not in ncbi.get_rank([taxid])[taxid]:
        return [taxid]
    elif not tree[taxid]:
        return []
    list_of_taxid = []
    for child in tree[taxid]:
        list_of_taxid + check_leaves_species(child, tree)
    return list_of_taxid


def make_tree(index, verbose=False):
    """
    Construct a tree hierarchy using conventional ranks from the genomes index.

    Parameters
    ----------
    index: dict
        The all genomes index
    verbose: boolean
        Print additional output.

    Returns
    -------
    (dict, dict)
        1. A dictionary object that represents the subranks below each rank.
        2. A dictionary that gives a list of missing ranks for species.

    """
    rank_list = ["superkingdom", "kingdom", "phylum", "class", "order",
                 "family", "genus", "species"]
    missing_ranks = defaultdict(list)
    looked_at = {}
    # ranks = {"superkingdom":defaultdict(set), "kingdom":defaultdict(set),
    #          "phylum":defaultdict(set), "class":defaultdict(set),
    #          "order":defaultdict(set), "family":defaultdict(set),
    #          "genus":defaultdict(set)}
    ranks = defaultdict(set)
    print("Creating the initial tree.", file=sys.stderr)
    for accession in tqdm(index['genomes']):
        try:
            taxid = int(index['genomes'][accession]['taxid'])
        except KeyError:
            print("Could not find a key for: {}".format(accession),
                  file=sys.stderr)
            raise KeyError
        if taxid in looked_at:
            continue
        if not index['taxids'][taxid]:
            continue
        lineage_rank = ncbi.get_rank(ncbi.get_lineage(taxid))
        lineage_rank = {lineage_rank[key]: key for key in lineage_rank if
                        lineage_rank[key] != 'no rank'}
        looked_at[taxid] = True
        text_translator = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
        exclude = False
        for text in text_translator.values():
            if (len(lineage_rank) <= 2 or
                    # 'unclassified' in text or  # Include unclassified?
                    'uncultured' in text or
                    'RNA virus' in text or
                    'satellite RNA' in text or
                    'environmental samples' in text):
                exclude = True
                break
        if exclude:
            continue
        for i, rank in enumerate(rank_list):
            if i == len(rank_list) - 1:
                break
            if rank not in lineage_rank:
                continue
            for j in range(i+1, len(rank_list)):
                if rank_list[j] in lineage_rank:
                    ranks[lineage_rank[rank]].add(
                        lineage_rank[rank_list[j]])
                    break
    for taxid in ranks:
        ranks[taxid] = list(ranks[taxid])
    unclassified = []
    print("Removing species with 'unclassified' in their lineage.",
          file=sys.stderr)
    for parent in tqdm(ranks):
        for taxid in ranks[parent]:
            if 'species' in ncbi.get_rank([taxid])[taxid]:
                text_translator = ncbi.get_taxid_translator(
                    ncbi.get_lineage(taxid))
                exclude = False
                for text in text_translator.values():
                    if 'unclassified' in text.lower():
                        if verbose:
                            print(text_translator, file=sys.stderr)
                        exclude = True
                        break
                if exclude:
                    unclassified.append(taxid)
    unclassified_to_remove = len(unclassified)
    for taxid in ranks:
        for species in unclassified:
            if species in ranks[taxid]:
                ranks[taxid].remove(species)
    # Add the root node for eukaryotes, bacteria, viruses, archaea
    ranks[1] = [2759, 2, 10239, 2157]
    print("Removing nodes with only one child.", file=sys.stderr)
    children_removed = remove_only_children(ranks, [1, 12884])
    print("Removing empties.", file=sys.stderr)
    empties_remove = 0
    to_delete = []
    for taxid in ranks:
        if not ranks[taxid]:
            to_delete.append(taxid)
            empties_remove += 1
    for taxid in to_delete:
        del ranks[taxid]
    print("Counting the levels and producing the taxonomic id list.",
          file=sys.stderr)
    counter_dict, taxid_list = count_levels(ranks, ranks[1])
    taxid_list = [1] + taxid_list
    shuffle(taxid_list)
    return (ranks, missing_ranks,
            unclassified_to_remove,
            empties_remove, children_removed, counter_dict, taxid_list)


"""
Split a list of genomes into smaller sequences.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def split_genomes(accessions_list, length, index, output,
                  includeNs=False, window_length=50):
    """
    Partition the DNA strings in a list of genomes and write to a fasta file.

    Parameters
    ----------
    accessions_list: iterable
        A list of genome accessions.
    length: int
        The length of the partition to take.
    index: dict
        The genomes index.
    output: str or writable
        The location of the output file as a string path or a writable object.

    Returns
    -------
    int
        The number of records written.

    """
    if window_length < 0:
        window_length = length
    if isinstance(output, str):
        if not output:
            output = sys.stdout
        else:
            output = open(output, 'w')
    if isinstance(output, SeqWriter):
        writer = output
    else:
        writer = SeqWriter(output, file_type='fasta')
    number_written = 0
    for accession in accessions_list:
        location = GENOMES + index['genomes'][accession]['location']
        # print(location, file=sys.stderr)
        only_files = [f for f in os.listdir(location) if
                      os.path.isfile(os.path.join(location, f))
                      and f.endswith('fna')]
        if len(only_files) == 0:
            print("Skipping {}.  Fasta file not found.".format(accession),
                  file=sys.stderr)
            continue
        location = os.path.join(location, only_files[0])
        reader = SeqReader(location, file_type='fasta')
        for header, sequence in reader:
            for i in range(0, len(sequence), window_length):
                substring = sequence[i:i+length].upper()
                if (len(substring) == length and (
                        includeNs or 'N' not in substring)):
                    writer.write(("{}_{}_[{}:{}]".format(
                        accession, header, i, i+length), substring))
                    number_written += 1
    return number_written


#!/usr/bin/env python
"""
Produces a fasta file of randomly chosen sequences corresponding to the
given taxonomic id.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def binary_search(target, array, begin, end):
    """
    Binary search

    Parameters
    ----------
    target: comparable
        An item to compare to
    array: iterable, indexable
        An array to search for
    begin: int
        Beginning index for the array
    end: int
        Ending index for the array

    Returns
    -------
    An item that satisfies the binary search.

    """
    if array == []:
        return None
    if begin == end:
        return begin
    if begin == end-1:
        return begin if target < array[begin] else end
    mid = begin + int(round((end - begin) / 2.0))
    if target < array[mid]:
        return binary_search(target, array, begin, mid)
    return binary_search(target, array, mid, end)


def uniform_samples_at_rank(taxid, index, sublevels, number):
    """
    Get a count of samples from each genome under the taxids given by ranks.

    Parameters
    ----------
    taxid: int
        A taxonomic id number to query from
    index: dict
        A dictionary representing genome locations and genome taxid
        associations
    sublevels: iterable
        An iterable that gives ranks underneath taxid to sample from.
    number: int
        The number of random sequences to sample.

    Returns
    -------
    dict
        A dictionary keyed on taxids where the value is the number of sequences
        to sample for each taxid.

    """
    try:
        uniform_number = round(number / len(sublevels))
    except ZeroDivisionError:
        return None
    uniform_sample_counts = []
    for level in sublevels:
        uniform_sample_counts.append((level, uniform_samples(level, index,
                                                             uniform_number)))
    return uniform_sample_counts


def uniform_samples(taxid, index, number):
    """
    Get a count of samples from each genome under taxid.  Does not do any
    sampling.

    Parameters
    ----------
    taxid: int
        A taxonomic id number to query from
    index: dict
        A dictionary representing genome locations and genome taxid
        associations
    number: int
        The number of random sequences to sample.

    Returns
    -------
    dict
        A dictionary keyed on accession ids.  The value is the number of
        sequences to sample from the corresponding genome.

    """
    accessions = [accession for
                  accession in index['taxids'][taxid] if
                  index['taxids'][taxid][accession]]
    # if len(accessions) > GENOMES_TO_KEEP:
    #     random.shuffle(accessions)
    #     accessions = accessions[:GENOMES_TO_KEEP]
    genome_lengths = [index['genomes'][accession]['base_length']
                      for accession in
                      accessions]
    try:
        name = ncbi.get_taxid_translator([taxid])[taxid]
    except KeyError:
        name = ""
    except sqlite3.DatabaseError:
        name = ""
    sys.stderr.write("For {}:{}, sampling from "
                     "{} genomes.\n".format(taxid, name, len(accessions)))
    sys.stderr.flush()
    for i, _ in enumerate(genome_lengths):
        if i == 0:
            continue
        genome_lengths[i] += genome_lengths[i-1]
    try:
        total_length = genome_lengths[-1]
    except IndexError as ie:
        print("Genome lengths are empty for taxid: {}".format(taxid),
              file=open("/home/jsporter/Sampling_Out/index_error_"+
                        str(taxid)+ ".out", 'w'))
        raise ie
    accession_counts = defaultdict(int)
    for _ in range(number):
        i = binary_search(random.randint(0, total_length - 1),
                          genome_lengths, 0, len(genome_lengths) - 1)
        accession_counts[accessions[i]] += 1
    return accession_counts


def get_fasta(accession_counts_list, length, index,
              output, taxid_file, window_length=50, verbose=False,
              thresholding=False, amino_acid=False, temp_dir='/tmp/'):
    """
    Save randomly sampled sequences in a fasta file written to output.

    Parameters
    ----------
    accession_counts_list: list
        A list of tuples.  The first element is the taxid associated with the
        genomes in the second element.  The second element is a dictionary
        keyed on accession number, and the value is the number of samples to
        take from that genome.
    length: int
        The number of bases to sample.  (The length of the string.)
    index: dict
        A dictionary representing information about genomes and taxids.
    output: writable
        A writable object to store fasta files too.  This could be an open
        file.
    taxid_file: writable
        A file-like object to store taxonomic ids that correspond to each
        sampled sequence.  Each taxonomic id will be on its own line.
    tmp: str
        A path to a directory to store temporary files.
    verbose: bool
        If True, print messages.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.


    Returns
    -------
    int
        A count of the number of fasta records written.

    """
    final_file = SeqWriter(output, file_type='fasta')
    fasta_record_count = 0
    for taxid, accession_counts in accession_counts_list:
        for accession in accession_counts.keys():
            if accession_counts[accession] <= 0:
                continue

            rand_string = "".join(random.choices(
                string.ascii_letters + string.digits, k=8))
            my_fasta = os.path.join(temp_dir,
                                    accession + "_" +
                                    rand_string +  "_random.fa")
            if verbose:
                sys.stderr.write("Writing fasta records for {}:\t".
                                 format(accession))
            # A thresholding feature.  If the genome is too small, use the
            # whole genome.
            if thresholding and accession_counts[accession] > float(
                index['genomes'][accession]['base_length']) / length:
                records_written = split_genomes([accession], length,
                                                index, final_file)
                for _ in range(records_written):
                    taxid_file.write(str(taxid) + "\n")
                fasta_record_count += records_written
                continue
            accession_location = os.path.join(GENOMES_LOCATION +
                                              index['genomes']
                                              [accession]
                                              ['location'])
            onlyfiles = [f for f in os.listdir(accession_location) if
                         os.path.isfile(os.path.join(accession_location, f))]
            for f in onlyfiles:
                if (f.endswith(".fna") or
                        f.endswith(".fasta") or
                        f.endswith(".fa") or
                        f.endswith(".faa")):
                    fasta_location = os.path.join(accession_location, f)
                elif f.endswith(".2bit"):
                    twobit_location = os.path.join(accession_location, f)
                elif f.endswith(".fai"):
                    fai_location = os.path.join(accession_location, f)
            bedtools_file = os.path.join(temp_dir, accession + "_" +
                                         rand_string + "_random.bed")
            get_random = True
            accession_cnt = 0
            while get_random:
                accession_number = accession_counts[accession] - accession_cnt
                bed_2bit_counts = get_random_bed_fast(accession_number,
                                                      length,
                                                      taxid,
                                                      accession,
                                                      fai_location,
                                                      bedtools_file,
                                                      fasta_location,
                                                      my_fasta,
                                                      taxid_file,
                                                      final_file,
                                                      amino_acid)
                get_random = not bed_2bit_counts[0]
                accession_cnt += bed_2bit_counts[1]
                if verbose:
                    sys.stderr.write(str(bed_2bit_counts[1]) +
                                     "  " +
                                     str(bed_2bit_counts[2] /
                                         (bed_2bit_counts[1] +
                                          bed_2bit_counts[2])) +
                                     " ")
                    sys.stderr.flush()
            fasta_record_count += accession_cnt
            if verbose:
                sys.stderr.write("\n")
                sys.stderr.flush()
    return fasta_record_count


def get_random_bed_2bit(number, length, taxid, accession, fai_location,
                        bedtools_file, twobit_location, my_fasta,
                        taxid_file, final_file):
    """
    Get random nucleotide sequences from a bed file and a 2bit file.  Exclude
    sequences with N's in them.

    Parameters
    ----------
    number: int
        The number of samples to take.
    length: int
        The length of the nucleotide sequence to get.
    taxid: int
        The taxonomic id to sample from
    accession: str
        A string indicating the accession id for a genome
    fai_location: str
        The location of the .fai faidx samtools index for the genome
    bedtoods_file: str
        The location to store the bedtools output file too.
    twobit_location: str
        The location of the 2bit genome file
    my_fasta: str
        The location to store a temporary fasta file for random sampling.
    taxid_file: writable
        A writable object where taxids are written for each random sample.
    final_file: writable
        A writable object that stores the fasta records.  This file represents
        the end product of all of the random sampling.

    Returns
    -------
    (boolean, int, int)
        A boolean to indicated if the records written is equal to number.
        An integer counting the records written,
        and an integer of the records that had indeterminate N characters.

    """
    raise NotImplementedError
    taxid = str(taxid)
    if number <= NUMBER_CUTOFF:
        my_sample = NUMBER_SAMPLE
    else:
        my_sample = int(number * SAMPLE_MULTIPLIER)
    bedtools_fd = open(bedtools_file, 'w')
    subprocess.run([BEDTOOLS_DIR + "bedtools", "random", "-l",
                    str(length), "-n",
                    str(my_sample), "-g",
                    fai_location], stdout=bedtools_fd)
    subprocess.run([UCSC_DIR + "twoBitToFa", "-bed="+bedtools_file,
                    "-bedPos", twobit_location, my_fasta])
    intermediate_fasta_file = SeqReader(my_fasta, file_type='fasta')
    records_with_n = 0
    records_written = 0
    records_malformed = 0
    for fasta_record in intermediate_fasta_file:
        if records_written >= number:
            break
        record_id, record_seq = fasta_record
        record_seq = record_seq.upper()
        if ">" in record_seq:
            records_malformed += 1
            continue
        if "N" in record_seq:
            records_with_n += 1
            continue
        record_id = accession + ":" + taxid + ":" + record_id
        final_file.write((record_id, record_seq))
        taxid_file.write(taxid + "\n")
        records_written += 1
    bedtools_fd.close()
    os.remove(bedtools_file)
    intermediate_fasta_file.close()
    os.remove(my_fasta)
    return (records_written >= number, records_written, records_with_n)


def get_random_bed_fast(number, length, taxid, accession, fai_location,
                        bedtools_file, fasta_location, my_fasta,
                        taxid_file, final_file, amino_acid=False):
    """
    Get random nucleotide sequences from a bed file and a fasta file.  Exclude
    sequences with N's in them.

    Parameters
    ----------
    number: int
        The number of samples to take.
    length: int
        The length of the nucleotide sequence to get.
    taxid: int
        The taxonomic id to sample from
    accession: str
        A string indicating the accession id for a genome
    fai_location: str
        The location of the .fai faidx samtools index for the genome
    bedtoods_file: str
        The location to store the bedtools output file too.
    fasta_location: str
        The location of the fasta genome file
    my_fasta: str
        The location to store a temporary fasta file for random sampling.
    taxid_file: writable
        A writable object where taxids are written for each random sample.
    final_file: writable
        A writable object that stores the fasta records.  This file represents
        the end product of all of the random sampling.

    Returns
    -------
    (boolean, int, int)
        A boolean to indicated if the records written is equal to number.
        An integer counting the records written,
        and an integer of the records that had indeterminate N characters.

    """
    taxid = str(taxid)
    if number <= NUMBER_CUTOFF:
        my_sample = NUMBER_SAMPLE
    else:
        my_sample = int(number * SAMPLE_MULTIPLIER)
    bedtools_fd = open(bedtools_file, 'w')
    subprocess.run([BEDTOOLS_DIR + "bedtools", "random", "-l",
                    str(length), "-n",
                    str(my_sample), "-g",
                    fai_location], stdout=bedtools_fd)
    subprocess.run([BEDTOOLS_DIR + "bedtools", "getfasta", "-fi",
                    fasta_location, "-bed", bedtools_file, "-fo", my_fasta])
    intermediate_fasta_file = SeqReader(my_fasta, file_type='fasta')
    records_with_n = 0
    records_written = 0
    for fasta_record in intermediate_fasta_file:
        if records_written >= number:
            break
        record_id, record_seq = fasta_record
        record_seq = record_seq.upper()
        if not amino_acid and "N" in record_seq:
            records_with_n += 1
            continue
        record_id = accession + ":" + taxid + ":" + record_id
        final_file.write((record_id, record_seq))
        taxid_file.write(taxid + "\n")
        records_written += 1
    bedtools_fd.close()
    os.remove(bedtools_file)
    intermediate_fasta_file.close()
    os.remove(my_fasta)
    return (records_written >= number, records_written, records_with_n)

"""
Randomly permutes a fasta file and a taxid file.  Can split the data into
training, validation, and testing partitions.

:Authors:
        Jacob Porter <jsporter@vt.edu>
"""

def write_fasta_taxid(both_records, fasta_writable, taxid_writable):
    """
    Write a fasta file and a taxid file from a list of both taxid and fasta
    records.

    Parameters
    ----------
    both_records: list
        A list of taxid records and fasta records where an entry in the list
        has the form (fasta_record, taxid)
    fasta_writable: writable
        A writable object where a fasta file can be written.
    taxid_writable: writable
        A writable object where a list of taxids can be written.

    Returns
    -------
    int
        A count of records written.

    """
    fasta_writer = SeqWriter(fasta_writable)
    taxid_writer = taxid_writable
    count = 0
    for record in both_records:
        fasta_writer.write(record[0])
        taxid_writer.write(record[1] + "\n")
        count += 1
    fasta_writer.flush()
    taxid_writer.flush()
    fasta_writer.close()
    taxid_writer.close()
    return count


def randomly_permute_fasta_taxid(fasta_file, taxid_file, fasta_out,
                                 taxid_out, split=True, split_amount=0.8):
    """
    Randomly permute the fasta and taxid files.  May split the input into
    training, validation, and testing data.

    Parameters
    ----------
    fasta_file: str
        The location of a fasta file.
    taxid_file: str
        The location of a taxid file.
    fasta_out: str
        The location to store the output fasta file.
    taxid_out: str
        The location to store the output taxid file.
    split: bool
        If True, split the data.
    split_amount: tuple of floats or float
        Percentage(s) to split into.

    Returns
    -------
    int
        The number of records written.

    """
    both_records = []
    for fasta, taxid in zip(SeqReader(fasta_file,
                                      file_type='fasta'),
                            open(taxid_file, 'r')):
        both_records.append((fasta, taxid.strip()))
    shuffle(both_records)
    print("The total number of records to write: {}".format(len(both_records)),
          file=sys.stderr)
    if split:
        print("Splitting the output.", file=sys.stderr)
        amounts = []
        if isinstance(split_amount, float):
            amounts.append(int(split_amount * len(both_records)))
        else:
            # print(list(map(float, split_amount.split(","))))
            # print(len(both_records))
            if sum(list(map(float, split_amount.split(",")))) > 1.0:
                print("The sum of the split amounts is larger than 1.0",
                      file=sys.stderr)
            for amount in list(map(float, split_amount.split(","))):
                amounts.append(round(amount * len(both_records)))
        begin = 0
        count = 0
        extensions = [".train", ".validate", ".test"]
        for i, amount in enumerate(amounts):
            ext = extensions[i] if i < len(extensions) else str(i)
            end = (begin+amount if
                   begin+amount <= len(both_records) else
                   len(both_records))
            count += write_fasta_taxid(both_records[begin:end],
                                       open(fasta_out + ext, 'w'),
                                       open(taxid_out + ext, 'w'))
            begin = end
        return count
    return write_fasta_taxid(both_records,
                             open(fasta_out, 'w'),
                             open(taxid_out, 'w'))


"""
Runs multiple taxonomic id sampling instances in parallel.
Creates training, validation,
and testing data and puts them in directories that Plinko expects.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def get_sample(taxid, sublevels, index_dir, number, length, data_dir,
               split_amount='0.8,0.1,0.1', thresholding=False,
               window_length=100, amino_acid=False, temp_dir="/tmp/"):
    """
    Get a random sample.  Create training, validation, and testing data sets
    and put them in the appropriate folders.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id
    rank: str
        A rank in the taxonomic system.  Examples: genus, family, etc.
    sublevels: iterable
        A set of taxonomic ids below the taxid to sample from.
    index_dir: str
        A path to the pickled genomes index object.
    number: int
        The number of samples to take.
    length: int
        The number of bases for each sample
    data_dir: str
        The path to the data directory where fasta files will be written.
    split_amount: str
        A comma seperated list of floats representing the percentage of the
        data to be used for training, validation, and testing data
    tmp_dir: str
        A path to write temporary files to.  This could be on the local hard
        drive for better speed.

    Returns
    -------
    (int, int)
        A tuple of fasta records sampled and permuted records written.

    """
    index = pickle.load(open(index_dir, 'rb'))
    accession_counts = uniform_samples_at_rank(taxid, index, sublevels, number)
    if not accession_counts:
        print("{} has no sublevels.".format(taxid), file=sys.stderr)
        return (0, 0)
    fasta_path = os.path.join(temp_dir, str(taxid) + ".fasta")
    taxid_path = os.path.join(temp_dir, str(taxid) + ".taxid")
    fasta_file = open(fasta_path, 'w')
    taxid_file = open(taxid_path, "w")
    fasta_records_count = get_fasta(accession_counts, length,
                                    index, fasta_file,
                                    taxid_file, window_length=window_length,
                                    temp_dir=temp_dir,
                                    thresholding=thresholding,
                                    amino_acid=amino_acid)
    fasta_file.close()
    taxid_file.close()
    permute_count = randomly_permute_fasta_taxid(fasta_path,
                                                 taxid_path,
                                                 fasta_path,
                                                 taxid_path,
                                                 split=True,
                                                 split_amount=split_amount)
    for ext, ml_path in [(".train", "train"),
                         (".validate", "validate"),
                         (".test", "test")]:
        save_dir = os.path.join(data_dir, str(taxid), ml_path)
        if not os.path.exists(save_dir):
            try:
                os.makedirs(save_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        shutil.move(fasta_path + ext, os.path.join(save_dir,
                                                   str(taxid) + ".fasta"))
        shutil.move(taxid_path + ext, os.path.join(save_dir,
                                                   str(taxid) + ".taxid"))
    if os.path.isfile(fasta_path):
        os.remove(fasta_path)
    if os.path.isfile(taxid_path):
        os.remove(taxid_path)
    return fasta_records_count, permute_count


def create_directories(data_dir):
    """
    Create directories to put the data if they do not exist.

    Parameters
    ----------
    data_dir: str
        A path to the top level of the directory to put the data.

    Returns
    -------
    None

    """
    # ranks = ["superkingdom", "kingdom", "phylum",
    #         "class", "order", "family", "genus"]
    # data_sets = ["train", "test", "validate"]
    # for rank in ranks:
    # for data_set in data_sets:
    #     path = os.path.join(data_dir, rank, data_set)
    #     if not os.path.exists(path):
    #         try:
    #             os.makedirs(path)
    #         except OSError as e:
    #             if e.errno != errno.EEXIST:
    #                 raise


def parallel_sample(taxid_list, ranks, index_dir, number, length,
                    data_dir, split_amount, processes,
                    thresholding=False, window_length=100, amino_acid=False,
                    temp_dir="/tmp"):
    """
    Get samples of data in parallel and writes them into files and a data
    directory that Plinko expects.

    Parameters
    ----------
    taxid_list: list<int>
        A list of ints representing taxonomic ids.
    ranks: dict
        The ranks object giving taxonomic ids at every rank.
    index_dir: str
        The path to the genomes index object.
    number: int
        The number of samples to take.
    length: int
        The number of bases for each sample.
    data_dir: str
        The path to the data directory where fasta files will be written.
    split_amount: str
        A comma seperated list of floating point values that represent how to
        split the training, validation, and test data sets.
    processes: int
        The number of processes to use.  This should be at or less than the
        number of physical cores that the CPU has.

    Returns
    -------
    list<(int, int)>
        A list of fasta records written by each process in the same order as
        the taxid_list.

    """
    # create_directories(data_dir)
    with Pool(processes=processes) as pool:
        process_list = []
        for taxid in taxid_list:
            sublevels = ranks[taxid]
            # my_rank = ncbi.get_rank([taxid])[taxid]
            process_list.append(pool.apply_async(get_sample,
                                                 args=(taxid,
                                                       sublevels, index_dir,
                                                       number,
                                                       length, data_dir,
                                                       split_amount,
                                                       thresholding,
                                                       amino_acid,
                                                       temp_dir)))
        output = []
        for taxid, process_desc in zip(taxid_list, process_list):
            counts = process_desc.get()
            output.append(counts)
            print("{}: {} samples drawn, {} samples written".
                  format(taxid, counts[0], counts[1]), file=sys.stderr)
        return output
    
    
class ArgClass:
    """ So that I don't have to duplicate argument info when
        the same thing is used in more than one mode."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    index = ArgClass("-i", "--index",
                     help="The location of the genomes index.",
                     default="./index.pck")
    verbose = ArgClass("-v", "--verbose", action='store_true', default=False)
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_download = subparsers.add_parser("download",
                                       help=("Download genomes from the NCBI."),
                                       formatter_class=argparse.
                                       ArgumentDefaultsHelpFormatter)
    ncbi_argument_parser(version=ncbi_version, parser=p_download)
    p_fai = subparsers.add_parser("fai",
                                  help=("Create fai files for each genome "
                                        "in the data store."),
                                  formatter_class=argparse.
                                  ArgumentDefaultsHelpFormatter)
    p_index = subparsers.add_parser("index",
                                    help=("Make a genomes and taxid index "
                                          "necessary for sampling."),
                                    formatter_class=argparse.
                                    ArgumentDefaultsHelpFormatter)
    p_tree = subparsers.add_parser("tree",
                                   help=('Construct an object '
                                   'representing a tree of '
                                   'taxonomic ranks.  '
                                   'Print a list of '
                                   'non leaf non viroid '
                                   'taxids to stdout.'),
                                   formatter_class=argparse.
                                   ArgumentDefaultsHelpFormatter)
    p_tree.add_argument("index", help=("The location of the pickled "
                                       "index file."))
    p_tree.add_argument("--output", "-o", type=str,
                        help=("The file to write the ranks tree to."),
                        default="./tree")
    p_tree.add_argument("--output_missing", "-m", type=str,
                        help=("The file location to write missing rank "
                              "information too."),
                        default="./missing.pck")
    p_tree.add_argument(*verbose.args, **verbose.kwargs)
    p_select = subparsers.add_parser("select",
                                     help=("Select which genomes to "
                                           "sample from.  The genomes"
                                           " index must be created."),
                                    formatter_class=argparse.
                                    ArgumentDefaultsHelpFormatter)
    p_sample = subparsers.add_parser("sample",
                                     help=("Create random samples of "
                                           "kmers from the genomes storage."
                                           "  Requires the taxonomic tree "
                                           "and the genomes index."),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    p_sample.add_argument(*index.args, **index.kwargs)
    p_permute = subparsers.add_parser("permute",
                                      help=("Permute and split a fasta file."),
                                      formatter_class=argparse.
                                      ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stdout.flush()
    mode = args.mode
    if mode == "download":
        ret = run_ncbi(args)
    elif mode == "fai":
        pass
    elif mode == "index":
        pass
    elif mode == "tree":
        print("Creating the taxonomic tree.", file=sys.stderr)
        (ranks, missing_ranks, unclassified,
         empties, children_removed,
         counter_dict, taxid_list) = make_tree(pickle.load(open(args.index, 'rb')))
        pickle.dump(ranks, open(args.output + '.pck', 'wb'))
        json.dump(ranks, open(args.output + '.json', 'w'), indent=4)
        pickle.dump(missing_ranks, open(args.output_missing, 'wb'))
        print("Unclassified: {}, empties: {}, children removed: {}".
              format(unclassified, empties, children_removed),
              file=sys.stderr)
        print("Counts: {}".format(counter_dict), file=sys.stderr)
        total_taxid = 0
        singleton_taxid = 0
        for taxid in counter_dict:
            total_taxid += counter_dict[taxid][0]
            singleton_taxid += counter_dict[taxid][1]
        print("Total taxid: {}, Singletons: {}, Multiclass:{}".
              format(total_taxid, singleton_taxid, total_taxid - singleton_taxid),
              file=sys.stderr)
        leaves_not_species = check_leaves_species(1, ranks)
        if leaves_not_species:
            print("There were leaves that were not species:\n{}".format(
                leaves_not_species), file=sys.stderr)
        for taxid in taxid_list:
            print(taxid)
    elif mode == "select":
        pass
    elif mode == "sample":
        pass
    elif mode == "permute":
        pass
    tock = datetime.datetime.now()
    print("Radogest {} took time: {}".format(mode, tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()