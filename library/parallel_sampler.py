#!/usr/bin/env python
"""
Runs multiple taxonomic id sampling instances in parallel.
Creates training, validation,
and testing data and puts them in directories that Plinko expects.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import os
import errno
import pickle
import shutil
import datetime
import sys
# import sqlite3
from multiprocessing import Pool
# from ete3 import NCBITaxa
from random_genome_sampler import uniform_samples_at_rank, get_fasta
from permute_split_fasta_taxid import randomly_permute_fasta_taxid

# ncbi = NCBITaxa()


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


def main():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description=('Creates random samples '
                                                  'from a list of taxonomic '
                                                  'ids for '
                                                  'Plinko training.  Random '
                                                  'samples are taken in '
                                                  'parallel.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("taxid", type=str,
                        help=("A file location where taxids are on each line"),
                        default="./taxids.txt")
    parser.add_argument("--tree", "-t", type=str,
                        help=("The location of the tree pickled "
                              "object."),
                        default="./tree.pck")
    parser.add_argument("--index", "-i", type=str,
                        help=("The location of the genomes index object."),
                        default="./genomes.withcounts.selected.index.pck")
    parser.add_argument("--length", "-l", type=int,
                        help=("The length in base pairs of samples to take."),
                        default=100)
    parser.add_argument("--window_length", "-w", type=int,
                        help=("The window length to use when using "
                              "thresholding."),
                        default=50)
    parser.add_argument("--number", "-n", type=int,
                        help=("The number of samples to take."),
                        default=3200000)
    parser.add_argument("--range", "-r", type=int, nargs=2,
                        help=("The subselection of taxonomic ids from the "
                              "taxid file to use. This is useful for "
                              "distributing the work across mutliple "
                              "compute nodes using the same taxid file.  "
                              "Use -1 to denote the end of the file. "),
                        default=[0, -1])
    parser.add_argument("--thresholding", "-g",
                        action='store_true', default=False)
    parser.add_argument("--split_amount", "-s", type=str,
                        help=("A comma seperated list of three percentages."
                              "The percentages to split the data into "
                              "for training, validation, and test sets."),
                        default="0.8,0.10,0.10")
    parser.add_argument("--amino_acid", "-a", action="store_true",
                        help=("Turn this switch on when sampling amino "
                              "acids."),
                        default=False)
    parser.add_argument("--data_dir", "-d", type=str,
                        help=("The directory to store data."),
                        default="./Data/")
    parser.add_argument("--temp_dir", "-u", type=str,
                        help=("The temporary directory to store "
                              "intermediate files such as BED files and "
                              "fasta files or sampling from individual "
                              " genomes. "),
                        default='/localscratch')
    parser.add_argument("--processes", "-p", type=int,
                        help=("The number of processes to use."),
                        default=16)
    args = parser.parse_args()
    now = datetime.datetime.now()
    sys.stderr.write("The arguments for parallel sampling "
                     "are:\n{}\n".format(str(args)))
    if not os.path.isfile(args.tree) or not os.path.isfile(args.index):
        parser.error("The tree object or the index object could not be "
                     "found.  Check the paths.\n")
    ranks = pickle.load(open(args.tree, "rb"))
    taxid_list = [int(taxid.strip()) for taxid in open(args.taxid, "r")]
    begin, end = args.range
    if end == -1:
        end = len(taxid_list)
    taxid_list = taxid_list[begin:end]
    if not os.path.isdir(args.temp_dir) or not os.path.exists(args.temp_dir):
        parser.error("The temporary directory could not be found or does not "
                     "exist.")
    parallel_sample(taxid_list,
                    ranks, args.index, args.number, args.length, args.data_dir,
                    args.split_amount, args.processes, args.thresholding,
                    args.window_length, args.amino_acid,
                    args.temp_dir)
    later = datetime.datetime.now()
    sys.stderr.write("The process took time: {}\n".format(later - now))


if __name__ == "__main__":
    main()
