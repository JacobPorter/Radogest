#!/usr/bin/env python
"""
Radogest: random genome sampler for trees.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


# TODO: Get a set of maximally distant genomes as a sampling strategy.  Per Andrew.  Handle case where a taxid has only one genome?  See 8892 for vertebrate refseq data.
# TODO: Provide annotation mapping (Coding domain, etc. info from gbff file)?
# TODO: Allow for all file types to be downloaded into the same directory?  Need to include file type information in the index.
# TODO: Add better parallelism to sampling, index creation.  PySpark?  Process pool?
# TODO: Add utility commands?  permute, split, chop, reverse_complement.
# TODO: When finished with code, check and update comments and README documentation.
# TODO: Write and submit a paper.

import argparse
import datetime
import sys
import json
import pickle
import os

from collections import defaultdict

from library.ncbi_genome_download.ncbi_genome_download.core import argument_parser as ncbi_argument_parser
from library.ncbi_genome_download.ncbi_genome_download import __version__ as ncbi_version
from library.ncbi_genome_download.ncbi_genome_download.__main__ import run_ncbi

_JSON_INDENT = 4


class ArgClass:

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


def write_ds(ds, path):
    """Write the data structure to the path."""
    with open(path, 'wb') as ds_file:
        ret = pickle.dump(ds, ds_file)
    return ret


def read_ds(path):
    """Read the data structure from the path."""
    with open(path, 'rb') as ds_file:
        ds = pickle.load(ds_file)
    return ds


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    print("Radogest was started at {}.".format(tick), file=sys.stderr)
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    index = ArgClass("--index", "-i", type=str,
                     help="The location of the genomes index.",
                     default="./index.pck")
    verbose = ArgClass("--verbose", "-v", type=int,
                       help=("Controls the verbosity of the output."),
                       default=0)
    genomes = ArgClass("--genomes", "-g", type=str,
                       help=("The location of the genomes directory."),
                       default="./")
    tree = ArgClass("--tree", "-r", type=str,
                    help="The location of the taxonomic tree data structure.",
                    default="./tree.pck")
    taxid = ArgClass("--taxid", "-t",
                     help=("The file location that lists taxonomic ids."
                           "One id per file."),
                     default="./taxid_list.txt")
    leave_compressed = ArgClass("--leave_compressed", "-l",
                                action='store_true',
                                help="Leave gzip files compressed.",
                                default=False)
#     fasta = ArgClass("fasta", type=str,
#                      help="The location of a fasta file.")
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_download = subparsers.add_parser("download",
                                       help=("Download genomes from "
                                             "the NCBI."),
                                       formatter_class=argparse.
                                       ArgumentDefaultsHelpFormatter)
    ncbi_argument_parser(version=ncbi_version, parser=p_download)
    p_fai = subparsers.add_parser("faidx",
                                  help=("Create fai files for each genome "
                                        "in the data store."),
                                  formatter_class=argparse.
                                  ArgumentDefaultsHelpFormatter)
    p_fai.add_argument(*genomes.args, **genomes.kwargs)
    p_fai.add_argument(*leave_compressed.args, **leave_compressed.kwargs)
    p_fai.add_argument(*verbose.args, **verbose.kwargs)
    p_fai.add_argument("--processes", "-p", type=int, default=2,
                       help=("The number of processes to use."))
    p_index = subparsers.add_parser("index",
                                    help=("Make a genomes and taxid index "
                                          "necessary for sampling."),
                                    formatter_class=argparse.
                                    ArgumentDefaultsHelpFormatter)
    p_index.add_argument("json", type=str,
                         help=("The path to a directory of json files "
                               "or to a single json file."))
    p_index.add_argument(*index.args, **index.kwargs)
    p_index.add_argument(*genomes.args, **genomes.kwargs)
    p_index.add_argument(*verbose.args, **verbose.kwargs)
    p_tree = subparsers.add_parser("tree",
                                   help=('Construct an object '
                                         'representing a tree of '
                                         'taxonomic ranks.  '
                                         'Print a list of '
                                         'non leaf non viroid '
                                         'taxids to stdout.'),
                                   formatter_class=argparse.
                                   ArgumentDefaultsHelpFormatter)
    p_tree.add_argument(*index.args, **index.kwargs)
    p_tree.add_argument(*tree.args, **tree.kwargs)
    p_tree.add_argument(*taxid.args, **taxid.kwargs)
    p_tree.add_argument(*verbose.args, **verbose.kwargs)
    p_select = subparsers.add_parser("select",
                                     help=("Select which genomes to "
                                           "sample from.  The genomes"
                                           " index must be created."),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    p_select.add_argument(*index.args, **index.kwargs)
    p_select.add_argument(*tree.args, **tree.kwargs)
    p_select.add_argument("taxid", type=int, help=('The taxonomic id for '
                                                   'the root.'))
    p_select.add_argument("--strategy", '-s', action='store',
                          choices=['PR', 'QST', 'QSL', 'GH', 'AG'],
                          help=('Choose the genome selection strategy.  '
                                'The choices are: ProportionalRandom (PR), '
                                'QualitySortingTree (QST), QualitySortingLeaf '
                                '(QSL), GenomeHoldout (GH), AllGenomes (AG)'),
                          default='PR')
    p_select.add_argument('--sample_amount', '-n', nargs='+', type=int,
                          help=('The number of genomes '
                                'to sample at each level. '
                                'Does not apply to AllGenomes.'),
                          default=[10])
    p_select.add_argument('--output', '-o', type=str,
                          help=('The location to store the index with '
                                'the down selected genomes.'),
                          default='./index_down.pck')
    p_sample = subparsers.add_parser("sample",
                                     help=("Create random samples of "
                                           "kmers from the genomes storage."
                                           "  Requires the taxonomic tree "
                                           "and the genomes index."),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    p_sample.add_argument("taxid",
                          help=("The file location that lists taxonomic ids."
                                "One id per file.  "
                                "Or, a taxonomic id when sampling "
                                "from one taxonomic id."))
    p_sample.add_argument(*index.args, **index.kwargs)
    p_sample.add_argument(*tree.args, **tree.kwargs)
    p_sample.add_argument(*genomes.args, **genomes.kwargs)
    p_sample.add_argument("--kmer_size", "-k", type=int,
                          help=("The length in base pairs i "
                                "of the kmers to take."),
                          default=100)
    p_sample.add_argument("--number", "-n", type=int,
                          help=("The number of samples to take."),
                          default=3200000)
    p_sample.add_argument("--range", "-e", type=int, nargs=2,
                          help=("The subselection of taxonomic ids from the "
                                "taxid file to use. This is useful for "
                                "distributing the work across mutliple "
                                "compute nodes using the same taxid file.  "
                                "Use -1 to denote the end of the file. "),
                          default=[0, -1])
    p_sample.add_argument("--include_wild", "-x",
                          help=("Include wild card characters "
                                "in the samples.  "
                                "When this is not used, samples with "
                                "wild card characters will be discarded."),
                          action='store_true', default=False)
    p_sample.add_argument("--thresholding", "-m",
                          help=("If the number of kmers requested is larger "
                                "than the number of kmers in the fasta "
                                "file, chop up the fasta file rather than "
                                "randomly sample from it."),
                          action='store_true', default=False)
    p_sample.add_argument("--chop", "-c", action="store_true", default=False,
                          help=("Chop up the set of genomes in "
                                "the sample set.  "
                                "This will NOT randomly sample kmers.  "
                                "The kmers in the sample set "
                                "may not be balanced."))
    p_sample.add_argument("--window_length", "-w", type=int,
                          help=("The window length to use when using "
                                "thresholding."),
                          default=50)
    p_sample.add_argument("--split",
                          help=("If set, splits the data into "
                                "training, test sets, etc."),
                          action='store_true', default=False)
    p_sample.add_argument("--split_amount", "-s", type=str,
                          help=("A comma seperated list of three percentages."
                                "The percentages to split the data into "
                                "for training, validation, and test sets."),
                          default="0.8,0.10,0.10")
    p_sample.add_argument("--prob", "-b", type=float,
                        help=("The probability that a sequence will be "
                              "converted to the reverse complement."),
                        default=0.5)
    p_sample.add_argument("--amino_acid", "-a", action="store_true",
                          help=("Turn this switch on when sampling amino "
                                "acids."),
                          default=False)
    p_sample.add_argument("--data_dir", "-d", type=str,
                          help=("The directory to store data."),
                          default="./Data/")
    p_sample.add_argument("--temp_dir", "-u", type=str,
                          help=("The temporary directory to store "
                                "intermediate files such as BED files and "
                                "fasta files or sampling from individual "
                                " genomes. "),
                          default='/localscratch/')
    p_sample.add_argument("--processes", "-p", type=int,
                          help=("The number of processes to use."),
                          default=16)
    p_sample.add_argument(*verbose.args, **verbose.kwargs)
#     p_permute = subparsers.add_parser("permute",
#                                       help=("Permutes a fasta file "
#                                             "and writes it to stdout"),
#                                       formatter_class=argparse.
#                                       ArgumentDefaultsHelpFormatter)
#     p_permute.add_argument(*fasta.args, **fasta.kwargs)
#     p_chop = subparsers.add_parser("chop",
#                                    help=("Chop a complete genome "
#                                          "into kmers and writes it "
#                                          "to stdout."),
#                                    formatter_class=argparse.
#                                    ArgumentDefaultsHelpFormatter)
#     p_chop.add_argument(*fasta.args, **fasta.kwargs)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    mode = args.mode
    index_write_error = ("Something went wrong writing the index.  "
                         "Check that the path is correct and writable.")
    if mode == "download":
        ret = run_ncbi(args)
        if ret == 0:
            success_string = "SUCCESSFUL"
        else:
            success_string = "FAILED"
        print("The download was {} with code: {}".format(
            success_string, ret), file=sys.stderr)
    elif mode == "faidx":
        from library.faidx import make_fai
        sys.stderr.write("Starting to go down {} to create fai files. "
                         " leave_compressed: {}, "
                         "verbose: {}.\n"
                         .format(args.genomes,
                                 args.leave_compressed,
                                 args.verbose))
        count = make_fai(args.genomes,
                         processes=args.processes,
                         leave_compressed=args.leave_compressed,
                         verbose=args.verbose)
        print(count, file=sys.stderr)
    elif mode == "index":
        from library.index import create_initial_index, update_index_root
        from library.index import ncbi
        index = {'taxids': defaultdict(dict), 'genomes': defaultdict(dict)}
        index["select"] = {"strategy": "INIT",
                           "sample_amount": None}
        try:
            write_ds(index, args.index)
        except IOError:
            parser.error(index_write_error)
        ncbi.update_taxonomy_database()
        if os.path.isfile(args.json):
            count = create_initial_index(index,
                                         json.load(open(args.json)),
                                         args.verbose)
        elif os.path.isdir(args.json):
            onlyfiles = [f for f in os.listdir(args.json)
                         if os.path.isfile(os.path.join(
                             args.json, f))]
            onlyfiles = [f for f in onlyfiles if f.endswith('.json')]
            count = 0
            for f_f in onlyfiles:
                sys.stderr.write("Processing file {}.\n".format(f_f))
                my_path = os.path.join(args.json, f_f)
                count += create_initial_index(index,
                                              json.load(open(my_path)),
                                              args.verbose)
        else:
            sys.stderr.write("The json parameter was not a "
                             "directory and was not "
                             "a file.  Please to check.\n")
            parser.print_help()
        index, count_update = update_index_root(index,
                                                args.genomes,
                                                args.verbose)
        print("The initial index counts were: {}".format(count),
              file=sys.stderr)
        print("The index update counts were: {}".format(count_update),
              file=sys.stderr)
        try:
            write_ds(index, args.index)
        except IOError:
            parser.error(index_write_error)
    elif mode == "tree":
        from library.tree import make_tree, check_leaves_species
        print("Creating the taxonomic tree.", file=sys.stderr)
        (tree, unclassified,
         empties, children_removed,
         counter_dict, taxid_list) = make_tree(
             read_ds(args.index))
        write_ds(tree, args.tree)
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
              format(total_taxid, singleton_taxid,
                     total_taxid - singleton_taxid),
              file=sys.stderr)
        leaves_not_species = check_leaves_species(1, tree)
        if leaves_not_species:
            print("There were leaves that were not species:\n{}".format(
                leaves_not_species), file=sys.stderr)
        with open(args.taxid, "w") as taxid_file:
            for taxid in taxid_list:
                print(taxid, file=taxid_file)
    elif mode == "select":
        from library.genome_selection.strategy import ProportionalRandom
        from library.genome_selection.strategy import QualitySortingTree
        from library.genome_selection.strategy import QualitySortingLeaf
        from library.genome_selection.strategy import AllGenomes
        from library.genome_selection.strategy import GenomeHoldout
        from library.genome_selection.strategy import EXCLUDED_GENOMES
        from library.genome_selection.traversal import StrategyNotFound
        from library.genome_selection.traversal import TaxTreeTraversal
        strategy_string = args.strategy.upper()
        index = read_ds(args.index)
        tree = read_ds(args.tree)
        sample_amount = args.sample_amount
        if not min(list(map(lambda x: x > 0, sample_amount))):
            parser.error('The sample amount needs to be a positive integer.')
        if strategy_string == 'PR':
            strategy = ProportionalRandom(index, sample_amount[0])
        elif strategy_string == 'QST':
            strategy = QualitySortingTree(index, sample_amount[0])
        elif strategy_string == 'QSL':
            strategy = QualitySortingLeaf(index, sample_amount[0])
        elif strategy_string == 'GH':
            strategy = GenomeHoldout(index, sample_amount)
        elif strategy_string == 'AG':
            strategy = AllGenomes(index)
        else:
            raise StrategyNotFound()
        index["select"] = {"strategy": strategy_string,
                           "sample_amount": sample_amount}
        traversal = TaxTreeTraversal(tree, strategy)
        levels_visited = traversal.select_genomes(args.taxid)
        for accession in EXCLUDED_GENOMES:
            print("WARNING: {} excluded because {}.".format(
                accession, EXCLUDED_GENOMES[accession]), file=sys.stderr)
        print("Levels visited: {}".format(levels_visited), file=sys.stderr)
        print("Creating a pickled index at {}.".format(args.output),
              file=sys.stderr)
        write_ds(index, args.output)
    elif mode == "sample":
        from library.sample import parallel_sample
        if not os.path.isfile(args.tree) or not os.path.isfile(args.index):
            parser.error("The tree object or the index object could not be "
                         "found.  Check the paths.\n")
        tree = read_ds(args.tree)
        try:
            taxid_list = [int(args.taxid)]
        except ValueError:
            taxid_list = [int(taxid.strip()) for taxid in open(args.taxid, "r")]
        begin, end = args.range
        if end == -1:
            end = len(taxid_list)
        taxid_list = taxid_list[begin:end]
        if (not os.path.isdir(args.temp_dir) or not
                os.path.exists(args.temp_dir)):
            parser.error("The temporary directory could not be found "
                         "or does not exist.")
        if args.prob < 0.0 or args.prob > 1.0:
            parser.error("The reverse complement probability {} "
                         "is not a valid probability.".format(args.prob))
        parallel_sample(taxid_list, args.genomes,
                        tree, args.index, args.number,
                        args.kmer_size, args.data_dir,
                        args.split, args.split_amount,
                        args.processes,
                        args.include_wild,
                        args.prob, 
                        args.thresholding,
                        args.chop,
                        args.window_length, args.amino_acid,
                        args.temp_dir,
                        args.verbose)
#     elif mode == "permute":
#         pass
#     elif mode == "chop":
#         pass
    else:
        parser.print_usage()
        print("There was no command specified.", file=sys.stderr)
        sys.exit(1)
    tock = datetime.datetime.now()
    print("Radogest {} took time: {}".format(mode, tock - tick),
          file=sys.stderr)


if __name__ == "__main__":
    main()
