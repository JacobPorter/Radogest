#!/usr/bin/env python
"""
Radogest: random genome sampler for trees.
Works with genomic information downloaded from
NCBI (National Center for Biotechnology Information)

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

# TODO: genome holdout genome methods are not splitting the appropriate number of genomes between test and train sets.  Example: taxid 89373.
# TODO: Exclude species, levels? where there are too few genomes.

import argparse
import datetime
import json
import os
import pickle
import sys
from collections import defaultdict

from library.ncbi_genome_download.ncbi_genome_download import \
    __version__ as ncbi_version
from library.ncbi_genome_download.ncbi_genome_download.__main__ import run_ncbi
from library.ncbi_genome_download.ncbi_genome_download.core import \
    argument_parser as ncbi_argument_parser

_JSON_INDENT = 4


class ArgClass:
    """A class for reusing command line parameters."""
    def __init__(self, *args, **kwargs):
        """Initialize the field members."""
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
        formatter_class=argparse.RawTextHelpFormatter, description=__doc__)
    index = ArgClass("--index",
                     "-i",
                     type=str,
                     help="The location of the genomes index.",
                     default="./index.pck")
    verbose = ArgClass("--verbose",
                       "-v",
                       type=int,
                       help=("Controls the verbosity of the output."),
                       default=0)
    genomes = ArgClass("--genomes",
                       "-g",
                       type=str,
                       help=("The location of the genomes directory."),
                       default="./")
    tree = ArgClass("--tree",
                    "-r",
                    type=str,
                    help="The location of the taxonomic tree data structure.",
                    default="./tree.pck")
    taxid = ArgClass("--taxid",
                     "-t",
                     help=("The file location that lists taxonomic ids, "
                           "one id per file."),
                     default="./taxid_list.txt")
    leave_compressed = ArgClass("--leave_compressed",
                                "-l",
                                action='store_true',
                                help="Leave gzip files compressed.",
                                default=False)
    window_length = ArgClass("--window_length",
                             "-w",
                             type=int,
                             help=("The window length to use when using "
                                   "thresholding or chopping."),
                             default=50)
    kmer_size = ArgClass("--kmer_size",
                         "-k",
                         type=int,
                         help=("The length in base pairs i "
                               "of the kmers to take."),
                         default=100)
    split = ArgClass("--split",
                     help=("If set, splits the data into "
                           "training, test sets, etc."),
                     action='store_true',
                     default=False)
    split_amount = ArgClass("--split_amount",
                            "-s",
                            type=str,
                            help=("A comma separated list of "
                                  "three percentages.  "
                                  "The percentages to split the data into "
                                  "for training, validation, and test sets."),
                            default="0.8,0.10,0.10")
    include_wild = ArgClass("--include_wild",
                            "-x",
                            help=("Include wild card characters "
                                  "in the samples.  "
                                  "When this is not used, samples with "
                                  "wild card characters will be discarded."),
                            action='store_true',
                            default=False)
    prob = ArgClass("--prob",
                    "-b",
                    type=float,
                    help=("The probability that a sequence will be "
                          "converted to the reverse complement."),
                    default=0.5)
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_download = subparsers.add_parser(
        "download",
        help=("Download genomes from "
              "the NCBI. "
              "If used by Radogest, "
              "each data type "
              "(nucleotide, coding domain, "
              "amino acid) should be "
              "downloaded into its own "
              "directory."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ncbi_argument_parser(version=ncbi_version, parser=p_download)
    p_fai = subparsers.add_parser(
        "faidx",
        help=("Create fai files for each genome "
              "in the data store."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_fai.add_argument(*genomes.args, **genomes.kwargs)
    p_fai.add_argument(*leave_compressed.args, **leave_compressed.kwargs)
    p_fai.add_argument(*verbose.args, **verbose.kwargs)
    p_fai.add_argument("--processes",
                       "-p",
                       type=int,
                       default=2,
                       help=("The number of processes to use."))
    p_index = subparsers.add_parser(
        "index",
        help=("Make a genomes and taxid index "
              "necessary for sampling."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_index.add_argument("json",
                         type=str,
                         help=("The path to a directory of json files "
                               "or to a single json file."))
    p_index.add_argument(*index.args, **index.kwargs)
    p_index.add_argument(*genomes.args, **genomes.kwargs)
    p_index.add_argument(*verbose.args, **verbose.kwargs)
    p_tree = subparsers.add_parser(
        "tree",
        help=('Construct an object '
              'representing a tree of '
              'taxonomic ranks.  '
              'Print a list of '
              'non leaf non viroid '
              'taxids to stdout.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_tree.add_argument(*index.args, **index.kwargs)
    p_tree.add_argument(*tree.args, **tree.kwargs)
    p_tree.add_argument(*taxid.args, **taxid.kwargs)
    p_tree.add_argument(*verbose.args, **verbose.kwargs)
    p_select = subparsers.add_parser(
        "select",
        help=("Select which genomes to "
              "sample from.  The genomes"
              " index must be created."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_select.add_argument(*index.args, **index.kwargs)
    p_select.add_argument(*tree.args, **tree.kwargs)
    p_select.add_argument("taxid",
                          type=int,
                          help=('The taxonomic id for '
                                'the root.'))
    p_select.add_argument(
        "--strategy",
        '-s',
        action='store',
        choices=['PR', 'QST', 'QSL', 'GHST', 'GHSL', 'GHGT', 'GHGL', 'AG'],
        help=('Choose the genome selection strategy.  '
              'The choices are: ProportionalRandom (PR), '
              'QualitySortingTree (QST), '
              'QualitySortingLeaf (QSL),'
              'GenomeHoldoutSpeciesTree (GHST), '
              'GenomeHoldoutSpeciesLeaf (GHSL), '
              'GenomeHoldoutGenomeTree (GHGT), '
              'GenomeHoldoutGenomeLeaf (GHGL), '
              'AllGenomes (AG)'),
        default='PR')
    p_select.add_argument('--select_amount',
                          '-n',
                          nargs='+',
                          type=int,
                          help=('The number of genomes '
                                'to select at each level. '
                                'Does not apply to AllGenomes.  '
                                'Genome holdout requires '
                                'two amounts: '
                                'the first one for training '
                                'data, and the second one for '
                                'testing data.'),
                          default=[10])
    p_select.add_argument('--random',
                          '-d',
                          action='store_true',
                          help="Select genomes at random rather than sort "
                          "them by quality when selecting genomes for "
                          "genome holdout strategies.  This option only "
                          "applies to genome holdout strategies.",
                          default=False)
    p_select.add_argument('--output',
                          '-o',
                          type=str,
                          help=('The location to store the index with '
                                'the down selected genomes.'),
                          default='./index_down.pck')
    p_select.add_argument(*verbose.args, **verbose.kwargs)
    p_sample = subparsers.add_parser(
        "sample",
        help=("Create random samples of "
              "kmers from the genomes storage."
              "  Requires the taxonomic tree "
              "and the genomes index."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_sample.add_argument("taxid",
                          help=("The file location that lists taxonomic ids, "
                                "one id per file.  "
                                "Or, a taxonomic id when sampling "
                                "from one taxonomic id."))
    p_sample.add_argument(*index.args, **index.kwargs)
    p_sample.add_argument(*tree.args, **tree.kwargs)
    p_sample.add_argument(*genomes.args, **genomes.kwargs)
    p_sample.add_argument(*kmer_size.args, **kmer_size.kwargs)
    p_sample.add_argument("--number",
                          "-n",
                          type=int,
                          help=("The number of samples to take."),
                          default=3200000)
    p_sample.add_argument("--range",
                          "-e",
                          type=int,
                          nargs=2,
                          help=("The subselection of taxonomic ids from the "
                                "taxid file to use. This is useful for "
                                "distributing the work across mutliple "
                                "compute nodes using the same taxid file.  "
                                "Use -1 to denote the end of the file. "),
                          default=[0, -1])
    p_sample.add_argument("--thresholding",
                          "-m",
                          help=("If the number of kmers requested is larger "
                                "than the number of kmers in the fasta "
                                "file, chop up the fasta file rather than "
                                "randomly sample from it."),
                          action='store_true',
                          default=False)
    p_sample.add_argument('--thresholds',
                          '-l',
                          nargs='+',
                          type=int,
                          help=('The total amount of genomic information '
                                'to get in kilobases.  '
                                'This option only applies to '
                                'genome holdout strategies when '
                                'multiple threholds are specified.  '
                                'A single value can be used for '
                                'other strategies.  '
                                'If None, then all of the genomes '
                                'in the set will be used.  '
                                'This option has nothing to do with '
                                'the thresholding option.'),
                          default=None)
    p_sample.add_argument("--chop",
                          "-c",
                          action="store_true",
                          default=False,
                          help=("Chop up the set of genomes in "
                                "the sample set.  "
                                "This will NOT randomly sample kmers.  "
                                "The kmers in the sample set "
                                "may not be balanced."))
    p_sample.add_argument("--save_genomes",
                          action="store_true",
                          default=False,
                          help="Save genomes in a separate folder when "
                          "using genome holdout strategies.")
    p_sample.add_argument(*window_length.args, **window_length.kwargs)
    p_sample.add_argument(*split.args, **split.kwargs)
    p_sample.add_argument(*split_amount.args, **split_amount.kwargs)
    p_sample.add_argument(*include_wild.args, **include_wild.kwargs)
    p_sample.add_argument(*prob.args, **prob.kwargs)
    p_sample.add_argument("--amino_acid",
                          "-a",
                          action="store_true",
                          help=("Turn this switch on when sampling amino "
                                "acids."),
                          default=False)
    p_sample.add_argument("--data_dir",
                          "-d",
                          type=str,
                          help=("The directory to store data."),
                          default="./Data/")
    p_sample.add_argument("--temp_dir",
                          "-u",
                          type=str,
                          help=("The temporary directory to store "
                                "intermediate files such as BED files and "
                                "fasta files or sampling from individual "
                                " genomes. "),
                          default='/tmp/')
    p_sample.add_argument("--processes",
                          "-p",
                          type=int,
                          help=("The number of processes to use."),
                          default=16)
    p_sample.add_argument(*verbose.args, **verbose.kwargs)
    input_fasta = ArgClass("input_fasta",
                           type=str,
                           help=("The input fasta file."))
    input_taxid = ArgClass("input_taxid",
                           type=str,
                           help=("The input taxid file."))
    output_fasta = ArgClass("output_fasta",
                            type=str,
                            help=("The output fasta file."))
    output_taxid = ArgClass("output_taxid",
                            type=str,
                            help=("The output taxid file."))
    p_permute = subparsers.add_parser(
        "util_permute",
        help=("Permute a fasta file "
              "and a taxid file."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_permute.add_argument(*input_fasta.args, **input_fasta.kwargs)
    p_permute.add_argument(*input_taxid.args, **input_taxid.kwargs)
    p_permute.add_argument(*output_fasta.args, **output_fasta.kwargs)
    p_permute.add_argument(*output_taxid.args, **output_taxid.kwargs)
    p_permute.add_argument(*split.args, **split.kwargs)
    p_permute.add_argument(*split_amount.args, **split_amount.kwargs)
    p_chop = subparsers.add_parser(
        "util_chop",
        help=("Chop a complete genome "
              "into kmers and write it out"),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_chop.add_argument("taxid",
                        type=int,
                        help=("The taxonomic id to include in the "
                              "fasta record id."))
    p_chop.add_argument(*output_fasta.args, **output_fasta.kwargs)
    p_chop.add_argument("accessions",
                        type=str,
                        nargs='+',
                        help=("The accession id of the genome."))
    p_chop.add_argument(*index.args, **index.kwargs)
    p_chop.add_argument(*genomes.args, **genomes.kwargs)
    p_chop.add_argument(*kmer_size.args, **kmer_size.kwargs)
    p_chop.add_argument(*include_wild.args, **include_wild.kwargs)
    p_chop.add_argument(*window_length.args, **window_length.kwargs)
    p_rc = subparsers.add_parser(
        "util_rc",
        help=("Randomly take the reverse "
              "complement of DNA sequences "
              "in a fasta file."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_rc.add_argument(*input_fasta.args, **input_fasta.kwargs)
    p_rc.add_argument(*output_fasta.args, **output_fasta.kwargs)
    p_rc.add_argument("--exclude_wild",
                      "-x",
                      action="store_true",
                      help=("Remove DNA sequences with wildcard characters."),
                      default=False)
    p_rc.add_argument(*prob.args, **prob.kwargs)
    p_rc.add_argument(*verbose.args, **verbose.kwargs)
    p_subtree = subparsers.add_parser(
        "util_subtree",
        help=("Get a list of all taxonomic "
              "ids underneath and including "
              "the given taxonomic id."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_subtree.add_argument("taxid",
                           help=("The taxid for the root of the subtree."))
    p_subtree.add_argument("--species",
                           "-s",
                           action="store_true",
                           help=("Include species taxonomic ids."),
                           default=False)
    p_subtree.add_argument(*tree.args, **tree.kwargs)
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
        print("The download was {} with code: {}".format(success_string, ret),
              file=sys.stderr)
    elif mode == "faidx":
        from library.faidx import make_fai
        sys.stderr.write("Starting to go down {} to create fai files. "
                         " leave_compressed: {}, "
                         "verbose: {}.\n".format(args.genomes,
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
        index["select"] = {"strategy": "INIT", "select_amount": None}
        try:
            write_ds(index, args.index)
        except IOError:
            parser.error(index_write_error)
        ncbi.update_taxonomy_database()
        if os.path.isfile(args.json):
            count = create_initial_index(index, json.load(open(args.json)),
                                         args.verbose)
        elif os.path.isdir(args.json):
            onlyfiles = [
                f for f in os.listdir(args.json)
                if os.path.isfile(os.path.join(args.json, f))
            ]
            onlyfiles = [f for f in onlyfiles if f.endswith('.json')]
            count = 0
            for f_f in onlyfiles:
                sys.stderr.write("Processing file {}.\n".format(f_f))
                my_path = os.path.join(args.json, f_f)
                count += create_initial_index(index, json.load(open(my_path)),
                                              args.verbose)
        else:
            sys.stderr.write("The json parameter was not a "
                             "directory and was not "
                             "a file.  Please to check.\n")
            parser.print_help()
        index, count_update = update_index_root(index, args.genomes,
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
        (tree, unclassified, empties, children_removed, counter_dict,
         taxid_list) = make_tree(read_ds(args.index))
        write_ds(tree, args.tree)
        print("Unclassified: {}, empties: {}, children removed: {}".format(
            unclassified, empties, children_removed),
              file=sys.stderr)
        print("Counts: {}".format(counter_dict), file=sys.stderr)
        total_taxid = 0
        singleton_taxid = 0
        for taxid in counter_dict:
            total_taxid += counter_dict[taxid][0]
            singleton_taxid += counter_dict[taxid][1]
        print("Total taxid: {}, Singletons: {}, Multiclass:{}".format(
            total_taxid, singleton_taxid, total_taxid - singleton_taxid),
              file=sys.stderr)
        leaves_not_species = check_leaves_species(1, tree)
        if leaves_not_species:
            print("There were leaves that were not species:\n{}".format(
                leaves_not_species),
                  file=sys.stderr)
        with open(args.taxid, "w") as taxid_file:
            for taxid in taxid_list:
                print(taxid, file=taxid_file)
    elif mode == "select":
        from library.genome_selection.strategy import ProportionalRandom
        from library.genome_selection.strategy import QualitySortingTree
        from library.genome_selection.strategy import QualitySortingLeaf
        from library.genome_selection.strategy import AllGenomes
        from library.genome_selection.strategy import GHSpeciesLeaf
        from library.genome_selection.strategy import GHSpeciesTree
        from library.genome_selection.strategy import GHGenomeLeaf
        from library.genome_selection.strategy import GHGenomeTree
        from library.genome_selection.strategy import EXCLUDED_GENOMES
        from library.genome_selection.traversal import StrategyNotFound
        from library.genome_selection.traversal import TaxTreeTraversal
        strategy_string = args.strategy.upper()
        index = read_ds(args.index)
        tree = read_ds(args.tree)
        select_amount = args.select_amount
        if not min(list(map(lambda x: x > 0, select_amount))):
            parser.error('The sample amount needs to be a positive integer.')
        if strategy_string == 'PR':
            strategy = ProportionalRandom(index, select_amount[0])
        elif strategy_string == 'QST':
            strategy = QualitySortingTree(index, select_amount[0])
        elif strategy_string == 'QSL':
            strategy = QualitySortingLeaf(index, select_amount[0])
        elif strategy_string.startswith('GH'):
            if len(select_amount) != 2:
                p_select.error("There must be two and only two sample amounts "
                               "when using the genome holdout strategy: "
                               "The first for the training set, "
                               "the second for the testing set.")
            warning = ("Species holdout strategies may not make sense "
                       "for selecting genera.")
            if strategy_string == 'GHSL':
                strategy = GHSpeciesLeaf(index,
                                         select_amount,
                                         random=args.random)
                print(warning, file=sys.stderr)
            elif strategy_string == 'GHST':
                strategy = GHSpeciesTree(index,
                                         select_amount,
                                         random=args.random)
                print(warning, file=sys.stderr)
            elif strategy_string == 'GHGL':
                strategy = GHGenomeLeaf(index,
                                        select_amount,
                                        random=args.random)
            elif strategy_string == 'GHGT':
                strategy = GHGenomeTree(index,
                                        select_amount,
                                        random=args.random)
            else:
                raise StrategyNotFound()
        elif strategy_string == 'AG':
            strategy = AllGenomes(index)
        else:
            raise StrategyNotFound()
        index["select"] = {
            "strategy": strategy_string,
            "select_amount": select_amount
        }
        traversal = TaxTreeTraversal(tree, strategy)
        levels_visited = traversal.select_genomes(args.taxid)
        if args.verbose >= 1:
            for accession in EXCLUDED_GENOMES:
                print("WARNING: {} excluded because {}.".format(
                    accession, EXCLUDED_GENOMES[accession]),
                      file=sys.stderr)
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
            taxid_list = [
                int(taxid.strip()) for taxid in open(args.taxid, "r")
            ]
        begin, end = args.range
        if end == -1:
            end = len(taxid_list)
        taxid_list = taxid_list[begin:end]
        if (not os.path.isdir(args.temp_dir)
                or not os.path.exists(args.temp_dir)):
            parser.error("The temporary directory could not be found "
                         "or does not exist.")
        if args.prob < 0.0 or args.prob > 1.0:
            parser.error("The reverse complement probability {} "
                         "is not a valid probability.".format(args.prob))
        parallel_sample(taxid_list, args.genomes, tree, args.index,
                        args.number, args.kmer_size, args.data_dir, args.split,
                        args.split_amount, args.processes, args.include_wild,
                        args.prob, args.thresholding, args.chop,
                        args.save_genomes, args.window_length, args.amino_acid,
                        args.thresholds, args.temp_dir, args.verbose)
    elif mode == "util_permute":  # permute and split
        from library.permute import randomly_permute_fasta_taxid
        permute_count = randomly_permute_fasta_taxid(
            args.input_fasta,
            args.input_taxid,
            args.output_fasta,
            args.output_taxid,
            split=args.split,
            split_amount=args.split_amount)
        print("There were {} records written.".format(permute_count),
              file=sys.stderr)
    elif mode == "util_chop":
        from library.chop import chop_genomes
        index = read_ds(args.index)
        locations = [
            args.genomes + index['genomes'][accession]['location']
            for accession in args.accessions
        ]
        number_written = chop_genomes(args.accessions, args.kmer_size,
                                      locations, args.taxid, args.output_fasta,
                                      None, args.include_wild,
                                      args.window_length)
        print("There were {} records written.".format(number_written),
              file=sys.stderr)
    elif mode == "util_rc":  # reverse complement
        from library.sample import get_rc_fasta
        read_counter, write_counter = get_rc_fasta(args.input_fasta,
                                                   args.output_fasta,
                                                   prob=args.prob,
                                                   remove=args.exclude_wild,
                                                   verbose=args.verbose)
        print("There were {} records read and {} records written.".format(
            read_counter, write_counter),
              file=sys.stderr)
    elif mode == "util_subtree":
        from library.taxid_subtree import get_subtree
        root = int(args.taxid)
        tree = pickle.load(open(args.tree, "rb"))
        taxid_list = []
        get_subtree(tree, root, taxid_list, args.species)
        for t in taxid_list:
            print(t, file=sys.stdout)
    else:
        parser.print_usage()
        print("There was no command specified.", file=sys.stderr)
        sys.exit(1)
    tock = datetime.datetime.now()
    print("Radogest {} took time: {}".format(mode, tock - tick),
          file=sys.stderr)


if __name__ == "__main__":
    main()
