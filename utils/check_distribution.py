#!/usr/bin/env python
"""
Print out the number of each taxid in a
Radogest generated fasta file.

:Authors:
    Jacob Porter <jsporter@virginia.edu>

"""
import argparse
import datetime
import json
import pickle
import sys
from collections import defaultdict

from ete3 import NCBITaxa
from tqdm import tqdm

from SeqIterator import SeqReader

ncbi = NCBITaxa()


def check_distribution(fasta_file, index_file, rank, g_list=False):
    """
    Get the distribution of the taxids in the fasta file.

    Parameters
    ----------
    fasta_file: str
        The locaton of a Radogest fasta file.

    Returns
    -------
    sk_dict: dict
        A dictionary giving the counds of taxids.

    """
    index = pickle.load(open(index_file, "rb"))
    if g_list:
        seq_reader = open(fasta_file, "r")
    else:
        seq_reader = SeqReader(fasta_file, file_type="fasta")
    rank_dict = defaultdict(int)
    total = 0
    for record in tqdm(seq_reader):
        if g_list:
            genome_acc = record.strip().split()[0]
        else:
            genome_acc = record[0].split(":")[0]
        total += 1
        try:
            taxid = int(index["genomes"][genome_acc]['species_taxid'])
            ranks = ncbi.get_rank(ncbi.get_lineage(taxid))
        except (ValueError, KeyError) as e:
            rank_dict[-1] += 1
            continue
        for tid in ranks:
            if ranks[tid] == rank:
                rank_dict[tid] += 1
                break
        else:
            rank_dict[0] += 1
    return rank_dict, total


def main():
    """Parse the arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__)
    parser.add_argument(
        "file",
        type=str,
        help=("The Radogest fasta file or a list of genomes to check."))
    parser.add_argument("--index",
                        "-i",
                        type=str,
                        help=("The location of the index."),
                        default="./index_initial.pck")
    parser.add_argument("--rank",
                        "-r",
                        type=str,
                        help=("The rank to calculate a distribution for."),
                        default="superkingdom")
    parser.add_argument(
        "--list",
        "-l",
        action="store_true",
        help="Check the distribution from a list of genome accessoins.",
        default=False)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    rank_dict, total = check_distribution(args.file, args.index, args.rank,
                                          args.list)
    print(total, file=sys.stdout)
    print(json.dumps(rank_dict), file=sys.stdout)
    toc = datetime.datetime.now()
    print("The process took time: {}.".format(toc - tic), file=sys.stderr)


if __name__ == "__main__":
    main()
