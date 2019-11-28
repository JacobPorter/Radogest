#!/usr/bin/env python
"""
Print out the number of each taxid in a
Radogest generated fasta file.

:Authors:
    Jacob Porter <jsporter@virginia.edu>

"""
import argparse
import datetime
import pickle
import sys
from collections import defaultdict

from ete3 import NCBITaxa

from SeqIterator import SeqIterator

ncbi = NCBITaxa()


def check_distribution(fasta_file, index_file, rank):
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
    seq_reader = SeqIterator(fasta_file, file_type="fasta")
    rank_dict = defaultdict(int)
    total = 0
    for record in seq_reader:
        genome_acc = record[0].split(":")[0]
        total += 1
        try:
            taxid = int(index["genomes"][genome_acc]['species_taxid'])
            ranks = ncbi.get_rank(ncbi.get_lineage(taxid))
        except (ValueError, KeyError) as e:
            continue
        for tid in ranks:
            if ranks[tid] == rank:
                rank_dict[tid] += 1
    return rank_dict


def main():
    """Parse the arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__)
    parser.add_argument("fasta_file",
                        type=str,
                        help=("The Radogest fasta file to check."))
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
    args = parser.parse_args()
    print(args, file=sys.stderr)
    counts = check_distribution(args.fasta_file, args.index, args.rank)
    print(counts, file=sys.stdout)
    toc = datetime.datetime.now()
    print("The process took time: {}.".format(toc - tic), file=sys.stderr)


if __name__ == "__main__":
    main()
