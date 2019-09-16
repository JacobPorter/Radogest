#!/usr/bin/env python
"""
Print out the number of each taxid in a
Radogest generated fasta file.

:Authors:
    Jacob Porter <jsporter@virginia.edu>

"""
import argparse
import datetime
import sys
from SeqIterator import SeqIterator
from collections import defaultdict


def check_distribution(fasta_file):
    """
    Get the distribution of the taxids in teh fasta file.

    Parameters
    ----------
    fasta_file: str
        The locaton of a Radogest fasta file.

    Returns
    -------
    sk_dict: dict
        A dictionary giving the counds of taxids.

    """
    seq_reader = SeqIterator(fasta_file, file_type="fasta")
    sk_dict = defaultdict(int)
    for record in seq_reader:
        sk = int(record[0].split(":")[1])
        sk_dict[sk] += 1
    return sk_dict


def main():
    """Parse the arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__)
    parser.add_argument("fasta_file", type=str,
                        help=("The Radogest fasta file to check."))
    args = parser.parse_args()
    print(args, file=sys.stderr)
    counts = check_distribution(args.fasta_file)
    print(counts, file=sys.stdout)
    toc = datetime.datetime.now()
    print("The process took time: {}.".format(toc - tic), file=sys.stderr)


if __name__ == "__main__":
    main()
