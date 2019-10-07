#!/usr/bin/env python
"""
Get a list of taxids that have data sizes above some threshold.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""
import argparse
import datetime
import sys
from ast import literal_eval


def process_file(file_location, threshold=10000000):
    """
    Process a balance_data.py file looking for taxids that have too much data.

    Parameters
    ----------
    file_location: str
        The location of a balance_data.py file.
    threshold: int
        The maximum data size to look for.

    Returns
    -------
    too_large: list
        A list of taxonomic ids

    """
    too_large = []
    with open(file_location) as fdr:
        for line in fdr:
            if line.startswith("("):
                line = line.strip()
                taxid_tuple = literal_eval(line)
                if taxid_tuple[1] == 'train' and taxid_tuple[2] > threshold:
                    too_large.append((taxid_tuple[0], taxid_tuple[2]))
    return too_large


def main():
    """Parse the arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        description=('Examine a balance_data.py output file and '
                     'look for taxids with data sizes that are too large.'))
    parser.add_argument("file",
                        type=str,
                        help=("The output of balance_data.py"))
    parser.add_argument("--threshold",
                        "-t",
                        type=int,
                        help=("The maximum size of a data set."),
                        default=10000000)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    output = process_file(args.file, args.threshold)
    for taxid in output:
        print(taxid[0], file=sys.stdout)
        print(taxid, file=sys.stderr)
    toc = datetime.datetime.now()
    print("The process took time {}.".format(toc - tic), file=sys.stderr)


if __name__ == "__main__":
    main()
