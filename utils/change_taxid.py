#!/usr/bin/env python
"""Change a taxid into another taxid in a taxid file."""
import argparse
import sys
import datetime


def change(file_location, taxids, change):
    """Change the taxids from a taxid file."""
    fd = open(file_location, 'r')
    count = 0
    for line in fd:
        count += 1
        input_taxid = int(line)
        if input_taxid in taxids:
            print(change, file=sys.stdout)
        else:
            print(input_taxid, file=sys.stdout)
    return count


def main():
    """Parse the arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Change one taxid into '
                                     'another in a taxid file.'))
    parser.add_argument("taxids", type=int, nargs='+',
                        help=("The taxids to change."))
    parser.add_argument("--change", "-c", type=int, default=-1,
                        help=("Taxid to change into."))
    parser.add_argument("--input", "-i", type=str, default="1.taxid",
                        help=("The location of the taxid to change."))
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    count = change(args.input, args.taxids, args.change)
    print("There were {} records processed.".format(count), file=sys.stderr)
    toc = datetime.datetime.now()
    print("The process took time: {}".format(toc - tic), file=sys.stderr)


if __name__ == '__main__':
    main()
