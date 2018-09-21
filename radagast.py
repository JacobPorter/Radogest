#!/usr/bin/env python
"""
Random genome sampler for trees.

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
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_download = subparsers.add_parser("download",
                                       help=("Download genomes from NCBI."),
                                       formatter_class=argparse.
                                       ArgumentDefaultsHelpFormatter)
    p_parallel_sample = subparsers.add_parser("parallel_sample",
                                              help=(""),
                                              formatter_class=argparse.
                                              ArgumentDefaultsHelpFormatter)
    p_parallel_sample.add_argument(*index.args, **index.kwargs)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stdout.flush()
    mode = args.mode
    if mode == "download":
        pass
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
