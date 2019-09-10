#!/bin/env python
"""
Calculate the number of overlapping sequences in Radogest generated data.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>

"""
import argparse
import datetime
import sys
from SeqIterator import SeqReader
from collections import defaultdict


def compute_overlap(fasta_file, fasta_compare=None, exact=False):
    """
    Calculate the number of overlapping sequences.

    Parameters
    ----------
    fasta_file: str
        The location of a fasta file.

    Returns
    -------
    float, int, int
        The percentage of sequences that overlap.
        The count of sequences in the fasta_file
        The count of sequences analyzed.

    """
    reader = SeqReader(fasta_file)
    seq_dict = defaultdict(list)
    count1 = 0
    for record in reader:
        fasta_id = record[0].split(":")
        accession = fasta_id[0]
        contig = fasta_id[2]
        location = tuple(map(int, fasta_id[3].split("-")))
        seq_dict[accession + ":" + contig].append(location)
        count1 += 1
    reader.close()
    if fasta_compare:
        reader_compare = SeqReader(fasta_compare)
        self_compare = False
    else:
        reader_compare = SeqReader(fasta_file)
        self_compare = True
    count2 = 0
    overlap = 0
    unique_overlap = 0
    unique_overlap_exact = 0
    overlap_exact = 0
    stop_print = 5
    for record in reader_compare:
        fasta_id = record[0].split(":")
        accession = fasta_id[0]
        contig = fasta_id[2]
        location = tuple(map(int, fasta_id[3].split("-")))
        loc_list = seq_dict[accession + ":" + contig]
        count2 += 1
        if self_compare:
            switch = 2
            switch_exact = 2
        else:
            switch = 1
            switch_exact = 1
        if loc_list:
            for loc in loc_list:
                if ((location[0] >= loc[0] and location[0] <= loc[1]) or
                        (location[1] >= loc[0] and location[1] <= loc[1])):
                    overlap += 1
                    if stop_print:
                        print(accession, contig, location, loc,
                              file=sys.stderr)
                        stop_print -= 1
                    if switch:
                        unique_overlap += 1
                        switch -= 1
                    if location[0] == loc[0]:
                        overlap_exact += 1
                        if switch_exact:
                            unique_overlap_exact += 1
                            switch_exact -= 1
    return (overlap, unique_overlap, overlap_exact, unique_overlap_exact,
            count1, count2)


def main():
    """Parse arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=(''))
    parser.add_argument("fasta_file", type=str,
                        help=("A fasta file from Radogest."))
    parser.add_argument("--compare_fasta", "-c", type=str,
                        help=("Another fasta file from Radogest to "
                              "compare to. "
                              "If comparing the fasta_file with itself, "
                              "do not use this option."),
                        default=None)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    stats = compute_overlap(args.fasta_file, args.compare_fasta)
    overlap, unique_overlap, overlap_ex, unique_overlap_ex, cnt1, cnt2 = stats
    if not args.compare_fasta:
        overlap = overlap - cnt1
        unique_overlap = unique_overlap - cnt1
        overlap_ex = overlap_ex - cnt1
        unique_overlap_ex = unique_overlap_ex - cnt1
    print("There were {},{},{} instances of overlap of non-unique "
          "sequences.".format(overlap, overlap/cnt1, overlap/cnt2),
          file=sys.stdout)
    print("There were {},{},{} unique sequences "
          "that overlap.".format(unique_overlap,
                                 unique_overlap/cnt1,
                                 unique_overlap/cnt2),
          file=sys.stdout)
    print("There were {},{},{} instances of exact overlapping non-unique "
          "sequences".format(overlap_ex, overlap_ex/cnt1, overlap_ex/cnt2),
          file=sys.stdout)
    print("There were {},{},{} unique sequences "
          "that exactly overlap.".format(unique_overlap_ex,
                                         unique_overlap_ex/cnt1,
                                         unique_overlap_ex/cnt2),
          file=sys.stdout)
    print("There were {} sequences in the fasta_file and {} sequences "
          "in the comparison file.".format(cnt1, cnt2),
          file=sys.stdout)
    toc = datetime.datetime.now()
    print("The process took time: {}".format(toc - tic), file=sys.stderr)


if __name__ == '__main__':
    main()
