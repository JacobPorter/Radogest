#!/usr/bin/env python
"""
Split a list of genomes into smaller sequences.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""
import argparse
import datetime
import pickle
import sys
import os
from SeqIterator import SeqIterator, SeqWriter
from get_paths import get_paths

paths = get_paths()
GENOMES = paths['GENOMES']


def split_genomes(accessions_list, length, index, output,
                  includeNs=False, window_length=50):
    """
    Partition the DNA strings in a list of genomes and write to a fasta file.

    Parameters
    ----------
    accessions_list: iterable
        A list of genome accessions.
    length: int
        The length of the partition to take.
    index: dict
        The genomes index.
    output: str or writable
        The location of the output file as a string path or a writable object.

    Returns
    -------
    int
        The number of records written.

    """
    if window_length < 0:
        window_length = length
    if isinstance(output, str):
        if not output:
            output = sys.stdout
        else:
            output = open(output, 'w')
    if isinstance(output, SeqWriter):
        writer = output
    else:
        writer = SeqWriter(output, file_type='fasta')
    number_written = 0
    for accession in accessions_list:
        location = GENOMES + index['genomes'][accession]['location']
        # print(location, file=sys.stderr)
        only_files = [f for f in os.listdir(location) if
                      os.path.isfile(os.path.join(location, f))
                      and f.endswith('fna')]
        if len(only_files) == 0:
            print("Skipping {}.  Fasta file not found.".format(accession),
                  file=sys.stderr)
            continue
        location = os.path.join(location, only_files[0])
        reader = SeqIterator(location, file_type='fasta')
        for header, sequence in reader:
            for i in range(0, len(sequence), window_length):
                substring = sequence[i:i+length].upper()
                if (len(substring) == length and (
                        includeNs or 'N' not in substring)):
                    writer.write(("{}_{}_[{}:{}]".format(
                        accession, header, i, i+length), substring))
                    number_written += 1
    return number_written


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Produces a fasta file of '
                                                  'all of the dna substrings '
                                                  'from a list of genomes.  '
                                                  'The fasta file is given '
                                                  'to stdout.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("accessions_list", type=str,
                        help=(),
                        default='./accessions_list.txt')
    parser.add_argument("--index", "-i", type=str,
                        help=("The location of the index that "
                              "maps taxonomic ids to genomes and "
                              "their locations, etc."),
                        default='./index.pck')
    parser.add_argument("--include_Ns", "-n", action='store_true',
                        help=("Include sequences with Ns in them."),
                        default=False)
    parser.add_argument("--length", "-l", type=int,
                        help=("The length of the substring to take."),
                        default=100)
    parser.add_argument("--window_length", "-w", type=int,
                        help=("The offset to use for sliding sampling over.  "
                              "If this equals the length, then samples will "
                              "not overlap.  If this equals 1, then all "
                              "kmers will be produced."),
                        default=50)
    parser.add_argument("--output", "-o",
                        help="The location of the output file.  "
                        "Defaults to stdout.",
                        default="")
    args = parser.parse_args()
    accessions_list = [line.strip() for line in open(args.accessions_list)]
    number_written = split_genomes(accessions_list, args.length,
                                   pickle.load(open(args.index, 'rb')),
                                   args.output, args.include_Ns,
                                   args.window_length)
    tock = datetime.datetime.now()
    print("Records written: {}".format(number_written), file=sys.stderr)
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == '__main__':
    main()
