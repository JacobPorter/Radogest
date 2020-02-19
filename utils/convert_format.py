#!/usr/bin/env python
"""
Converts Radogest/Plinko fasta and taxid files to a different format.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>
"""
import argparse
import csv
import datetime
import sys
from tqdm import tqdm

from SeqIterator import SeqReader, SeqWriter


def get_files_vm(fasta_file, taxid_file, output_file, taxid=10239):
    """
    Get sequences, etc. from a Radogest/Pinko fasta file in a ViraMiner file.
    """
    csvfile = open(output_file, 'w') if isinstance(output_file,
                                                   str) else output_file
    vm_writer = csv.writer(csvfile, delimiter=',')
    fasta_reader = SeqReader(fasta_file)
    fasta_records = [record for record in fasta_reader]
    fasta_reader.close()
    with open(taxid_file, "r") as taxid_fd:
        taxids = [int(taxid.strip()) for taxid in taxid_fd]
    for i, record in tqdm(enumerate(fasta_records)):
        if taxids[i] == taxid:
            vm_writer.writerow([record[0], record[1], 1])
        else:
            vm_writer.writerow([record[0], record[1], 0])


def get_files_dvf(fasta_file,
                  taxid_file,
                  extract_file,
                  remainder_file,
                  extract_taxid=10239):
    """
    Get the files in DeepVirFinder and VirFinder format.
    """
    extract_w = SeqWriter(extract_file)
    remainder_w = SeqWriter(open(remainder_file, "w"))
    with open(fasta_file, "r") as fasta_fd:
        fasta_reader = SeqReader(fasta_fd)
        fasta_records = [record for record in fasta_reader]
    with open(taxid_file, "r") as taxid_fd:
        taxids = [int(taxid.strip()) for taxid in taxid_fd]
    for i, record in enumerate(fasta_records):
        if taxids[i] == extract_taxid:
            extract_w.write(record)
        else:
            remainder_w.write(record)


def main():
    """Parse arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Copy data.'))
    parser.add_argument("fasta_file",
                        type=str,
                        help=("The location of a Radogest fasta file."))
    parser.add_argument("taxid_file",
                        type=str,
                        help=("The location of a taxid file "
                              "corresponding to the fasta file."))
    parser.add_argument("--output",
                        "-o",
                        type=str,
                        help=("Where to print the output file."),
                        default=sys.stdout)
    parser.add_argument("--remainder",
                        "-r",
                        type=str,
                        help=("The remainder file location if any.  "
                              "This is for DeepVirFinder"),
                        default=None)
    parser.add_argument("--type",
                        type=str,
                        choices=['viraminer', 'deepvirfinder'],
                        help=("Choose the output format."),
                        default='deepvirfinder')
    parser.add_argument("--taxid",
                        "-t",
                        type=int,
                        help=("The taxid to extract."),
                        default=10239)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    if args.type == 'viraminer':
        get_files_vm(args.fasta_file,
                     args.taxid_file,
                     args.output,
                     taxid=args.taxid)
    else:
        get_files_dvf(args.fasta_file,
                      args.taxid_file,
                      args.output,
                      args.remainder,
                      extract_taxid=args.taxid)
    toc = datetime.datetime.now()
    print("The process took time: {}".format(toc - tic), file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == "__main__":
    main()
