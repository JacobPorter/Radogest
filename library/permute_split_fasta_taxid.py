#!/usr/bin/env python
"""
Randomly permutes a fasta file and a taxid file.  Can split the data into
training, validation, and testing partitions.

:Authors:
        Jacob Porter <jsporter@vt.edu>
"""

import argparse
import datetime
import sys
from random import shuffle
import SeqIterator


def write_fasta_taxid(both_records, fasta_writable, taxid_writable):
    """
    Write a fasta file and a taxid file from a list of both taxid and fasta
    records.

    Parameters
    ----------
    both_records: list
        A list of taxid records and fasta records where an entry in the list
        has the form (fasta_record, taxid)
    fasta_writable: writable
        A writable object where a fasta file can be written.
    taxid_writable: writable
        A writable object where a list of taxids can be written.

    Returns
    -------
    int
        A count of records written.

    """
    fasta_writer = SeqIterator.SeqWriter(fasta_writable)
    taxid_writer = taxid_writable
    count = 0
    for record in both_records:
        fasta_writer.write(record[0])
        taxid_writer.write(record[1] + "\n")
        count += 1
    fasta_writer.flush()
    taxid_writer.flush()
    fasta_writer.close()
    taxid_writer.close()
    return count


def randomly_permute_fasta_taxid(fasta_file, taxid_file, fasta_out,
                                 taxid_out, split=True, split_amount=0.8):
    """
    Randomly permute the fasta and taxid files.  May split the input into
    training, validation, and testing data.

    Parameters
    ----------
    fasta_file: str
        The location of a fasta file.
    taxid_file: str
        The location of a taxid file.
    fasta_out: str
        The location to store the output fasta file.
    taxid_out: str
        The location to store the output taxid file.
    split: bool
        If True, split the data.
    split_amount: tuple of floats or float
        Percentage(s) to split into.

    Returns
    -------
    int
        The number of records written.

    """
    both_records = []
    for fasta, taxid in zip(SeqIterator.SeqIterator(fasta_file,
                                                    file_type='fasta'),
                            open(taxid_file, 'r')):
        both_records.append((fasta, taxid.strip()))
    shuffle(both_records)
    print("The total number of records to write: {}".format(len(both_records)),
          file=sys.stderr)
    if split:
        print("Splitting the output.", file=sys.stderr)
        amounts = []
        if isinstance(split_amount, float):
            amounts.append(int(split_amount * len(both_records)))
        else:
            # print(list(map(float, split_amount.split(","))))
            # print(len(both_records))
            if sum(list(map(float, split_amount.split(",")))) > 1.0:
                print("The sum of the split amounts is larger than 1.0",
                      file=sys.stderr)
            for amount in list(map(float, split_amount.split(","))):
                amounts.append(round(amount * len(both_records)))
        begin = 0
        count = 0
        extensions = [".train", ".validate", ".test"]
        for i, amount in enumerate(amounts):
            ext = extensions[i] if i < len(extensions) else str(i)
            end = (begin+amount if
                   begin+amount <= len(both_records) else
                   len(both_records))
            count += write_fasta_taxid(both_records[begin:end],
                                       open(fasta_out + ext, 'w'),
                                       open(taxid_out + ext, 'w'))
            begin = end
        return count
    return write_fasta_taxid(both_records,
                             open(fasta_out, 'w'),
                             open(taxid_out, 'w'))


def main():
    """Parse arguments for randomly permuting the fasta and taxid files."""
    parser = argparse.ArgumentParser(description=('Randomly permute fasta and '
                                                  'taxid records.  Partitions '
                                                  'the data into training, '
                                                  'validation, and testing '
                                                  'files.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta_file", type=str, help=('The location of the '
                                                      'fasta file.'))
    parser.add_argument("taxid_file", type=str, help=('The location of the '
                                                      'taxid file.'))
    parser.add_argument("--split", "-s", action='store_true',
                        help=("Split the data into training, validation, and "
                              "testing sets."),
                        default=False)
    parser.add_argument("--split_amount", "-a", type=str,
                        help=("The percentage of the input to use for "
                              "training, validation, and testing "
                              "seperated by commas."),
                        default='0.8,0.1,0.1')
    args = parser.parse_args()
    now = datetime.datetime.now()
    sys.stderr.write("Permuting the fasta file {} "
                     "and the taxid file {}.\n".format(args.fasta_file,
                                                       args.taxid_file))
    if args.split:
        sys.stderr.write("Splitting:{}\nAmount:{}\n".format(args.split,
                                                            args.split_amount))
    count = randomly_permute_fasta_taxid(args.fasta_file,
                                         args.taxid_file,
                                         args.fasta_file ,
                                         args.taxid_file,
                                         args.split,
                                         args.split_amount)
    sys.stderr.write("{} records written.\n".format(count))
    later = datetime.datetime.now()
    sys.stderr.write("The process took time:\t{}.\n".format(later - now))


if __name__ == "__main__":
    main()
