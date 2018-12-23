#!/usr/bin/env python
"""
Randomly perform the reverse complement on some DNA strings from a fasta file.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""
import argparse
import datetime
import sys
import random
from SeqIterator.SeqIterator import SeqReader, SeqWriter

# The reverse complement mapping.
reverse_mapping = {"A": "T",
                   "T": "A",
                   "C": "G",
                   "G": "C",
                   "N": "N",
                   "Y": "R",
                   "R": "Y",
                   "W": "W",
                   "S": "S",
                   "K": "M",
                   "M": "K",
                   "D": "H",
                   "V": "B",
                   "H": "D",
                   "B": "V",
                   "X": "X",
                   "-": "-",
                   }


def get_reverse_complement(seq):
    """
    Calculate the reverse complement of a DNA string.

    Parameters
    ----------
    seq: iterable
        An ordered list of DNA bases.

    Returns
    -------
    A string of DNA bases that represent the reverse complement of the input.

    """
    return "".join([reverse_mapping[base.upper()] for base in reversed(seq)])


def get_rc_fasta(filename_input, filename_output, prob=0.5):
    """
    Randomly do the reverse comlement for DNA sequences in a fasta file.

    Parameters
    ----------
    filename_input: str
        The location of the filename.
    filename_output: str
        The location of the output.  If None or False, sys.stdout will be used.
    prob: float
        A number between 0 and 1.0 that represents the probability that the
        reverse complement will be done.

    Returns
    -------
    counter: int
        A count of the records processed.

    """
    input_iterator = SeqReader(filename_input, file_type="fasta")
    output_fd = open(filename_output, 'w') if filename_output else sys.stdout
    counter = 0
    output_writer = SeqWriter(output_fd, file_type="fasta")
    for seq_id, seq_seq in input_iterator:
        add_id = ""
        if random.random() < prob:
            try:
                seq_seq = get_reverse_complement(seq_seq)
            except KeyError:
                print(counter, seq_id, seq_seq)
                raise KeyError
            add_id = "_RC"
        output_writer.write((seq_id + add_id, seq_seq))
        counter += 1
    return counter


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Randomly converts '
                                                  'sequences to the reverse '
                                                  'complement from an input '
                                                  'fasta file.  Output is '
                                                  'a fasta file.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta_filename", type=str,
                        help=("The name of the input fasta filename."))
    parser.add_argument("-p", "--prob", type=float,
                        help=("The probability that a sequence will be "
                              "converted to the reverse complement."),
                        default=0.5)
    parser.add_argument("-o", "--output", type=str,
                        help=("The output file location to write the "
                              "fasta file to. "
                              "Defaults to stdout."),
                        default=None)
    args = parser.parse_args()
    number_written = get_rc_fasta(args.fasta_filename, args.output,
                                  prob=args.prob)
    tock = datetime.datetime.now()
    print("Records written: {}".format(number_written), file=sys.stderr)
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == '__main__':
    main()
