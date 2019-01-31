import random
import sys
from SeqIterator.SeqIterator import SeqReader, SeqWriter
from collections import defaultdict

_RC_PROB = 0.5


reverse_mapping_init = {"A": "T",
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

reverse_mapping = defaultdict(lambda: "N")

for key in reverse_mapping_init:
    reverse_mapping[key] = reverse_mapping_init[key]


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


def random_reverse(seq_id, seq, prob=_RC_PROB):
    """
    Randomly compute the reverse complement of a sequence and modify the id.
    
    Parameters
    ----------
    seq_id: str
        The id of the sequence for a fasta record.
    seq: str
        The DNA sequence.
    prob: float 0 
        The probability of taking the reverse complement.
        0.0 <= prob <= 1.0
        
    Returns
    -------
    (str, str)
        A tuple representing a sequence id followed by a DNA sequence.
    
    """
    if random.random() >= prob:
        return (seq_id + ":+", seq)
    else:
        try:
            rc_seq = get_reverse_complement(seq)
        except KeyError:
            print((seq_id, seq), file=sys.stderr)
            raise KeyError
        return (seq_id + ":-", rc_seq)


def get_rc_fasta(filename_input, 
                 filename_output, 
                 prob=_RC_PROB, 
                 remove=False, 
                 verbose=0):
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
    remove: boolean
        If True, remove sequences with N's in them.

    Returns
    -------
    counter: int
        A count of the records processed.

    """
    input_iterator = SeqReader(filename_input, file_type="fasta")
    output_fd = open(filename_output, 'w') if filename_output else sys.stdout
    read_counter = 0
    write_counter = 0
    output_writer = SeqWriter(output_fd, file_type="fasta")
    for seq_id, seq_seq in input_iterator:
        read_counter += 1
        if remove and ("n" in seq_seq or "N" in seq_seq):
            continue
        output_writer.write(random_reverse(seq_id, seq_seq, prob=prob))
        write_counter += 1
        if verbose > 1 and write_counter % 1000 == 0:
            print("Reverse complement: {} records written so far.".format(
                write_counter),
                 file=sys.stderr)
    return read_counter, write_counter