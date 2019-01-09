"""
Split a list of genomes into smaller sequences.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import sys
import os
from SeqIterator.SeqIterator import SeqReader, SeqWriter

# from config import GENOMES_NT, GENOMES_AA, GENOMES_CD


def split_genomes(accessions_list, 
                  length, 
                  index, 
                  genomes_dir, 
                  output,
                  include_wild=False, 
                  window_length=50):
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
        location = genomes_dir + index['genomes'][accession]['location']
        # print(location, file=sys.stderr)
        only_files = [f for f in os.listdir(location) if
                      os.path.isfile(os.path.join(location, f))
                      and f.endswith('fna')]
        if len(only_files) == 0:
            print("Skipping {}.  Fasta file not found.".format(accession),
                  file=sys.stderr)
            continue
        location = os.path.join(location, only_files[0])
        reader = SeqReader(location, file_type='fasta')
        for header, sequence in reader:
            for i in range(0, len(sequence), window_length):
                substring = sequence[i:i+length].upper()
                if (len(substring) == length and (
                        include_wild or 'N' not in substring)):
                    writer.write(("{}_{}_[{}:{}]".format(
                        accession, header, i, i+length), substring))
                    number_written += 1
    return number_written