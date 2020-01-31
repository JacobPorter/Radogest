"""
Chop genomes into smaller sequences with a possibly overlapping window.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import operator
import os
import random
import sys

from SeqIterator.SeqIterator import SeqReader, SeqWriter

# The default seed for the random number generator for subsampling.
_SEED = 42


def chop_a_genome(location,
                  accession,
                  length,
                  taxid,
                  species_taxid,
                  writer,
                  queue=None,
                  include_wild=False,
                  window_length=50,
                  subsample=0,
                  sub_cutoff=0.2,
                  seed=_SEED):
    """
    Cut up a single genome into kmers and write to a file.

    Parameters
    ----------
    location: str
        The location of a fasta file.
    accession: str
        The accession id of a genome in the genomes directory.
    length: int
        The length of the kmer to get.
    taxid: int
        The taxid to use for the fasta record.
    species_taxid: int
        The species taxid associated with the genome accession.
    writer: SeqWriter
        The SeqWriter object to write fasta records to.
    queue: Queue
        A shared queue to write output to.
        The queue handles writing to files for multiprocessing.
    include_wild: boolean
        Determines whether to include kmers with wildcard characters.
    window_length: int
        The amount to slide the window when taking kmers.
    subsample: int
        Determine whether random subsampling should occur.
        If 0, then do not subsample.
        If 1, then subsample where samples are taken
        if the random number is >= sub_cutoff
        If -1, then subsample where samples are taken
        if a random number is < sub_cutoff
    sub_cutoff: float
        A number between [0, 1.0).
        This determines the percentage of kmers to sample.
    seed: str or float or int
        The object to seed the random number generator.
        This allows for a deterministic sampling procedure.


    Returns
    -------
    int
        The number of records written.

    """
    if subsample:
        random.seed(seed)
    if subsample > 0:
        op = operator.ge
    elif subsample < 0:
        op = operator.lt
    number_written = 0
    only_files = [
        f for f in os.listdir(location)
        if os.path.isfile(os.path.join(location, f)) and f.endswith('fna')
    ]
    if len(only_files) == 0:
        print("Skipping {}.  Fasta file not found.".format(accession),
              file=sys.stderr)
        return 0
    location = os.path.join(location, only_files[0])
    reader = SeqReader(location, file_type='fasta')
    for header, sequence in reader:
        for i in range(0, len(sequence), window_length):
            substring = sequence[i:i + length].upper()
            if (len(substring) == length
                    and (include_wild or 'N' not in substring)):
                seq_id = "{}:{}:{}:{}:{}-{}".format(accession, taxid,
                                                    species_taxid,
                                                    header.split()[0], i,
                                                    i + length)
                keep = True
                if subsample and not op(random.random(), sub_cutoff):
                    keep = False
                if queue and keep:
                    queue.put((seq_id, substring, str(taxid)))
                if writer and keep:
                    writer.write((seq_id, substring))
                if keep:
                    number_written += 1
    return number_written


def get_output_writer(output):
    """
    Get a SeqWriter object.

    Parameters
    ----------
    output: str, sys.stdout
        The type of writable object to use.
        Either a file location or sys.stdout.

    Returns
    -------
    A SeqWriter object.

    """
    if isinstance(output, str):
        if not output:
            output = sys.stdout
        else:
            output = open(output, 'w')
    if isinstance(output, SeqWriter):
        writer = output
    else:
        writer = SeqWriter(output, file_type='fasta')
    return writer


def chop_genomes(accessions_list,
                 length,
                 locations,
                 taxid,
                 species_taxid,
                 output,
                 queue=None,
                 include_wild=False,
                 window_length=50,
                 subsample=0,
                 sub_cutoff=0.2,
                 seed=_SEED):
    """
    Chop the DNA strings from fasta files from a list of genomes and
    write them to a fasta file.  This is a higher level function.

    Parameters
    ----------
    accessions_list: iterable
        A list of genome accessions.
    length: int
        The length of the partition to take.
    locations: list<str>
        A list of locations of the accession information
        in the same order as the accessions_list
    taxid: int
        The taxid to use for the fasta record.
    species_taxid: int
        The species taxid associated with the genome accession.
    output: str or writable
        The location of the output file as a string path or a writable object.
    queue: Queue
        A shared queue to write output to.
        The queue handles writing to files for multiprocessing.
    include_wild: boolean
        Determines whether to include kmers with wildcard characters.
    window_length: int
        The amount to slide the window when taking kmers.
    subsample: int
        Determine whether random subsampling should occur.
        If 0, then do not subsample.
        If 1, then subsample where samples are taken
        if the random number is >= sub_cutoff
        If -1, then subsample where samples are taken
        if a random number is < sub_cutoff
    sub_cutoff: float
        A number between [0, 1.0).
        This determines the percentage of kmers to sample.
    seed: str or float or int
        The object to seed the random number generator.
        This allows for a deterministic sampling procedure.

    Returns
    -------
    int
        The number of records written.

    """
    if window_length < 0:
        window_length = length
    if queue:
        writer = None
    else:
        writer = get_output_writer(output)
    number_written = 0
    for location, accession in zip(locations, accessions_list):
        number_written += chop_a_genome(location,
                                        accession,
                                        length,
                                        taxid,
                                        species_taxid,
                                        writer,
                                        queue=queue,
                                        include_wild=include_wild,
                                        window_length=window_length,
                                        subsample=subsample,
                                        sub_cutoff=sub_cutoff,
                                        seed=seed)
    return number_written
