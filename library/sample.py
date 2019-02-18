"""
Produces a fasta file of randomly chosen kmer sequences corresponding to the
given taxonomic id(s).

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import os
import errno
import pickle
import shutil
import sys
import random
import string
import subprocess
import sqlite3
import tempfile
import operator
from multiprocessing import Pool
from collections import defaultdict

from ete3 import NCBITaxa

from SeqIterator.SeqIterator import SeqReader, SeqWriter
from library.permute import randomly_permute_fasta_taxid
from library.chop import chop_genomes

from config import BEDTOOLS

ncbi = NCBITaxa()

# If the sample falls below NUMBER_CUTOFF, get NUMBER_SAMPLE samples instead.
NUMBER_CUTOFF = 150
NUMBER_SAMPLE = 300
# The actual number of random samples taken is multipled by this multiplier.
# This is done so that samples with N's in them can be excluded.
SAMPLE_MULTIPLIER = 1.2

# Standard deviations to take for kmer size checking
_STD_DEV = 1.5

# KMERS above this size will be checked for their wildcard percentage.
_WILDCARD_KMER_T = 1000
# The number of kmers to sample when checking the wildcard percentage.
_WILDCARD_SAMPLE_NUM = 500
# The percentage of kmers that have wildcards before the genome is discarded.
_WILDCARD_PERCENT_T = 0.70

# The default probability of taking the reverse complement of a DNA sequence.
_RC_PROB = 0.5

# For the genome holdout selection strategy, these constants indicate
# the data set that a genome will be in.
_TRAIN = 1
_TEST = 2

# The base for the thresholds amount.  1000 means the value is in kilobases.
THRESHOLD_BASE = 1000

# The reverse complement mapping.
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
    verbose: int
        Controls the verbosity.

    Returns
    -------
    counter: int, int
        A count of the records read and written.

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


def binary_search(target, array, begin, end):
    """
    Perform a generic binary search for an array.

    Parameters
    ----------
    target: comparable
        An item to compare to
    array: iterable, indexable
        An array to search for
    begin: int
        Beginning index for the array
    end: int
        Ending index for the array

    Returns
    -------
    An item that satisfies the binary search.

    """
    if array == []:
        return None
    if begin == end:
        return begin
    if begin == end-1:
        return begin if target < array[begin] else end
    mid = begin + int(round((end - begin) / 2.0))
    if target < array[mid]:
        return binary_search(target, array, begin, mid)
    return binary_search(target, array, mid, end)


def uniform_samples_at_rank(index, sublevels, genomes_dir,
                            number, kmer_length,
                            include_wild, amino_acid,
                            temp_dir, include_list, threshold=None):
    """
    Get a count of samples from each genome under the taxids given by ranks.

    Parameters
    ----------
    index: dict
        A dictionary representing genome locations and genome taxid
        associations
    sublevels: iterable
        An iterable that gives ranks underneath taxid to sample from.
    genomes_dir: str
        The location of the directory where the genomes are stored.
    number: int
        The number of random sequences to sample.
    kmer_length: int
        The length of the kmer to sample.
    include_wild: boolean
        True if wildcard DNA characters are desired.  False otherwise.
    amino_acid: boolean
        Set to True if sampling amino acid data.
    temp_dir: str
        The location of a directory to write temporary files.
    include_list: list
        Determines the value from the genomes index that determines
        if the genome will be included for consideration.  Ordinarily,
        this will be set to [True], but for the genome holdout strategy,
        it could be [1] or [2].
    threshold: int
        A value that controls how much genomic content to include.

    Returns
    -------
    dict
        A dictionary keyed on taxids where the value is the number of sequences
        to sample for each taxid.

    """
    taxid_accessions = {}
    for taxid in sublevels:
        my_accessions = [accession
                         for accession in index['taxids'][taxid]
                         if include_accession(accession,
                                              taxid,
                                              index,
                                              genomes_dir,
                                              kmer_length,
                                              include_wild,
                                              amino_acid,
                                              temp_dir,
                                              include_list)]
        my_sums = [index['genomes'][accession]['contig_sum']
                   for accession in my_accessions]
        if threshold and my_accessions:
            threshold_actual = threshold * THRESHOLD_BASE
            genomic_content = 0
            for i in range(len(my_accessions)):
                if genomic_content > threshold_actual:
                    break
                genomic_content += my_sums[i]
            my_accessions = my_accessions[0:i]
            my_sums = my_sums[0:i]
        if my_accessions:
            taxid_accessions[taxid] = [my_accessions,
                                       my_sums]
    if not taxid_accessions:
        return None
    try:
        uniform_number = round(number / len(taxid_accessions))
    except ZeroDivisionError:
        return None
    uniform_sample_counts = []
    for taxid in taxid_accessions:
        uniform_sample_counts.append((taxid,
                                      uniform_samples(taxid,
                                                      taxid_accessions[taxid],
                                                      uniform_number)))
    return uniform_sample_counts


def include_accession(accession, taxid, index, genomes_dir,
                      kmer_length, include_wild,
                      amino_acid, temp_dir="/localscratch",
                      include_list=[True]):
    """
    Determine whether to include a genome in the sampling

    Parameters
    ----------
    accession: str
        The genome accession id
    taxid: int
        The taxonomic id where sampling is desired.
        The taxid must be in the genome's species lineage.
    genomes_dir: str
        The location of the directory where the genomes are stored.
    index: dict
        The genomes index created
    kmer_length: int
        A positive integer representing the kmer length desired.
    include_wild: boolean
        True if wildcard DNA characters are desired.  False otherwise.
    amino_acid: boolean
        Set to True if sampling amino acid data.
    temp_dir: str
        The location of a directory to write temporary files.
    include_list: list
        Determines the value from the genomes index that determines
        if the genome will be included for consideration.  Ordinarily,
        this will be set to [True], but for the genome holdout strategy,
        it could be [1] or [2].


    Returns
    -------
    boolean
        Return True if the genome should be included.
        Return False if the genome should not be included.

    """
    if not index['taxids'][taxid][accession] in include_list:
        return False
    mean = index['genomes'][accession]['contig_mean']
    std = index['genomes'][accession]['contig_std']
    mx = index['genomes'][accession]['contig_max']
    cnt = index['genomes'][accession]['contig_count']
    if cnt > 5:
        inside_std = (kmer_length >= mean - _STD_DEV * std and
                      kmer_length <= mean + _STD_DEV * std and
                      kmer_length <= mx)
    else:
        inside_std = kmer_length <= mx
    if (kmer_length > _WILDCARD_KMER_T and not include_wild
            and not amino_acid and inside_std):
        file_locations_d = file_locations(accession, genomes_dir,
                                          index, temp_dir)
        fai_location = file_locations_d["fai"]
        fasta_location = file_locations_d["fasta_location"]
        taxid_file = None
        final_file = None
        (_,
         records_written,
         records_with_n) = get_random_bed_fast(_WILDCARD_SAMPLE_NUM,
                                               kmer_length,
                                               taxid,
                                               accession,
                                               fai_location,
                                               fasta_location,
                                               taxid_file,
                                               final_file,
                                               include_wild=include_wild,
                                               amino_acid=amino_acid,
                                               temp_dir=temp_dir)
        return (records_with_n / records_written < _WILDCARD_PERCENT_T)
    return inside_std


def uniform_samples(taxid, accession_sum, number):
    """
    Get a count of samples from each genome under taxid.  Does not do any
    sampling.

    Parameters
    ----------
    taxid: int
        A taxonomic id number to query from
    accession_sum: iterable
        This data structure contains two items: a list of accessions
        and a corresponding list of genome lengths.
    number: int
        The number of random sequences to sample.

    Returns
    -------
    dict
        A dictionary keyed on accession ids.  The value is the number of
        sequences to sample from the corresponding genome.

    """
    accessions = accession_sum[0]
    genome_lengths = accession_sum[1]
    try:
        name = ncbi.get_taxid_translator([taxid])[taxid]
    except KeyError:
        name = ""
    except sqlite3.DatabaseError:
        name = ""
    sys.stderr.write("For {}:{}, will sample from "
                     "{} fasta files.\n".format(taxid, name, len(accessions)))
    sys.stderr.flush()
    for i, _ in enumerate(genome_lengths):
        if i == 0:
            continue
        genome_lengths[i] += genome_lengths[i-1]
    try:
        total_length = genome_lengths[-1]
    except IndexError as ie:
        print("Genome lengths are empty for taxid: {}".format(taxid),
              file=sys.stderr)
        raise ie
    accession_counts = defaultdict(int)
    for _ in range(number):
        i = binary_search(random.randint(0, total_length - 1),
                          genome_lengths, 0, len(genome_lengths) - 1)
        accession_counts[accessions[i]] += 1
    return accession_counts


def file_locations(accession, genomes_dir, index, temp_dir):
    """
    Get file locations for sampling from the fasta file.

    Parameters
    ----------
    accession: str
        The accession id of the genome.
    genomes_dir: str
        The location where the genomes or fasta files are stored.
    index: dict
        The genomes dictionary index created by Radogest.
    temp_dir: str
        The temporary directory to store files

    Returns
    -------
    dict
        A dictionary of file locations.

    """
    rand_string = "".join(random.choices(
        string.ascii_letters + string.digits, k=12))
    my_fasta = os.path.join(temp_dir,
                            accession + "_" +
                            rand_string + "_random.fa")
    accession_location = os.path.join(genomes_dir +
                                      index['genomes']
                                      [accession]
                                      ['location'])
    onlyfiles = [f for f in os.listdir(accession_location) if
                 os.path.isfile(os.path.join(accession_location, f))]
    fasta_location = None
    twobit_location = None
    fai_location = None
    for f in onlyfiles:
        if (f.endswith(".fna") or
                f.endswith(".fasta") or
                f.endswith(".fa") or
                f.endswith(".faa")):
            fasta_location = os.path.join(accession_location, f)
        elif f.endswith(".2bit"):
            twobit_location = os.path.join(accession_location, f)
        elif f.endswith(".fai"):
            fai_location = os.path.join(accession_location, f)
    bedtools_file = os.path.join(temp_dir, accession + "_" +
                                 rand_string + "_random.bed")
    return {"my_fasta": my_fasta, "fasta_location": fasta_location,
            "twobit": twobit_location, "fai": fai_location,
            "bed": bedtools_file}


def get_fasta(accession_counts_list, length, index, genomes_dir,
              output, taxid_file, include_wild=False,
              window_length=50,
              thresholding=False, chop=False,
              amino_acid=False,
              temp_dir='/localscratch/', verbose=0):
    """
    Save randomly sampled sequences in a fasta file written to output.

    Parameters
    ----------
    accession_counts_list: list
        A list of tuples.  The first element is the taxid associated with the
        genomes in the second element.  The second element is a dictionary
        keyed on accession number, and the value is the number of samples to
        take from that genome.
    length: int
        The number of bases to sample.  (The length of the string.)
    index: dict
        A dictionary representing information about genomes and taxids.
    genomes_dir: str
        The location of the root of where the fasta files are stored.
    output: writable
        A writable object to store fasta files too.  This could be an open
        file.
    taxid_file: writable
        A file-like object to store taxonomic ids that correspond to each
        sampled sequence.  Each taxonomic id will be on its own line.
    include_wild: boolean
        Determines if sequences with wild card characters will be kept.
    window_length: int
        The length of the offset for the sliding window
        if thresholding or chopping is used.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.
    chop: bool
        If True, chop up the genome into kmers based on a sliding window.
    amino_acid: bool
        If True, the data is amino acid data.
    temp_dir: str
        A path to a directory to store temporary files.
    verbose: bool
        If True, print messages.

    Returns
    -------
    int
        A count of the number of fasta records written.

    """
    final_file = SeqWriter(output, file_type='fasta')
    fasta_record_count = 0
    for taxid, accession_counts in accession_counts_list:
        for accession in accession_counts.keys():
            if accession_counts[accession] <= 0:
                continue
            if verbose > 1:
                sys.stderr.write("Writing fasta records for "
                                 "taxid {} from genome accession {}:\t".
                                 format(taxid, accession))
            # A thresholding and chopping feature.
            # If thresholding and the genome is too small,
            # use the whole genome.
            # Or chop the genome up.
            if (thresholding and accession_counts[accession] > float(
                    index['genomes'][
                        accession]['contig_sum']) / length) or chop:
                records_written = chop_genomes([accession], length,
                                               index, genomes_dir,
                                               taxid,
                                               final_file,
                                               include_wild=include_wild,
                                               window_length=window_length)
                for _ in range(records_written):
                    taxid_file.write(str(taxid) + "\n")
                fasta_record_count += records_written
                continue
            file_locations_d = file_locations(accession, genomes_dir,
                                              index, temp_dir)
            fai_location = file_locations_d["fai"]
            fasta_location = file_locations_d["fasta_location"]
            get_random = True
            accession_cnt = 0
            while get_random:
                accession_number = accession_counts[accession] - accession_cnt
                bed_2bit_counts = get_random_bed_fast(accession_number,
                                                      length,
                                                      taxid,
                                                      accession,
                                                      fai_location,
                                                      fasta_location,
                                                      taxid_file,
                                                      final_file,
                                                      include_wild,
                                                      amino_acid,
                                                      temp_dir)
                get_random = not bed_2bit_counts[0]
                accession_cnt += bed_2bit_counts[1]
                if verbose > 1:
                    sys.stderr.write("amount: " + str(bed_2bit_counts[1]) +
                                     "  " +
                                     "N freq: " + str(bed_2bit_counts[2] /
                                                      (bed_2bit_counts[1] +
                                                       bed_2bit_counts[2])) +
                                     " ")
                    sys.stderr.flush()
            fasta_record_count += accession_cnt
            if verbose:
                sys.stderr.write("\n")
                sys.stderr.flush()
    return fasta_record_count


def get_random_bed_fast(number, length, taxid, accession, fai_location,
                        fasta_location, taxid_file, final_file,
                        include_wild=False, amino_acid=False,
                        temp_dir="/localscratch/"):
    """
    Get random nucleotide sequences from a bed file and a fasta file.  Exclude
    sequences with N's in them.

    Parameters
    ----------
    number: int
        The number of samples to take.
    length: int
        The length of the nucleotide sequence to get.
    taxid: int
        The taxonomic id to sample from
    accession: str
        A string indicating the accession id for a genome
    fai_location: str
        The location of the .fai faidx samtools index for the genome
    fasta_location: str
        The location of the fasta genome file
    taxid_file: writable
        A writable object where taxids are written for each random sample.
    final_file: writable
        A writable object that stores the fasta records.  This file represents
        the end product of all of the random sampling.
    include_wild: boolean
        Determines if sequences with wild card characters will be kept.
    amino_acid: bool
        If True, the data is amino acid data.
    temp_dir: str
        A path to a directory to store temporary files.

    Returns
    -------
    (boolean, int, int)
        A boolean to indicate if the records written is equal to number.
        An integer counting the records written,
        and an integer of the records that had 'N' characters.

    """
    taxid = str(taxid)
    prefix = taxid + "_" + accession + "_"
    if number <= NUMBER_CUTOFF:
        my_sample = NUMBER_SAMPLE
    else:
        my_sample = int(number * SAMPLE_MULTIPLIER)
    with tempfile.NamedTemporaryFile(
        mode="w+",
        suffix=".bed",
        prefix=prefix,
        dir=temp_dir) as bedtools_fd, tempfile.NamedTemporaryFile(
            mode="w+",
            suffix=".fasta",
            prefix=prefix,
            dir=temp_dir) as my_fasta_fd:
        bedtools_file = bedtools_fd.name
        subprocess.run([BEDTOOLS + "bedtools", "random", "-l",
                        str(length), "-n",
                        str(my_sample), "-g",
                        fai_location], stdout=bedtools_fd)
        subprocess.run([BEDTOOLS + "bedtools", "getfasta", "-fi",
                        fasta_location, "-bed", bedtools_file],
                       stdout=my_fasta_fd)
        intermediate_fasta_file = SeqReader(my_fasta_fd.name,
                                            file_type='fasta')
        records_with_n = 0
        records_written = 0
        for fasta_record in intermediate_fasta_file:
            if records_written >= number:
                break
            record_id, record_seq = fasta_record
            record_seq = record_seq.upper()
            if not amino_acid and "N" in record_seq:
                records_with_n += 1
                if not include_wild:
                    continue
            record_id = accession + ":" + taxid + ":" + record_id
            if final_file:
                final_file.write((record_id, record_seq))
            if taxid_file:
                taxid_file.write(taxid + "\n")
            records_written += 1
        intermediate_fasta_file.close()
    return (records_written >= number, records_written, records_with_n)


def get_sample(taxid, sublevels, index_dir, genomes_dir,
               number, length, data_dir,
               split=True, split_amount='0.8,0.1,0.1',
               include_wild=False,
               prob=_RC_PROB,
               thresholding=False,
               chop=False,
               window_length=50,
               amino_acid=False,
               thresholds=None,
               temp_dir="/localscratch/",
               verbose=0):
    """
    Get a random sample.  Create training, validation, and testing data sets
    and put them in the appropriate folders.
    Assesses the genome sampling strategy.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id
    sublevels: iterable
        A set of taxonomic ids below the taxid to sample from.
    index: dict
        The genomes index object.
    genomes_dir: str
        The location of the root of where the fasta files are stored.
    number: int
        The number of samples to take.
    length: int
        The number of bases for each sample
    data_dir: str
        The path to the data directory where fasta files will be written.
    split: bool
        Determine whether to split the data or not.
    split_amount: str
        A comma seperated list of floats representing the percentage of the
        data to be used for training, validation, and testing data
    include_wild: boolean
        Determines if sequences with wild card characters will be kept.
    prob: float
        A number between 0 and 1.0 that represents the probability that the
        reverse complement will be done.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.
    chop: bool
        If True, chop up the genome into kmers based on a sliding window.
    window_length: int
        The length of the offset for the sliding window
        if thresholding or chopping is used.
    amino_acid: bool
        If True, the data is amino acid data.
    thresholds: list<int>
        Values that control how much genomic content to include.  
        Multiple values are only valid for genome holdout strategies.
    temp_dir: str
        A path to write temporary files to.  This could be on the local hard
        drive for better speed.
    verbose: int
        Determines the verbosity.

    Returns
    -------
    (int, int)
        A tuple of the number of fasta records sampled
        and permuted records written.

    """
    index = pickle.load(open(index_dir, 'rb'))
    if not thresholds:
        thresholds = [None] * 3
    try:
        strategy = index['select']['strategy']
    except KeyError:
        strategy = None
    print("The index selection strategy is {}".format(strategy),
          file=sys.stderr)
    if strategy.startswith("GH"):
        print("Getting the testing data with genome holdout.",
              file=sys.stderr)
        test_count = get_sample_worker(taxid, sublevels, index, genomes_dir,
                                       number, length, data_dir,
                                       split=False, split_amount=split_amount,
                                       include_wild=include_wild, prob=prob,
                                       thresholding=thresholding,
                                       chop=chop,
                                       window_length=window_length,
                                       amino_acid=amino_acid,
                                       temp_dir=temp_dir,
                                       include_list=[_TEST],
                                       threshold=thresholds[_TEST-1],
                                       verbose=verbose)
        shutil.move(os.path.join(data_dir, str(taxid), "train"),
                    os.path.join(data_dir, str(taxid), "test"))
        print("Getting the training data with genome holdout.",
              file=sys.stderr)
        train_count = get_sample_worker(taxid, sublevels, index, genomes_dir,
                                        number, length, data_dir,
                                        split=False,
                                        split_amount=split_amount,
                                        include_wild=include_wild, prob=prob,
                                        thresholding=thresholding,
                                        chop=chop,
                                        window_length=window_length,
                                        amino_acid=amino_acid,
                                        temp_dir=temp_dir,
                                        include_list=[_TRAIN],
                                        threshold=thresholds[_TRAIN-1],
                                        verbose=verbose)
        return tuple(map(operator.add, test_count, train_count))
    else:
        return get_sample_worker(taxid, sublevels, index, genomes_dir,
                                 number, length, data_dir,
                                 split=split, split_amount=split_amount,
                                 include_wild=include_wild, prob=prob,
                                 thresholding=thresholding,
                                 chop=chop,
                                 window_length=window_length,
                                 amino_acid=amino_acid, temp_dir=temp_dir,
                                 threshold=thresholds[0],
                                 verbose=verbose)


def get_sample_worker(taxid, sublevels, index, genomes_dir,
                      number, length, data_dir,
                      split=True, split_amount='0.8,0.1,0.1',
                      include_wild=False,
                      prob=_RC_PROB, thresholding=False, chop=False,
                      window_length=50,
                      amino_acid=False, temp_dir="/localscratch/",
                      include_list=[True],
                      threshold=None,
                      verbose=0):
    """
    Get a random sample.  Create training, validation, and testing data sets
    and put them in the appropriate folders.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id
    sublevels: iterable
        A set of taxonomic ids below the taxid to sample from.
    index: dict
        The genomes index object.
    genomes_dir: str
        The location of the root of where the fasta files are stored.
    number: int
        The number of samples to take.
    length: int
        The number of bases for each sample
    data_dir: str
        The path to the data directory where fasta files will be written.
    split: bool
        Determine whether to split the data or not.
    split_amount: str
        A comma seperated list of floats representing the percentage of the
        data to be used for training, validation, and testing data
    include_wild: boolean
        Determines if sequences with wild card characters will be kept.
    prob: float
        A number between 0 and 1.0 that represents the probability that the
        reverse complement will be done.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.
    chop: bool
        If True, chop up the genome into kmers based on a sliding window.
    window_length: int
        The length of the offset for the sliding window
        if thresholding or chopping is used.
    amino_acid: bool
        If True, the data is amino acid data.
    temp_dir: str
        A path to write temporary files to.  This could be on the local hard
        drive for better speed.
    include_list: list
        Determines the value from the genomes index that determines
        if the genome will be included for consideration.  Ordinarily,
        this will be set to [True], but for the genome holdout strategy,
        it could be [1] or [2].
    threshold: int
        A value that controls how much genomic content to include.
    verbose: int
        Determines the verbosity.

    Returns
    -------
    (int, int)
        A tuple of the number of fasta records sampled
        and permuted records written.

    """
    print("Determining accessions to sample from.", file=sys.stderr)
    sys.stderr.flush()
    accession_counts = uniform_samples_at_rank(index, sublevels, genomes_dir,
                                               number, length,
                                               include_wild, amino_acid,
                                               temp_dir, include_list, threshold)
    if not accession_counts:
        print("{} has no sublevels.".format(taxid), file=sys.stderr)
        return (0, 0)
    print("Getting the kmer samples.", file=sys.stderr)
    sys.stderr.flush()
    fasta_path_init = os.path.join(temp_dir, str(taxid) + ".init.fasta")
    taxid_path = os.path.join(temp_dir, str(taxid) + ".taxid")
    fasta_file = open(fasta_path_init, "w")
    taxid_file = open(taxid_path, "w")
    fasta_records_count = get_fasta(accession_counts, length,
                                    index, genomes_dir, fasta_file,
                                    taxid_file, include_wild=include_wild,
                                    window_length=window_length,
                                    temp_dir=temp_dir,
                                    thresholding=thresholding,
                                    chop=chop,
                                    amino_acid=amino_acid,
                                    verbose=verbose)
    fasta_file.close()
    taxid_file.close()
    print("Finished getting the kmer samples.", file=sys.stderr)
    sys.stderr.flush()
    if not amino_acid:
        print("Getting the reverse complements.", file=sys.stderr)
        sys.stderr.flush()
        fasta_path = os.path.join(temp_dir, str(taxid) + ".fasta")
        _, _ = get_rc_fasta(fasta_path_init,
                            fasta_path,
                            prob=prob,
                            remove=False,
                            verbose=verbose)
        if os.path.isfile(fasta_path_init):
            os.remove(fasta_path_init)
    else:
        fasta_path = fasta_path_init
    print("Permuting the fasta records.", file=sys.stderr)
    sys.stderr.flush()
    permute_count = randomly_permute_fasta_taxid(fasta_path,
                                                 taxid_path,
                                                 fasta_path,
                                                 taxid_path,
                                                 split=split,
                                                 split_amount=split_amount)
    print("Writing the fasta file(s) to their final destination.",
          file=sys.stderr)
    sys.stderr.flush()
    for ext, ml_path in [(".train", "train"),
                         (".validate", "validate"),
                         (".test", "test")]:
        save_dir = os.path.join(data_dir, str(taxid), ml_path)
        if not os.path.exists(save_dir):
            try:
                os.makedirs(save_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        if not split:
            ext = ""
        shutil.copy(fasta_path + ext, os.path.join(save_dir,
                                                   str(taxid) + ".fasta"))
        shutil.copy(taxid_path + ext, os.path.join(save_dir,
                                                   str(taxid) + ".taxid"))
        try:
            os.remove(fasta_path + ext)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
        try:
            os.remove(taxid_path + ext)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
        if not split:
            break
    if os.path.isfile(fasta_path):
        os.remove(fasta_path)
    if os.path.isfile(taxid_path):
        os.remove(taxid_path)
    return fasta_records_count, permute_count


def parallel_sample(taxid_list, genomes_dir, ranks, index_dir, number, length,
                    data_dir, split, split_amount, processes,
                    include_wild=False, prob=_RC_PROB,
                    thresholding=False, chop=False,
                    window_length=100,
                    amino_acid=False, 
                    thresholds=None,
                    temp_dir="/localscratch/",
                    verbose=0):
    """
    Get samples of data in parallel and writes them into files and a data
    directory.

    Parameters
    ----------
    taxid_list: list<int>
        A list of ints representing taxonomic ids.
    genomes_dir: str
        The location of the root of where the fasta files are stored.
    ranks: dict
        The tree object giving taxonomic ids at every rank.
    index_dir: str
        The path to the genomes index object.
    number: int
        The number of samples to take.
    length: int
        The number of bases for each sample.
    data_dir: str
        The path to the data directory where fasta files will be written.
    split: bool
        Determine whether to split the data or not.
    split_amount: str
        A comma seperated list of floating point values that represent how to
        split the training, validation, and test data sets.
    processes: int
        The number of processes to use.  This should be at or less than the
        number of physical cores that the CPU has.
    include_wild: bool
        When true, samples will include wild card characters.
        When false, samples will not include wild card characters.
        prob: float
        A number between 0 and 1.0 that represents the probability that the
        reverse complement will be done.
    prob: float
        A number between 0 and 1.0 that represents the probability that the
        reverse complement will be done.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.
    chop: bool
        If True, chop up the genome into kmers based on a sliding window.
    window_length: int
        The length of the offset for the sliding window
        if thresholding or chopping is used.
    amino_acid: bool
        If True, the data is amino acid data.
    thresholds: list<int>
        Values that control how much genomic content to include.  
        Multiple values are only valid for genome holdout strategies.
    temp_dir: str
        A path to write temporary files to.  This could be on the local hard
        drive for better speed.
    verbose: int
        Determines the verbosity.


    Returns
    -------
    list<(int, int)>
        A list of the number of fasta records written by each process
        in the same order as the taxid_list.

    """
    with Pool(processes=processes) as pool:
        process_list = []
        for taxid in taxid_list:
            sublevels = ranks[taxid]
            process_list.append(pool.apply_async(get_sample,
                                                 args=(taxid,
                                                       sublevels, index_dir,
                                                       genomes_dir,
                                                       number,
                                                       length, data_dir,
                                                       split,
                                                       split_amount,
                                                       include_wild,
                                                       prob,
                                                       thresholding,
                                                       chop,
                                                       window_length,
                                                       amino_acid,
                                                       thresholds,
                                                       temp_dir,
                                                       verbose)))
        output = []
        for taxid, process_desc in zip(taxid_list, process_list):
            counts = process_desc.get()
            output.append(counts)
            print("{}: {} samples drawn, {} samples written".
                  format(taxid, counts[0], counts[1]), file=sys.stderr)
    return output
