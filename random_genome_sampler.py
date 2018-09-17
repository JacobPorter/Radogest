#!/usr/bin/env python
"""
Produces a fasta file of randomly chosen sequences corresponding to the
given taxonomic id.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import datetime
import pickle
import sys
import random
import string
import subprocess
import os
import sqlite3
from collections import defaultdict
from get_paths import get_paths
from split_genomes import split_genomes
from ete3 import NCBITaxa
import SeqIterator

ncbi = NCBITaxa()

# The location of the genomes and the tools that need to be used.
# GENOMES_LOCATION = "/groups/fungcat/datasets/current/fasta/Genomes/"
# GENOMES_LOCATION = "/scratch/fungcat/jsporter/Genomes/"
# BEDTOOLS_DIR = "/home/jsporter/Applications/bedtools2/bin/"
# UCSC_DIR = "/home/jsporter/Applications/UCSC/"
# If the sample falls below NUMBER_CUTOFF, get NUMBER_SAMPLE samples instead.
NUMBER_CUTOFF = 150
NUMBER_SAMPLE = 300
# The actual number of random samples taken is multipled by this multiplier.
# This is done so that samples with N's in them can be excluded.
SAMPLE_MULTIPLIER = 1.2
# If there are more than GENOMES_TO_KEEP genomes, randomly shuffle the genomes
# and keep only GENOMES_TO_KEEP
# GENOMES_TO_KEEP = 50

# Tools and genomes locations.
path_dict = get_paths()
GENOMES_LOCATION = path_dict['GENOMES']
BEDTOOLS_DIR = path_dict['BEDTOOLS']
UCSC_DIR = path_dict['UCSC']


def binary_search(target, array, begin, end):
    """
    Binary search

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


def uniform_samples_at_rank(taxid, index, sublevels, number):
    """
    Get a count of samples from each genome under the taxids given by ranks.

    Parameters
    ----------
    taxid: int
        A taxonomic id number to query from
    index: dict
        A dictionary representing genome locations and genome taxid
        associations
    sublevels: iterable
        An iterable that gives ranks underneath taxid to sample from.
    number: int
        The number of random sequences to sample.

    Returns
    -------
    dict
        A dictionary keyed on taxids where the value is the number of sequences
        to sample for each taxid.

    """
    try:
        uniform_number = round(number / len(sublevels))
    except ZeroDivisionError:
        return None
    uniform_sample_counts = []
    for level in sublevels:
        uniform_sample_counts.append((level, uniform_samples(level, index,
                                                             uniform_number)))
    return uniform_sample_counts


def uniform_samples(taxid, index, number):
    """
    Get a count of samples from each genome under taxid.  Does not do any
    sampling.

    Parameters
    ----------
    taxid: int
        A taxonomic id number to query from
    index: dict
        A dictionary representing genome locations and genome taxid
        associations
    number: int
        The number of random sequences to sample.

    Returns
    -------
    dict
        A dictionary keyed on accession ids.  The value is the number of
        sequences to sample from the corresponding genome.

    """
    accessions = [accession for
                  accession in index['taxids'][taxid] if
                  index['taxids'][taxid][accession]]
    # if len(accessions) > GENOMES_TO_KEEP:
    #     random.shuffle(accessions)
    #     accessions = accessions[:GENOMES_TO_KEEP]
    genome_lengths = [index['genomes'][accession]['base_length']
                      for accession in
                      accessions]
    try:
        name = ncbi.get_taxid_translator([taxid])[taxid]
    except KeyError:
        name = ""
    except sqlite3.DatabaseError:
        name = ""
    sys.stderr.write("For {}:{}, sampling from "
                     "{} genomes.\n".format(taxid, name, len(accessions)))
    sys.stderr.flush()
    for i, _ in enumerate(genome_lengths):
        if i == 0:
            continue
        genome_lengths[i] += genome_lengths[i-1]
    try:
        total_length = genome_lengths[-1]
    except IndexError as ie:
        print("Genome lengths are empty for taxid: {}".format(taxid),
              file=open("/home/jsporter/Sampling_Out/index_error_"+
                        str(taxid)+ ".out", 'w'))
        raise ie
    accession_counts = defaultdict(int)
    for _ in range(number):
        i = binary_search(random.randint(0, total_length - 1),
                          genome_lengths, 0, len(genome_lengths) - 1)
        accession_counts[accessions[i]] += 1
    return accession_counts


def get_fasta(accession_counts_list, length, index,
              output, taxid_file, window_length=50, verbose=False,
              thresholding=False, temp_dir='/tmp/'):
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
    output: writable
        A writable object to store fasta files too.  This could be an open
        file.
    taxid_file: writable
        A file-like object to store taxonomic ids that correspond to each
        sampled sequence.  Each taxonomic id will be on its own line.
    tmp: str
        A path to a directory to store temporary files.
    verbose: bool
        If True, print messages.
    thresholding: bool
        If True, use the whole genome if the samples requested is larger
        than the genome.  If False, and the genome is smaller than the
        samples requested, portions of the genome will be over represented
        because of random sampling.


    Returns
    -------
    int
        A count of the number of fasta records written.

    """
    final_file = SeqIterator.SeqWriter(output, file_type='fasta')
    fasta_record_count = 0
    for taxid, accession_counts in accession_counts_list:
        for accession in accession_counts.keys():
            if accession_counts[accession] <= 0:
                continue

            rand_string = "".join(random.choices(
                string.ascii_letters + string.digits, k=8))
            my_fasta = os.path.join(temp_dir,
                                    accession + "_" +
                                    rand_string +  "_random.fa")
            if verbose:
                sys.stderr.write("Writing fasta records for {}:\t".
                                 format(accession))
            # A thresholding feature.  If the genome is too small, use the
            # whole genome.
            if thresholding and accession_counts[accession] > float(
                index['genomes'][accession]['base_length']) / length:
                records_written = split_genomes([accession], length,
                                                index, final_file)
                for _ in range(records_written):
                    taxid_file.write(str(taxid) + "\n")
                fasta_record_count += records_written
                continue
            accession_location = os.path.join(GENOMES_LOCATION +
                                              index['genomes']
                                              [accession]
                                              ['location'])
            onlyfiles = [f for f in os.listdir(accession_location) if
                         os.path.isfile(os.path.join(accession_location, f))]
            for f in onlyfiles:
                if (f.endswith(".fna") or
                        f.endswith(".fasta") or
                        f.endswith(".fa")):
                    fasta_location = os.path.join(accession_location, f)
                elif f.endswith(".2bit"):
                    twobit_location = os.path.join(accession_location, f)
                elif f.endswith(".fai"):
                    fai_location = os.path.join(accession_location, f)
            bedtools_file = os.path.join(temp_dir, accession + "_" +
                                         rand_string + "_random.bed")
            get_random = True
            accession_cnt = 0
            while get_random:
                accession_number = accession_counts[accession] - accession_cnt
                bed_2bit_counts = get_random_bed_fast(accession_number,
                                                      length,
                                                      taxid,
                                                      accession,
                                                      fai_location,
                                                      bedtools_file,
                                                      fasta_location,
                                                      my_fasta,
                                                      taxid_file,
                                                      final_file)
                get_random = not bed_2bit_counts[0]
                accession_cnt += bed_2bit_counts[1]
                if verbose:
                    sys.stderr.write(str(bed_2bit_counts[1]) +
                                     "  " +
                                     str(bed_2bit_counts[2] /
                                         (bed_2bit_counts[1] +
                                          bed_2bit_counts[2])) +
                                     " ")
                    sys.stderr.flush()
            fasta_record_count += accession_cnt
            if verbose:
                sys.stderr.write("\n")
                sys.stderr.flush()
    return fasta_record_count


def get_random_bed_2bit(number, length, taxid, accession, fai_location,
                        bedtools_file, twobit_location, my_fasta,
                        taxid_file, final_file):
    """
    Get random nucleotide sequences from a bed file and a 2bit file.  Exclude
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
    bedtoods_file: str
        The location to store the bedtools output file too.
    twobit_location: str
        The location of the 2bit genome file
    my_fasta: str
        The location to store a temporary fasta file for random sampling.
    taxid_file: writable
        A writable object where taxids are written for each random sample.
    final_file: writable
        A writable object that stores the fasta records.  This file represents
        the end product of all of the random sampling.

    Returns
    -------
    (boolean, int, int)
        A boolean to indicated if the records written is equal to number.
        An integer counting the records written,
        and an integer of the records that had indeterminate N characters.

    """
    taxid = str(taxid)
    if number <= NUMBER_CUTOFF:
        my_sample = NUMBER_SAMPLE
    else:
        my_sample = int(number * SAMPLE_MULTIPLIER)
    bedtools_fd = open(bedtools_file, 'w')
    subprocess.run([BEDTOOLS_DIR + "bedtools", "random", "-l",
                    str(length), "-n",
                    str(my_sample), "-g",
                    fai_location], stdout=bedtools_fd)
    subprocess.run([UCSC_DIR + "twoBitToFa", "-bed="+bedtools_file,
                    "-bedPos", twobit_location, my_fasta])
    intermediate_fasta_file = SeqIterator.SeqIterator(my_fasta,
                                                      file_type='fasta')
    records_with_n = 0
    records_written = 0
    records_malformed = 0
    for fasta_record in intermediate_fasta_file:
        if records_written >= number:
            break
        record_id, record_seq = fasta_record
        record_seq = record_seq.upper()
        if ">" in record_seq:
            records_malformed += 1
            continue
        if "N" in record_seq:
            records_with_n += 1
            continue
        record_id = accession + ":" + taxid + ":" + record_id
        final_file.write((record_id, record_seq))
        taxid_file.write(taxid + "\n")
        records_written += 1
    bedtools_fd.close()
    os.remove(bedtools_file)
    intermediate_fasta_file.close()
    os.remove(my_fasta)
    return (records_written >= number, records_written, records_with_n)


def get_random_bed_fast(number, length, taxid, accession, fai_location,
                        bedtools_file, fasta_location, my_fasta,
                        taxid_file, final_file):
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
    bedtoods_file: str
        The location to store the bedtools output file too.
    fasta_location: str
        The location of the fasta genome file
    my_fasta: str
        The location to store a temporary fasta file for random sampling.
    taxid_file: writable
        A writable object where taxids are written for each random sample.
    final_file: writable
        A writable object that stores the fasta records.  This file represents
        the end product of all of the random sampling.

    Returns
    -------
    (boolean, int, int)
        A boolean to indicated if the records written is equal to number.
        An integer counting the records written,
        and an integer of the records that had indeterminate N characters.

    """
    taxid = str(taxid)
    if number <= NUMBER_CUTOFF:
        my_sample = NUMBER_SAMPLE
    else:
        my_sample = int(number * SAMPLE_MULTIPLIER)
    bedtools_fd = open(bedtools_file, 'w')
    subprocess.run([BEDTOOLS_DIR + "bedtools", "random", "-l",
                    str(length), "-n",
                    str(my_sample), "-g",
                    fai_location], stdout=bedtools_fd)
    subprocess.run([BEDTOOLS_DIR + "bedtools", "getfasta", "-fi",
                    fasta_location, "-bed", bedtools_file, "-fo", my_fasta])
    intermediate_fasta_file = SeqIterator.SeqIterator(my_fasta,
                                                      file_type='fasta')
    records_with_n = 0
    records_written = 0
    for fasta_record in intermediate_fasta_file:
        if records_written >= number:
            break
        record_id, record_seq = fasta_record
        record_seq = record_seq.upper()
        if "N" in record_seq:
            records_with_n += 1
            continue
        record_id = accession + ":" + taxid + ":" + record_id
        final_file.write((record_id, record_seq))
        taxid_file.write(taxid + "\n")
        records_written += 1
    bedtools_fd.close()
    os.remove(bedtools_file)
    intermediate_fasta_file.close()
    os.remove(my_fasta)
    return (records_written >= number, records_written, records_with_n)


def main():
    """Parse arguments"""
    parser = argparse.ArgumentParser(description=('Produces a fasta file of '
                                                  'randomly sampled sequences '
                                                  'given a taxonomic id.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("taxid", type=int,
                        help=("A taxonomic id to sample from."))
    parser.add_argument("--tree", "-r", type=str,
                        help=("The pickled ranks tree object."),
                        default="./tree.pck")
    parser.add_argument("--index", "-i", type=str,
                        help=("The location of the index that "
                              "maps taxonomic ids to genomes and "
                              "their locations, etc."),
                        default='./genomes.withcounts.selected.index.pck')
    parser.add_argument("--number", "-n", type=int,
                        help=("The number of sequences to sample."),
                        default=3200000)
    parser.add_argument("--length", "-l", type=int,
                        help=("The length of the sequences to sample."),
                        default=100)
    parser.add_argument("--thresholding", "-g", action='store_true',
                        default=False)
    parser.add_argument("--window_length", "-w", type=int,
                        help=("The length of the window size to use when "
                              "using thresholding."),
                        default=50)
    parser.add_argument("--output", "-o", type=str,
                        help=("The location to store the output fasta file."),
                        default="./random_sample.fasta")
    parser.add_argument("--tax_id_file", "-f", type=str,
                        help=("The location of the file to store the taxid "
                              "associated with each sequence."),
                        default="./random_sample.taxid")
    parser.add_argument("--temp_dir", "-t", type=str,
                        help=("The location of a directory to store "
                              "temporary files."),
                        default="/localscratch/")
    parser.add_argument("--verbose", "-v", action='store_true', default=False)
    args = parser.parse_args()
    now = datetime.datetime.now()
    sys.stderr.write("Generating {} random samples each of length {} for "
                     "{}.\n".format(args.number,
                                    args.length,
                                    args.taxid))
    sys.stderr.write("The output will be written to {}.\n".format(args.output))
    sys.stderr.write("Using index at {}.  Using tree ranks object "
                     "at {} and temp "
                     "directory at {}.\n".format(args.index,
                                                 args.tree,
                                                 args.temp))
    index = pickle.load(open(args.index, 'rb'))
    ranks = pickle.load(open(args.tree, 'rb'))
    sublevels = ranks[args.taxid]
    accession_counts = uniform_samples_at_rank(args.taxid, index, sublevels,
                                               args.number)
    if isinstance(args.tax_id_file, str):
        args.tax_id_file = open(args.tax_id_file, 'w')
    fasta_records_count = get_fasta(accession_counts, int(args.length),
                                    index, open(args.output, 'w'),
                                    args.tax_id_file, args.window_length,
                                    args.verbose, args.thresholding,
                                    args.temp_dir)
    sys.stderr.write("There were {} fasta records written.\n"
                     "The process took {} time.\n".format(
                         fasta_records_count, datetime.datetime.now() - now))


if __name__ == '__main__':
    main()
