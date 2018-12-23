"""
Produces a fasta file of randomly chosen sequences corresponding to the
given taxonomic id.
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
from multiprocessing import Pool
from collections import defaultdict

from ete3 import NCBITaxa

from SeqIterator.SeqIterator import SeqReader, SeqWriter
from library.permute import randomly_permute_fasta_taxid
from library.chop import split_genomes

ncbi = NCBITaxa()

# If the sample falls below NUMBER_CUTOFF, get NUMBER_SAMPLE samples instead.
NUMBER_CUTOFF = 150
NUMBER_SAMPLE = 300
# The actual number of random samples taken is multipled by this multiplier.
# This is done so that samples with N's in them can be excluded.
SAMPLE_MULTIPLIER = 1.2

from config import BEDTOOLS
# from config import GENOMES_NT, GENOMES_AA, GENOMES_CD


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


def get_fasta(accession_counts_list, length, index, genomes_dir,
              output, taxid_file, window_length=50, verbose=False,
              thresholding=False, amino_acid=False, 
              temp_dir='/localscratch/'):
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
    final_file = SeqWriter(output, file_type='fasta')
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
                                                index, genomes_dir, 
                                                final_file, 
                                                window_length=window_length)
                for _ in range(records_written):
                    taxid_file.write(str(taxid) + "\n")
                fasta_record_count += records_written
                continue
            accession_location = os.path.join(genomes_dir +
                                              index['genomes']
                                              [accession]
                                              ['location'])
            onlyfiles = [f for f in os.listdir(accession_location) if
                         os.path.isfile(os.path.join(accession_location, f))]
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
                                                      final_file,
                                                      amino_acid)
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


def get_random_bed_fast(number, length, taxid, accession, fai_location,
                        bedtools_file, fasta_location, my_fasta,
                        taxid_file, final_file, amino_acid=False):
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
    subprocess.run([BEDTOOLS + "bedtools", "random", "-l",
                    str(length), "-n",
                    str(my_sample), "-g",
                    fai_location], stdout=bedtools_fd)
    subprocess.run([BEDTOOLS + "bedtools", "getfasta", "-fi",
                    fasta_location, "-bed", bedtools_file, "-fo", my_fasta])
    intermediate_fasta_file = SeqReader(my_fasta, file_type='fasta')
    records_with_n = 0
    records_written = 0
    for fasta_record in intermediate_fasta_file:
        if records_written >= number:
            break
        record_id, record_seq = fasta_record
        record_seq = record_seq.upper()
        if not amino_acid and "N" in record_seq:
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



"""
Runs multiple taxonomic id sampling instances in parallel.
Creates training, validation,
and testing data and puts them in directories that Plinko expects.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


def get_sample(taxid, sublevels, index_dir, genomes_dir,
               number, length, data_dir,
               split=True,
               split_amount='0.8,0.1,0.1', thresholding=False,
               window_length=50, amino_acid=False, temp_dir="/localscratch/"):
    """
    Get a random sample.  Create training, validation, and testing data sets
    and put them in the appropriate folders.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id
    rank: str
        A rank in the taxonomic system.  Examples: genus, family, etc.
    sublevels: iterable
        A set of taxonomic ids below the taxid to sample from.
    index_dir: str
        A path to the pickled genomes index object.
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
    tmp_dir: str
        A path to write temporary files to.  This could be on the local hard
        drive for better speed.

    Returns
    -------
    (int, int)
        A tuple of fasta records sampled and permuted records written.

    """
    index = pickle.load(open(index_dir, 'rb'))
    accession_counts = uniform_samples_at_rank(taxid, index, sublevels, number)
    if not accession_counts:
        print("{} has no sublevels.".format(taxid), file=sys.stderr)
        return (0, 0)
    fasta_path = os.path.join(temp_dir, str(taxid) + ".fasta")
    taxid_path = os.path.join(temp_dir, str(taxid) + ".taxid")
    fasta_file = open(fasta_path, "w")
    taxid_file = open(taxid_path, "w")
    fasta_records_count = get_fasta(accession_counts, length,
                                    index, genomes_dir, fasta_file,
                                    taxid_file, window_length=window_length,
                                    temp_dir=temp_dir,
                                    thresholding=thresholding,
                                    amino_acid=amino_acid)
    fasta_file.close()
    taxid_file.close()
    permute_count = randomly_permute_fasta_taxid(fasta_path,
                                                 taxid_path,
                                                 fasta_path,
                                                 taxid_path,
                                                 split=split,
                                                 split_amount=split_amount)
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


def create_directories(data_dir):
    """
    Create directories to put the data if they do not exist.

    Parameters
    ----------
    data_dir: str
        A path to the top level of the directory to put the data.

    Returns
    -------
    None

    """
    # ranks = ["superkingdom", "kingdom", "phylum",
    #         "class", "order", "family", "genus"]
    # data_sets = ["train", "test", "validate"]
    # for rank in ranks:
    # for data_set in data_sets:
    #     path = os.path.join(data_dir, rank, data_set)
    #     if not os.path.exists(path):
    #         try:
    #             os.makedirs(path)
    #         except OSError as e:
    #             if e.errno != errno.EEXIST:
    #                 raise


def parallel_sample(taxid_list, genomes_dir, ranks, index_dir, number, length,
                    data_dir, split, split_amount, processes,
                    thresholding=False, window_length=100, amino_acid=False,
                    temp_dir="/tmp"):
    """
    Get samples of data in parallel and writes them into files and a data
    directory that Plinko expects.

    Parameters
    ----------
    taxid_list: list<int>
        A list of ints representing taxonomic ids.
    genomes_dir: str
        The location of the root of where the fasta files are stored.
    ranks: dict
        The ranks object giving taxonomic ids at every rank.
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

    Returns
    -------
    list<(int, int)>
        A list of fasta records written by each process in the same order as
        the taxid_list.

    """
    # create_directories(data_dir)
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
                                                       thresholding,
                                                       window_length,
                                                       amino_acid,
                                                       temp_dir)))
        output = []
        for taxid, process_desc in zip(taxid_list, process_list):
            counts = process_desc.get()
            output.append(counts)
            print("{}: {} samples drawn, {} samples written".
                  format(taxid, counts[0], counts[1]), file=sys.stderr)
        return output
