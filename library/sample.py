"""
Produces a fasta file of randomly chosen kmer sequences corresponding to the
given taxonomic id(s).

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import os
import errno
import shutil
import sys
import random
import string
import subprocess
import sqlite3
import tempfile
import operator
from multiprocessing import Pool, Value, Process, Manager
# from multiprocessing import Pipe
from collections import defaultdict

from ete3 import NCBITaxa

from SeqIterator.SeqIterator import SeqReader, SeqWriter
from library.permute import randomly_permute_fasta_taxid
from library.chop import chop_genomes
from library.util import which

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

# Length of the random string
RAND_LEN = 24

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


def file_locations(accession, accession_location, temp_dir):
    """
    Get file locations for sampling from the fasta file.

    Parameters
    ----------
    accession: str
        The accession id of the genome.
    genomes_dir: str
        The location where the genomes or fasta files are stored.
    accession_location: str
        The path to the accessions location.
    temp_dir: str
        The temporary directory to store files

    Returns
    -------
    dict
        A dictionary of file locations.

    """
    rand_string = "".join(random.choices(
        string.ascii_letters + string.digits, k=RAND_LEN))
    my_fasta = os.path.join(temp_dir,
                            accession + "_" +
                            rand_string + "_random.fa")
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


# def index_service(index_dir, pipes):
#     """
#     Access the genomes index object through pipes.
#     This allows for paralellism.
#
#     Parameters
#     ----------
#     index_dir: str
#         The location of the genomes index object.
#     pipes: list<Pipe>
#         A list of Pipe objects to send and receive
#         requests for information from the index.
#
#     Returns
#     -------
#     None
#
#     """
#     index = pickle.load(open(index_dir, 'rb'))
#     stay = True
#     while(stay):
#         ready_pipes = connection.wait(pipes)
#         for pipe in ready_pipes:
#             index_accessor = pipe.recv()
#             if index_accessor:
#                 d1 = index
#                 for accessor in index_accessor:
#                     d1 = d1[accessor]
#                 pipe.send(d1)
#             else:
#                 stay = False


def file_service(fasta_file_location, taxid_file_location,
                 record_count, queue, pills=1):
    """
    Write to a fasta file and a taxid file
    by pulling information off of a queue.

    Parameters
    ----------
    fasta_file_location: str
        The location of a fasta file to write to.
    taxid_file_location: str
        The location of a taxid file to write to.
    record_count: Value("i")
        A shared integer value for recording the records processed.
    queue: Queue
        A multiprocessing queue to read information from.
    pills: int
        The number of poison pills to swallow before quitting.
        This number should be the same number as the number of
        processes that are using this service.

    Returns
    -------
    count: int
        A count of the number of records processed.

    """
    fasta_file = SeqWriter(open(fasta_file_location, "w"), file_type='fasta')
    taxid_file = open(taxid_file_location, "w")
    swallowed = 0
    while(swallowed < pills):
        record = queue.get()
        if record:
            fasta_file.write(record[0:2])
            taxid = str(record[2])
            taxid_file.write(taxid + "\n")
            if taxid not in record_count:
                record_count[taxid] = 1
            else:
                temp_count = record_count[taxid]
                record_count[taxid] = temp_count + 1
        else:
            swallowed += 1
    fasta_file.flush()
    taxid_file.flush()
    fasta_file.close()
    taxid_file.close()


def get_fasta(accession_counts_list, length, index, genomes_dir,
              fasta_path, taxid_path,
              index_dir,
              include_wild=False,
              window_length=50,
              thresholding=False, chop=False,
              amino_acid=False,
              temp_dir='/localscratch/',
              processes=1,
              verbose=0):
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
    fasta_path: str
        A path for a fasta file.
    taxid_path: str
        A path for a taxid file.
    index_dir: str
        The path to the genomes index object.
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
    processes: int
        The number of processes to use.
    verbose: bool
        If True, print messages.

    Returns
    -------
    int
        A count of the number of fasta records written.

    """
    if processes == 1:
        final_file = SeqWriter(open(fasta_path, "w"), file_type='fasta')
        taxid_file = open(taxid_path, "w")
        fasta_record_count = {}
    else:
        pool = Pool(processes=processes)
#         pool_connections = []
#         service_connections = []
#         for _ in len(processes):
#             c1, c2 = Pipe()
#             pool_connections.append(c1)
#             service_connections.append(c2)
#         index_process = Process(target=index_service,
#                                 args=(index_dir, service_connections))
        # fasta_record_count = Value('i', 0)
        manager = Manager()
        fasta_record_count = manager.dict()
        queue = manager.Queue()
        file_process = Process(target=file_service,
                               args=(fasta_path, taxid_path,
                                     fasta_record_count,
                                     queue))
#         index_process.start()
        file_process.start()
    for taxid, accession_counts in accession_counts_list:
        my_fasta_record_count = 0
        for accession in accession_counts.keys():
            if accession_counts[accession] <= 0:
                continue
            if verbose > 1:
                sys.stderr.write("Writing fasta records for "
                                 "taxid {} from genome accession {}\n".
                                 format(taxid, accession))
            # A thresholding and chopping feature.
            # If thresholding and the genome is too small,
            # use the whole genome.
            # Or chop the genome up.
            if (thresholding and accession_counts[accession] > float(
                    index['genomes'][
                        accession]['contig_sum']) / length) or chop:
                location = genomes_dir + index['genomes'][accession][
                        'location']
                if processes == 1:
                    records_written = chop_genomes([accession],
                                                   length,
                                                   [location],
                                                   taxid,
                                                   final_file,
                                                   queue=None,
                                                   include_wild=include_wild,
                                                   window_length=window_length)
                    for _ in range(records_written):
                        taxid_file.write(str(taxid) + "\n")
                    my_fasta_record_count += records_written
                else:
                    pool.apply_async(chop_genomes,
                                     args=([accession],
                                           length,
                                           [location],
                                           taxid,
                                           None,
                                           queue,
                                           include_wild,
                                           window_length))
                continue
            accession_location = os.path.join(genomes_dir +
                                              index['genomes']
                                              [accession]
                                              ['location'])
            file_locations_d = file_locations(accession,
                                              accession_location,
                                              temp_dir)
            fai_location = file_locations_d["fai"]
            fasta_location = file_locations_d["fasta_location"]
            if processes == 1:
                my_fasta_record_count += random_bed_fast_worker(
                    accession_counts[accession],
                    length,
                    taxid,
                    accession,
                    fai_location,
                    fasta_location,
                    taxid_file,
                    final_file,
                    None,
                    include_wild,
                    amino_acid,
                    temp_dir,
                    verbose)
            else:
                pool.apply_async(
                    random_bed_fast_worker,
                    args=(accession_counts[accession],
                          length,
                          taxid,
                          accession,
                          fai_location,
                          fasta_location,
                          None,
                          None,
                          queue,
                          include_wild,
                          amino_acid,
                          temp_dir,
                          verbose))
        fasta_record_count[str(taxid)] = my_fasta_record_count
    if processes == 1:
        final_file.close()
        taxid_file.close()
    else:
        pool.close()
        pool.join()
#         index_process.join()
        queue.put(None)
        file_process.join()
        # index_process.close()
        # file_process.close()
    return fasta_record_count


def random_bed_fast_worker(total_accession_count,
                           length,
                           taxid,
                           accession,
                           fai_location,
                           fasta_location,
                           taxid_file,
                           final_file,
                           queue,
                           include_wild,
                           amino_acid,
                           temp_dir,
                           verbose=0):
    """
    Get random nucleotide sequences from a bed file and a fasta file.  Exclude
    sequences with N's in them.
    Iterate until there are total_accession_count records.

    Parameters
    ----------
    total_accession_count: int
        The total number of sequences to draw.
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
    queue: Queue
        If there is a queue object, use it to write to the taxid file
        and the fasta file.  This is for parallelism.
    include_wild: boolean
        Determines if sequences with wild card characters will be kept.
    amino_acid: bool
        If True, the data is amino acid data.
    temp_dir: str
        A path to a directory to store temporary files.
    verbose: bool
        If True, print messages.

    Returns
    -------
    accession_cnt: int
        The count of the number of samples drawn.

    """
    get_random = True
    accession_cnt = 0
    while get_random:
        accession_number = total_accession_count - accession_cnt
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
                                              temp_dir,
                                              queue)
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
    if verbose:
        sys.stderr.write("\n")
        sys.stderr.flush()
    return accession_cnt


def get_random_bed_fast(number, length, taxid, accession, fai_location,
                        fasta_location, taxid_file, final_file,
                        include_wild=False, amino_acid=False,
                        temp_dir="/localscratch/",
                        queue=None):
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
    queue: Queue
        If there is a queue object, use it to write to the taxid file
        and the fasta file.  This is for parallelism.

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
        bedtools_path_in = BEDTOOLS + "bedtools"
        bedtools_path_out = which(bedtools_path_in)
        if not bedtools_path_out:
            raise FileNotFoundError(bedtools_path_in)
        subprocess.run([bedtools_path_out, "random", "-l",
                        str(length), "-n",
                        str(my_sample), "-g",
                        fai_location], stdout=bedtools_fd)
        subprocess.run([bedtools_path_out, "getfasta", "-fi",
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
            if queue:
                queue.put((record_id, record_seq, taxid))
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
               processes=1,
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
    index_dir: str
        The path to the genomes index object.
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
    processes: int
        The number of processes to use.
    verbose: int
        Determines the verbosity.

    Returns
    -------
    (int, int)
        A tuple of the number of fasta records sampled
        and permuted records written.

    """
    from radogest import read_ds
    index = read_ds(index_dir)
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
                                       index_dir,
                                       split=False, split_amount=split_amount,
                                       include_wild=include_wild, prob=prob,
                                       thresholding=thresholding,
                                       chop=chop,
                                       window_length=window_length,
                                       amino_acid=amino_acid,
                                       temp_dir=temp_dir,
                                       include_list=[_TEST],
                                       threshold=thresholds[_TEST-1],
                                       processes=processes,
                                       verbose=verbose)
        shutil.move(os.path.join(data_dir, str(taxid), "train"),
                    os.path.join(data_dir, str(taxid), "test"))
        print("Getting the training data with genome holdout.",
              file=sys.stderr)
        train_count = get_sample_worker(taxid, sublevels, index, genomes_dir,
                                        number, length, data_dir,
                                        index_dir,
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
                                        processes=processes,
                                        verbose=verbose)
        return ((test_count[0], train_count[0]), 
                test_count[1] + train_count[1])
    else:
        return (get_sample_worker(taxid, sublevels, index, genomes_dir,
                                  number, length, data_dir,
                                  index_dir,
                                  split=split, split_amount=split_amount,
                                  include_wild=include_wild, prob=prob,
                                  thresholding=thresholding,
                                  chop=chop,
                                  window_length=window_length,
                                  amino_acid=amino_acid, temp_dir=temp_dir,
                                  threshold=thresholds[0],
                                  processes=processes,
                                  verbose=verbose), )


def get_sample_worker(taxid, sublevels, index, genomes_dir,
                      number, length, data_dir,
                      index_dir,
                      split=True, split_amount='0.8,0.1,0.1',
                      include_wild=False,
                      prob=_RC_PROB, thresholding=False, chop=False,
                      window_length=50,
                      amino_acid=False, temp_dir="/localscratch/",
                      include_list=[True],
                      threshold=None,
                      processes=1,
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
    index_dir: str
        The path to the genomes index object.
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
    processes: int
        The number of processes to use.
    verbose: int
        Determines the verbosity.

    Returns
    -------
    (int, int)
        A tuple of the number of fasta records sampled
        and permuted records written.

    """
    def get_random_string():
        return ''.join(random.choice(
            string.ascii_uppercase + string.digits + string.ascii_lowercase)
            for _ in range(RAND_LEN))
    print("Determining accessions to sample from.", file=sys.stderr)
    sys.stderr.flush()
    accession_counts = uniform_samples_at_rank(index, sublevels, genomes_dir,
                                               number, length,
                                               include_wild, amino_acid,
                                               temp_dir, include_list,
                                               threshold)
    if not accession_counts:
        print("{} has no sublevels.".format(taxid), file=sys.stderr)
        return (0, 0)
    print("Getting the kmer samples.", file=sys.stderr)
    sys.stderr.flush()
    random_str = get_random_string()
    fasta_path_init = os.path.join(temp_dir,
                                   str(taxid) + "." + random_str +
                                   ".init.fasta")
    taxid_path = os.path.join(temp_dir, str(taxid) + "." + random_str +
                              ".taxid")
    fasta_records_count = get_fasta(accession_counts, length,
                                    index, genomes_dir, fasta_path_init,
                                    taxid_path,
                                    index_dir,
                                    include_wild=include_wild,
                                    window_length=window_length,
                                    temp_dir=temp_dir,
                                    thresholding=thresholding,
                                    chop=chop,
                                    amino_acid=amino_acid,
                                    processes=processes,
                                    verbose=verbose)
    print("Finished getting the kmer samples.", file=sys.stderr)
    sys.stderr.flush()
    if not amino_acid:
        print("Getting the reverse complements.", file=sys.stderr)
        sys.stderr.flush()
        random_str = get_random_string()
        fasta_path = os.path.join(temp_dir, str(taxid) + "." +
                                  random_str + ".fasta")
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
        The number of processes to use.
        If the taxid_list has only one item, then all of the processes
        will be used to do sampling for that one taxid.
        If the taxid_list has multiple items, then one process
        will be used to do sampling per taxid.
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
    def count_string(counts):
        str1 = "\t"
        for d in counts:
            str1 += "{"
            d = dict(d)
            for key in d:
                str1 += str(key) + ":" + str(d[key]) + ", "
            str1 += "}\t"
        return str1
    def printer(taxid, counts):
        print("For taxid {}, drawn {} \n\t and written: {}.".
              format(taxid, count_string(counts[0]), counts[1]), file=sys.stderr)
    output = []
    if len(taxid_list) == 1:
        taxid = taxid_list[0]
        sublevels = ranks[taxid]
        counts = get_sample(taxid,
                            sublevels,
                            index_dir,
                            genomes_dir,
                            number,
                            length, data_dir, split,
                            split_amount,
                            include_wild,
                            prob,
                            thresholding,
                            chop,
                            window_length,
                            amino_acid,
                            thresholds,
                            temp_dir,
                            processes,
                            verbose)
        output.append(counts)
        printer(taxid, counts)
    else:
        with Pool(processes=processes) as pool:
            process_list = []
            for taxid in taxid_list:
                sublevels = ranks[taxid]
                process_list.append(pool.apply_async(get_sample,
                                                     args=(taxid,
                                                           sublevels,
                                                           index_dir,
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
                                                           1,
                                                           verbose)))
            for taxid, process_desc in zip(taxid_list, process_list):
                counts = process_desc.get()
                output.append(counts)
                printer(taxid, counts)
    return output
