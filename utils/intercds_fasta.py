#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create a fasta file of intergenic regions given a gff file and a
whole genome.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>
"""
from __future__ import print_function

import argparse
import datetime
import os
import sys
import operator
import subprocess
import mmap
import math
import re
from collections import defaultdict
from multiprocessing import Pool

try:
    import gffutils
except ImportError:
    print("gffutils could not be imported.  It may require Python 2.7.",
          file=sys.stderr)

from SeqIterator import SeqReader, SeqWriter

# GFF regions to exclude when creating fasta file.
EXCLUDED_REGIONS = ('gene', 'exon')

FASTA_ENDINGS = ['fasta', 'fa', 'fna']
GFF_ENDINGS = ['gff']
FAI_ENDINGS = ['fai']


def name_ends(name, endings, addition=""):
    """
    Determine if a string ends with something.

    Parameters
    ----------
    name: str
        The string to determine if it ends with something.
    endings: iterable
        A list of endings to check.
    addition: str
        A constant string to search for that is added to all endings.

    Returns
    -------
    True if name ends with something in endings plus the addition.

    """
    for end in endings:
        if name.endswith(end + addition):
            return True
    return False


def genes_on_chrom(chrom, db):
    """Yield genes on `chrom`, sorted by start position"""
    for g in db.features_of_type(EXCLUDED_REGIONS,
                                 order_by='start',
                                 limit=(chrom, 0, 1e12)):
        g.strand = '.'
        yield g


def intergenic(chroms, db):
    """Yield intergenic features"""
    for chrom in chroms:
        genes = genes_on_chrom(chrom, db)
        for intergenic in db.interfeatures(genes):
            yield intergenic


def process_gff(db):
    """
    Yield the exon, gene features from a gff file.

    Parameters
    ----------
    db: database
        gffutils database

    Returns
    -------
    yield gff features

    """
    chroms = [i['seqid'] for i in db.execute('SELECT DISTINCT seqid FROM '
                                             'features;')]
    for feature in intergenic(chroms, db):
        yield feature.astuple()


def get_db(gff_location):
    """
    Build or connect to a gffutils database.

    Parameters
    ----------
    gff_location: str
        The location of a gff file.

    Returns
    -------
    db: database
        An object representing a gffutils database.

    """
    if gff_location.endswith("gz"):
        subprocess.call(['gunzip', gff_location])
        gff_location = gff_location[0:len(gff_location)-3]
    if not os.path.isfile(gff_location + ".db"):
        db = gffutils.create_db(gff_location, gff_location + ".db",
                                id_spec={'gene': 'db_xref'})
    else:
        db = gffutils.FeatureDB(gff_location + ".db")
    return db


def get_intercds_cds(genome_location, fai_location,
                     cds_location, output, verbose=0):
    """
    Write a fasta file of intercds DNA.

    Parameters
    ----------
    genome_location: str
        The location of the whole genome fasta file.
    fai_location: str
        The location of the fai file for the whole genome file.
    cds_location: str
        The location of the cds fasta file.
    output: str or writable
        The place to write the output.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        The number of intercds fasta records.

    """
    def write_sequence(seq_id, beg, end):
        """Write the sequence to a file.  Returns 0 if successful."""
        try:
            fai_line = fasta_dict[seq_id]
        except KeyError:
            print("The sequence accession {} was not found for "
                  "the genome {}.  Further sequence writing will be "
                  "skipped.".format(seq_id, genome_location),
                  file=sys.stderr)
            return 1
        offset = fai_line[1]
        linebases = fai_line[2]
        if not end:
            end = fai_line[0]
        n_beg = offset + beg + (math.ceil(float(beg)/float(linebases))-1)
        n_end = offset + end + (math.ceil(float(end)/float(linebases))-1)
        try:
            inter_seq = str(mm[int(n_beg):int(n_end)],
                            'utf-8').replace("\n", "")
        except TypeError:
            inter_seq = str(mm[int(n_beg):int(n_end)]).replace("\n", "")
        inter_id = "{}:{}:{}:{}".format(seq_id, beg+1, end+1, end-beg)
        writer.write((inter_id, inter_seq))
        return 0

    if isinstance(output, str):
        output = open(output, "w")
    writer = SeqWriter(output, file_type="fasta")
    cds_gzip = True if cds_location.endswith("gz") else False
    cds_reader = SeqReader(cds_location, gzip_switch=cds_gzip)
    location_pattern = re.compile("location=[0-9|a-zA-z|)|(|..|,|>|<|=]+")
    range_pattern = re.compile("[0-9]+\.\.[0-9]+")
    loc_split = re.compile("\.+")
    # contig_id = re.compile("\|.+_cds_")
    cds_locations = defaultdict(list)
    for cds_record in cds_reader:
        cds_header = cds_record[0]
        try:
            seq_id = cds_header.split("_cds_")[0].split("|")[1]
        except IndexError:
            print("A sequence id could not be found for {}".format(cds_header),
                  file=sys.stderr)
            continue
        try:
            raw_locations = re.findall(range_pattern,
                                       re.findall(location_pattern,
                                                  cds_header)[0])
            location = [tuple(map(int, re.split(loc_split, loc)))
                        for loc in raw_locations]
            cds_locations[seq_id].extend(location)
        except (IndexError, ValueError):
            print("The locations could not be found for {}".format(cds_header),
                  file=sys.stderr)
            continue
    fasta_file = open(genome_location, "r")
    mm = mmap.mmap(fasta_file.fileno(), 0, access=mmap.ACCESS_READ)
    fasta_dict = {}
    for line in open(fai_location):
        fai_line = line.split()
        fasta_dict[fai_line[0]] = list(map(int, fai_line[1:]))
    count = 0
    for seq_id in cds_locations:
        locations = cds_locations[seq_id]
        locations.sort()
        cds_locations[seq_id] = locations
        end = 1
        for cds_loc in locations:
            if cds_loc[0] > end:
                if verbose:
                    print(cds_loc, seq_id, end, cds_loc[0], file=sys.stderr)
                if write_sequence(seq_id, end-1, cds_loc[0]-1):
                    return 0
                count += 1
            end = cds_loc[1]
        write_sequence(seq_id, end, False)
        count += 1
    return count


def get_intergenic_fasta(fasta_location, gff_location, fai_location,
                         output, remove=True, verbose=0):
    """
    Write a fasta file of intercds DNA.

    Parameters
    ----------
    fasta_location: str
        The location of a fasta file.
    gff_location: str
        The location of a gff file that corresponds to the fasta file.
    fai_location: str
        The location of the fai file corresponding to the fasta file.
    output: writable, str
        An object that represents where to write the intergenic fasta file.
    remove: boolean
        If True, remove the initial fasta files and the fai file if it exists.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        A count of the fasta records written.

    """
    def intergenic_boundary(region, begin=True):
        """
        Rewrite the boundary between chromosomes to
        exclude excluded regions.

        Parameters
        ----------
        feature: list
            A list of three elements: chromosome name, beginning position,
            ending position
        begin: boolean
            True if the region is at the beginning of the chromosome and
            False if the region is at the ending of the chromosome.

        Returns
        -------
        region: list
            A three element list comprising a region

        """
        if begin:
            boundary_region = db.region(seqid=region[0],
                                        start=region[1],
                                        end=region[2])
        else:
            boundary_region = db.region(seqid=region[0],
                                        start=region[1])
        selector = 4 if begin else 5
        bound = float('inf') if begin else -float('inf')
        cp = operator.lt if begin else operator.gt
        for feature in boundary_region:
            feature = feature.astuple()
            if feature[3] in EXCLUDED_REGIONS and cp(int(feature[selector]),
                                                     bound):
                bound = int(feature[selector])
        if begin:
            return [region[0], region[1], bound]
        else:
            return [region[0], bound, region[2]]

    def get_sequence(seq_id, beg, end=""):
        """Return a DNA sequence.  Use either fasta file or fai file."""
        if fai_location:
            fai_line = fasta_dict[seq_id]
            offset = fai_line[1]
            linebases = fai_line[2]
            if not end:
                end = fai_line[0]
            n_beg = offset + beg + (math.ceil(float(beg)/float(linebases))-1)
            n_end = offset + end + (math.ceil(float(end)/float(linebases))-1)
            return str(mm[int(n_beg):int(n_end)]).replace("\n", "")
        elif not fai_location and not end:
            return fasta_dict[seq_id][beg:]
        else:
            return fasta_dict[seq_id][beg:end]

    def write_feature(feature, begin_end="middle"):
        def fasta_id(region):
            return "{}_{}_{}".format(region[0], region[1], region[2])

        """Write a region of the genome to a fasta file."""
        if begin_end == "end":
            feature = intergenic_boundary(feature, begin=False)
            writer.write((fasta_id(feature),
                         get_sequence(feature[0], feature[1])))
            return
        if begin_end == "begin":
            feature = intergenic_boundary(feature, begin=True)
        writer.write((fasta_id(feature),
                      get_sequence(feature[0], feature[1], feature[2])))

    def chromosome_break(this_feature):
        """
        Handle the boundary between chromosomes.
        This requires special handling
        since gffutils does not handle this case.
        """
        this_count = 0
        if not this_feature[0] == last_feature[0]:
            if last_feature[0]:  # The end
                write_feature([last_feature[0], last_feature[2], "end"],
                              "end")
                this_count += 1
            if this_feature[0]:  # The beginning
                write_feature([this_feature[0], 0, this_feature[1]],
                              "begin")
                this_count += 1
        return this_feature, this_count

    count = 0
    if isinstance(output, str):
        output = open(output, "w")
    if fasta_location.endswith("gz"):
        gzip_switch = True
    else:
        gzip_switch = False
    writer = SeqWriter(output, file_type="fasta")
    fasta_dict = {}
    if not fai_location:
        reader = SeqReader(fasta_location, file_type="fasta",
                           gzip_switch=gzip_switch)
        if verbose:
            print("Loading {} into memory.".format(fasta_location),
                  file=sys.stderr)
        for fasta_record in reader:
            fasta_dict[fasta_record[0].split()[0]] = fasta_record[1]
    else:
        if verbose:
            print("Loading {} into memory.".format(fai_location),
                  file=sys.stderr)
        for line in open(fai_location):
            fai_line = line.split()
            fasta_dict[fai_line[0]] = list(map(int, fai_line[1:]))
        fasta_file = open(fasta_location, "r")
        mm = mmap.mmap(fasta_file.fileno(), 0, access=mmap.ACCESS_READ)
    if verbose:
        print("Creating or getting the gff database {}".format(gff_location),
              file=sys.stderr)
    last_feature = [False]
    db = get_db(gff_location)
    if verbose:
        print("Creating the output fasta file for {}.".format(fasta_location),
              file=sys.stderr)
    for feature in process_gff(db):
        this_feature = [feature[1], int(feature[4]) - 1, int(feature[5]) - 1]
        if this_feature[2] <= this_feature[1]:
            continue
        if (verbose >= 2 and last_feature[0] and
                last_feature[0] != this_feature[0]):
            print("Finishing with contig {} from {}".format(last_feature[0],
                                                            fasta_location),
                  file=sys.stderr)
        last_feature, this_count = chromosome_break(this_feature)
        write_feature(this_feature)
        count += this_count + 1
    last_feature, this_count = chromosome_break([False])
    if verbose:
        print("Finished with {}".format(fasta_location), file=sys.stderr)
    if remove:
        os.remove(fasta_location)
        if fai_location:
            os.remove(fai_location)
    return count + this_count


def process_directory(path, files, remove=True, verbose=0):
    """
    Create the intergenic fasta file a directory.

    Parameters
    ----------
    path: str
        The location of the directory.
    files: list
        A list of the files at the directory in path
    remove: boolean
        If True, remove the initial fasta files and the fai file if it exists.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        A count of the fasta files created.  Should be a 0 or a 1.

    """
    fasta_file = None
    gff_file = None
    fai_file = None
    count = 0
    for f in files:
        if (name_ends(f, FASTA_ENDINGS) or
                name_ends(f, FASTA_ENDINGS, ".gz")):
            fasta_file = f
        if (name_ends(f, GFF_ENDINGS) or
                name_ends(f, GFF_ENDINGS, ".gz")):
            gff_file = f
        if f.endswith(".fai"):
            fai_file = f
    if fasta_file and gff_file:
        output = os.path.join(path,
                              "{}.intergenic.fna".format(
                                  fasta_file.replace(".fna", "").replace(
                                      ".gz", "")))
        fai_location = os.path.join(path, fai_file) if fai_file else None
        fasta_records = get_intergenic_fasta(os.path.join(path,
                                                          fasta_file),
                                             os.path.join(path,
                                                          gff_file),
                                             fai_location,
                                             output,
                                             remove=remove,
                                             verbose=verbose)
        if fasta_records:
            count += 1
    return count


def scale_up(root_directory, processes=1, remove=True, verbose=0):
    """
    Create intergenic fasta files for all fasta files and GFF files.

    Parameters
    ----------
    root_directory: str
        The location of the directory where all of the GFF and fasta files
        are located.
    processes: int
        The number of processes to use.
    remove: boolean
        If True, remove the initial fasta files and the fai file if it exists.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        The number of fasta records created.

    """
    count = 0
    if processes == 1:
        for path, _, files in os.walk(root_directory):
            count += process_directory(path, files, remove, verbose)
    else:
        pd_list = []
        pool = Pool(processes=processes)
        for path, _, files in os.walk(root_directory):
            pd = pool.apply_async(process_directory, args=(path, files,
                                                           remove,
                                                           verbose))
            pd_list.append(pd)
        for pd in pd_list:
            count += pd.get()
        pool.close()
        pool.join()
    return count


def return_fasta(files, file_type='fasta'):
    """
    Get a fasta file from a list of files.

    Parameters
    ----------
    files: list
        A list of strings representing files
    file_type: str
        Either 'fasta' or 'fai'.  Gets those files.

    Returns
    -------
    f: str
        A fasta file if it exists.

    """
    ending = FASTA_ENDINGS
    if file_type == 'fai':
        ending = FAI_ENDINGS
    for f in files:
        if name_ends(f, ending):
            return f
    return None


def process_cd_directory(cds_directory, files,
                         genome_directory, intercds,
                         verbose=0):
    """
    Produce intercds fasta files for a CDs directory.

    cds_directory: str
        The cds directory ending with a genome accession.
    files: iterable
        A list of files in the cds_directory.
    genome_directory: str
        The base directory for whole genomes.
    intercds: str
        The base directory to write intercds fasta files.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        A count of the fasta files created.

    """
    f_cds = return_fasta(files)
    if cds_directory.endswith('/'):
        cds_directory = cds_directory[0:len(cds_directory)-1]
    accession = os.path.basename(cds_directory)
    genome_dir_accession = os.path.join(genome_directory, accession)
    f_genome = None
    f_fai = None
    if os.path.exists(genome_dir_accession):
        files_g = [f for f in os.listdir(genome_dir_accession)
                   if os.path.isfile(os.path.join(genome_dir_accession, f))]
        f_genome = return_fasta(files_g)
        f_fai = return_fasta(files_g, file_type='fai')
    if verbose >= 1:
        print("Getting intercds file for {} and {} with fai file {}.".format(
            f_cds, f_genome, f_fai), file=sys.stderr)
        sys.stderr.flush()
    if f_cds and f_genome:
        intercds_accession = os.path.join(intercds, accession)
        os.makedirs(intercds_accession)
        inter_filename = f_genome.replace(".fna", "") + "_intercds.fna"
        with open(os.path.join(intercds_accession, inter_filename), "w") as fd:
            cnt_cds = get_intercds_cds(os.path.join(genome_dir_accession,
                                                    f_genome),
                                       os.path.join(genome_dir_accession,
                                                    f_fai),
                                       os.path.join(cds_directory, f_cds),
                                       fd,
                                       verbose=1 if verbose >= 2 else 0)
            if cnt_cds > 0:
                return 1
    return 0


def scale_up_cds(cds_directory, genome_directory,
                 intercds, processes, verbose=0):
    """
    Create intercds files from cds directories.

    Parameters
    ----------
    cds_directory: str
        The directory where cds are located.
    genome_directory: str
        The directory where whole genomes are located.
    intercds: str
        The directory where the intercds fasta files will be written.
    processes: int
        The number of processes to use.
    verbose: int
        The level of verbosity.

    Returns
    -------
    count: int
        The number of fasta files created.

    """
    count = 0
    total = 0
    if processes == 1:
        for path, _, files in os.walk(cds_directory):
            count += process_cd_directory(path, files, genome_directory,
                                          intercds, verbose)
            total += 1
    else:
        pd_list = []
        pool = Pool(processes=processes)
        for path, _, files in os.walk(cds_directory):
            pd = pool.apply_async(process_cd_directory,
                                  args=(path, files,
                                        genome_directory, intercds, verbose))
            total += 1
            pd_list.append(pd)
        for pd in pd_list:
            count += pd.get()
        pool.close()
        pool.join()
    return count, total


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    print("Intergenic_fasta was started at {}.".format(tick), file=sys.stderr)
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_one_gff = subparsers.add_parser("one_gff",
                                      help=("Extract complementary regions "
                                            "from a single fasta file.  "
                                            "Use a gff file and a "
                                            "whole genome file."),
                                      formatter_class=argparse.
                                      ArgumentDefaultsHelpFormatter)
    p_one_gff.add_argument("fasta_file", type=str,
                           help="The location of the fasta file.")
    p_one_gff.add_argument("--keep", "-k",
                           help=("Keep the initial fasta files and the "
                                 "fai file "
                                 "if it exists.  Otherwise, these will be "
                                 "deleted."),
                           action='store_true', default=False)
    p_one_gff.add_argument("--verbose", "-v", type=int,
                           help="The level of verbosity.",
                           default=1)
    p_one_gff.add_argument("gff_file", type=str,
                           help="The location of the gff file.")
    p_one = subparsers.add_parser("one_cds",
                                  help="Extract complementary regions "
                                  "from a single fasta file using a CDs "
                                  "file and a whole genome file.",
                                  formatter_class=argparse.
                                  ArgumentDefaultsHelpFormatter)
    p_one.add_argument("genome_file", type=str,
                       help="The location of the whole genome file.")
    p_one.add_argument("fai_file", type=str,
                       help="The fai file for the whole genome file.")
    p_one.add_argument("cds_file", type=str,
                       help="The location of the CDs file.")
    p_one.add_argument("--verbose", "-v", type=int,
                       help="The level of verbosity.",
                       default=1)
    p_all_gff = subparsers.add_parser("all_gff",
                                      help=("Create complementary fasta files "
                                            "for all genomes with gff files."),
                                      formatter_class=argparse.
                                      ArgumentDefaultsHelpFormatter)
    p_all_gff.add_argument("gff_directory", type=str,
                           help="The location of the directory of gff files.")
    p_all_gff.add_argument("--processes", "-p", type=int,
                           help="The number of processes to use.",
                           default=1)
    p_all_gff.add_argument("--keep", "-k",
                           help=("Keep the initial fasta files and the fai "
                                 "file "
                                 "if it exists.  Otherwise, these will be "
                                 "deleted."),
                           action='store_true', default=False)
    p_all_gff.add_argument("--verbose", "-v", type=int,
                           help="The level of verbosity.",
                           default=0)
    p_all_cds = subparsers.add_parser("all_cds",
                                      help=("Create intercds fasta files "
                                            "for all genomes with cds files."),
                                      formatter_class=argparse.
                                      ArgumentDefaultsHelpFormatter)
    p_all_cds.add_argument("cds_directory", type=str,
                           help="The location of the directory where the cds "
                           "files are located.")
    p_all_cds.add_argument("genome_directory", type=str,
                           help="The location of the whole genomes directory.")
    p_all_cds.add_argument("--intercds", "-i", type=str,
                           help="The directory to create intercds "
                           "fasta files.",
                           default="./")
    p_all_cds.add_argument("--processes", "-p", type=int,
                           help="The number of processes to use.",
                           default=20)
    p_all_cds.add_argument("--verbose", "-v", type=int,
                           help="The level of verbosity.",
                           default=0)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    mode = args.mode
    if mode == "one_gff":
        count = get_intergenic_fasta(args.fasta_file, args.gff_file,
                                     sys.stdout, remove=not args.keep,
                                     verbose=args.verbose)
        print("There were {} fasta records written.".format(count),
              file=sys.stderr)
    elif mode == "one_cds":
        count = get_intercds_cds(args.genome_file, args.fai_file,
                                 args.cds_file,
                                 sys.stdout, verbose=args.verbose)
        print("There were {} fasta records written.".format(count),
              file=sys.stderr)
    elif mode == "all_gff":
        count = scale_up(args.gff_directory, args.processes,
                         remove=not args.keep,
                         verbose=args.verbose)
        print("There were {} fasta files created.".format(count),
              file=sys.stderr)
    elif mode == "all_cds":
        count = scale_up_cds(args.cds_directory, args.genome_directory,
                             args.intercds, args.processes, args.verbose)
        print("There were {} fasta files created out of {} "
              "directories.".format(count[0], count[1]),
              file=sys.stderr)
    else:
        parser.error("The mode was not recognized.  Please to check.")
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
