#!/usr/bin/env python
"""
Copy whole genome directories with matching CDs files into another directory.

:Authors:
    Jacob S. Porter <jsporter@vt.edu>

"""
import argparse
import datetime
import sys
import os
import shutil
import pickle
import re
import tqdm
from collections import defaultdict
from multiprocessing import Pool

from intercds_fasta import find_genome_files, name_ends, FASTA_ENDINGS
from SeqIterator import SeqReader

# The statuses of the DNA sequence.
CD = 0
INTERCD = 1
MIXED = 2
UNKNOWN = 3


def process_cd_directory(cds_directory, files,
                         genome_directory, new_dir,
                         verbose):
    """Copy a genome directory if there is a matching CD file."""
    genome_files = find_genome_files(files,
                                     cds_directory,
                                     genome_directory)
    f_cds, f_genome, f_fai, accession, genome_dir_accession = genome_files
    if f_cds and f_genome:
        os.makedirs(os.path.join(new_dir, accession))
        shutil.copy(os.path.join(genome_dir_accession, f_genome),
                    os.path.join(new_dir, accession))
        return 1
    return 0


def copy(genome_dir, cds_dir, new_dir, processes=1, verbose=0):
    """
    Copy directories with whole genomes and cds to a new location.

    Parameters
    ----------
    genome_dir: str
        The location of a directory with genomes.
    cds_dir: str
        The location of the cds directory.
    new_dir: str
        The directory to copy whole genome files to.
    processes: int
        The number of processes to use.
    verbose: int
        Controls the level of verbosity.

    Returns
    -------
    count: int
        The count of directories copied.
    total: int
        The number of directories examined.

    """
    count = 0
    total = 0
    if processes <= 1:
        for path, _, files in os.walk(cds_dir):
            count += process_cd_directory(path, files, genome_dir,
                                          new_dir, verbose)
            total += 1
    else:
        pd_list = []
        pool = Pool(processes=processes)
        for path, _, files in os.walk(cds_dir):
            pd = pool.apply_async(process_cd_directory,
                                  args=(path, files,
                                        genome_dir, new_dir, verbose))
            total += 1
            pd_list.append(pd)
        for pd in pd_list:
            count += pd.get()
        pool.close()
        pool.join()
    return count, total


def divide(fasta_file, cds_index, cds_dir, output, verbose=0):
    """Divide the fasta file based on the identity of the sequence."""
    def get_cds_list(cds_location, contig_count):
        cds_locations = defaultdict(list)
        if contig_count:
            files_g = [f for f in os.listdir(cds_location)
                       if os.path.isfile(os.path.join(cds_location, f)) and
                       (name_ends(f, FASTA_ENDINGS) or
                        name_ends(f, FASTA_ENDINGS, ".gz"))]
            cds_location = os.path.join(cds_location, files_g[0])
            cds_gzip = True if cds_location.endswith("gz") else False
            cds_reader = SeqReader(cds_location, gzip_switch=cds_gzip)
            for cds_record in cds_reader:
                cds_header = cds_record[0]
                try:
                    seq_id = cds_header.split("_cds_")[0].split("|")[1]
                except IndexError:
                    print("A sequence id could not be "
                          "found for {}".format(cds_header),
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
                    print("The locations could not be found "
                          "for {}".format(cds_header),
                          file=sys.stderr)
                    continue
            for seq_id in cds_locations:
                locations = cds_locations[seq_id]
                locations.sort()
                end = -1
                start = -1
                new_locations = []
                for location in locations:
                    if start == -1:
                        start = location[0]
                        end = location[1]
                        continue
                    if end < location[0]:
                        new_locations.append((start, end))
                        start = location[0]
                        end = location[1]
                    else:
                        end = location[1]
                if end != -1:
                    new_locations.append((start, end))
                cds_locations[seq_id] = new_locations
        return cds_locations

    def determine_identity(gen_accession, seq_accession, coords):
        if gen_accession not in intercds_dict:
            location = index["genomes"][gen_accession]['location']
            contig_count = index["genomes"][gen_accession]['contig_count']
            location = os.path.join(cds_dir, location)
            intercds_dict[gen_accession] = get_cds_list(location, contig_count)
        locations = intercds_dict[gen_accession]
        if locations:
            coords = tuple(map(int, coords.split("-")))
            for loc in locations:
                if coords[0] >= loc[0] and coords[1] < loc[1]:
                    return CD
                elif ((coords[0] >= loc[0] and coords[0] < loc[1]) or
                      (coords[1] >= loc[0] and coords[1] < loc[1]) or
                      (loc[0] >= coords[0] and loc[1] < coords[1])):
                    return MIXED
            return INTERCD
        else:
            return UNKNOWN
    if isinstance(output, str):
        output = open(output, "w")
    location_pattern = re.compile("location=[0-9|a-zA-z|)|(|..|,|>|<|=]+")
    range_pattern = re.compile("[0-9]+\.\.[0-9]+")
    loc_split = re.compile("\.+")
    index = pickle.load(open(cds_index, "rb"))
    reader = SeqReader(fasta_file)
    intercds_dict = {}
    type_counter = defaultdict(int)
    count = 0
    for record in tqdm(reader):
        coords = record[0].split(":")
        record_identity = determine_identity(coords[0], coords[2], coords[3])
        type_counter[record_identity] += 1
        print(record_identity, file=output)
        count += 1
    return type_counter


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_copy = subparsers.add_parser("copy",
                                   help="Copy whole genome directories "
                                   "that have cds.",
                                   formatter_class=argparse.
                                   ArgumentDefaultsHelpFormatter)
    p_copy.add_argument("genomes_dir", type=str,
                        help="The directory where whole genomes are "
                        "located.")
    p_copy.add_argument("cds_dir", type=str,
                        help="The directory where cds are located.")
    p_copy.add_argument("new_dir", type=str,
                        help="The directory to move genomes directories to.")
    p_copy.add_argument("--processes", "-p", type=int,
                        help="The number of processes to use.",
                        default=1)
    p_copy.add_argument("--verbose", "-v", type=int,
                        help="The level of verbosity.",
                        default=0)
    p_divide = subparsers.add_parser("divide",
                                     help="Partition a fasta file into "
                                     "cds, intercds, mixed, and unknown.",
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    p_divide.add_argument("fasta_file", type=str,
                          help="A fasta kmer file generated by Radogest.")
    p_divide.add_argument("cds_index", type=str,
                          help="The location of the Radogest CDs index.")
    p_divide.add_argument("cds_dir", type=str,
                          help="The location of the CDs directory.")
    p_divide.add_argument("--output", "-o", type=str,
                          help="A file to write the output.",
                          default=sys.stdout)
    p_divide.add_argument("--verbose", "-v", type=int,
                          help="The level of verbosity.",
                          default=0)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    mode = args.mode
    if mode == "copy":
        count, total = copy(args.genomes_dir, args.cds_dir,
                            args.new_dir, processes=args.processes,
                            verbose=args.verbose)
        print("There were {} directories copied out of {}.".format(count,
                                                                   total),
              file=sys.stderr)
    elif mode == "divide":
        print("Identity mappings: CD:{}, INTERCD:{}, "
              "MIXED:{}, UNKNOWN:{}".format(CD, INTERCD, MIXED, UNKNOWN),
              file=sys.stderr)
        counts = divide(args.fasta_file,
                        args.cds_index,
                        args.cds_dir,
                        args.output,
                        args.verbose)
        print("There were {} records.".format(counts),
              file=sys.stderr)
    else:
        parser.error("The mode was not recognized.  Please to check.")
    tock = datetime.datetime.now()
    print("The process took time {}.".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
