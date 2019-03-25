#!/usr/bin/env python2
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
from multiprocessing import Pool

import gffutils

from SeqIterator import SeqReader, SeqWriter

# GFF regions to exclude when creating fasta file.
EXCLUDED_REGIONS = ('gene', 'exon')

FASTA_ENDINGS = ['fasta', 'fa', 'fna']
GFF_ENDINGS = ['gff']


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
        subprocess.run(['gunzip', gff_location])
        gff_location = gff_location[0:len(gff_location)-3]
    if not os.path.isfile(gff_location + ".db"):
        db = gffutils.create_db(gff_location, gff_location + ".db",
                                id_spec={'gene': 'db_xref'})
    else:
        db = gffutils.FeatureDB(gff_location + ".db")
    return db


def get_intergenic_fasta(fasta_location, gff_location, output, verbose=0):
    """
    Write a fasta file of intergenic DNA.

    Parameters
    ----------
    fasta_location: str
        The location of a fasta file.
    gff_location: str
        The location of a gff file that corresponds to the fasta file.
    output: writable, str
        An object that represents where to write the intergenic fasta file.

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

    def write_feature(feature, begin_end="middle"):
        def fasta_id(region):
            return "{}_{}_{}".format(region[0], region[1], region[2])

        """Write a region of the genome to a fasta file."""
        if begin_end == "end":
            feature = intergenic_boundary(feature, begin=False)
            writer.write((fasta_id(feature),
                         fasta_dict[feature[0]][feature[1]:]))
            return
        if begin_end == "begin":
            feature = intergenic_boundary(feature, begin=True)
        writer.write((fasta_id(feature),
                      fasta_dict[feature[0]][feature[1]:feature[2]]))

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
    reader = SeqReader(fasta_location, file_type="fasta",
                       gzip_switch=gzip_switch)
    writer = SeqWriter(output, file_type="fasta")
    fasta_dict = {}
    for fasta_record in reader:
        fasta_dict[fasta_record[0].split()[0]] = fasta_record[1]
    last_feature = [False]
    db = get_db(gff_location)
    for feature in process_gff(db):
        this_feature = [feature[1], int(feature[4]) - 1, int(feature[5]) - 1]
        if this_feature[2] <= this_feature[1]:
            continue
        if verbose and last_feature[0] and last_feature[0] != this_feature[0]:
            print("Finishing with contig {}".format(last_feature[0]),
                  file=sys.stderr)
        last_feature, this_count = chromosome_break(this_feature)
        write_feature(this_feature)
        count += this_count + 1
    last_feature, this_count = chromosome_break([False])
    return count + this_count


def process_directory(path, files):
    """
    Create the intergenic fasta file a directory.

    Parameters
    ----------
    path: str
        The location of the directory.
    files: list
        A list of the files at the directory in path

    Returns
    -------
    count: int
        A count of the fasta files created.  Should be a 0 or a 1.

    """
    fasta_file = None
    gff_file = None
    for f in files:
        if (name_ends(f, FASTA_ENDINGS) or
                name_ends(f, FASTA_ENDINGS, ".gz")):
            fasta_file = f
        if (name_ends(f, GFF_ENDINGS) or
                name_ends(f, GFF_ENDINGS, ".gz")):
            gff_file = f
    if fasta_file and gff_file:
        output = os.path.join(path,
                              "{}.intergenic.fna".format(
                                  fasta_file.replace(".fna", "").replace(
                                      ".gz", "")))
        fasta_records = get_intergenic_fasta(os.path.join(path,
                                                          fasta_file),
                                             os.path.join(path,
                                                          gff_file),
                                             output)
        if fasta_records:
            return 1
    if fasta_file:
        os.remove(os.path.join(path, fasta_file))
    return 0


def scale_up(root_directory, processes=1):
    """
    Create intergenic fasta files for all fasta files and GFF files.

    Parameters
    ----------
    root_directory: str
        The location of the directory where all of the GFF and fasta files
        are located.

    Returns
    -------
    count: int
        The number of fasta records created.

    """
    count = 0
    if processes == 1:
        for path, _, files in os.walk(root_directory):
            count += process_directory(path, files)
    else:
        pool = Pool(processes=processes)
        pd_list = []
        for path, _, files in os.walk(root_directory):
            pd = pool.apply_async(process_directory, args=(path, files))
        for pd in pd_list:
            count += pd.get()
    return count


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    print("Intergenic_fasta was started at {}.".format(tick), file=sys.stderr)
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=__doc__)
    subparsers = parser.add_subparsers(help="sub-commands", dest="mode")
    p_one = subparsers.add_parser("one",
                                  help=("Extract complementary regions "
                                        "from a single fasta file."),
                                  formatter_class=argparse.
                                  ArgumentDefaultsHelpFormatter)
    p_one.add_argument("fasta_file", type=str,
                       help="The location of the fasta file.")
    p_one.add_argument("gff_file", type=str,
                       help="The location of the gff file.")
    p_all = subparsers.add_parser("all",
                                  help=("Create complementary fasta files "
                                        "for all genomes with gff files."),
                                  formatter_class=argparse.
                                  ArgumentDefaultsHelpFormatter)
    p_all.add_argument("gff_directory", type=str,
                       help="The location of the directory of gff files.")
    p_all.add_argument("--processes", "-p", type=int,
                       help="The number of processes to use.",
                       default=1)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    sys.stderr.flush()
    mode = args.mode
    if mode == "one":
        count = get_intergenic_fasta(args.fasta_file, args.gff_file,
                                     sys.stdout, verbose=1)
        print("There were {} fasta records written.".format(count),
              file=sys.stderr)
    elif mode == "all":
        count = scale_up(args.gff_directory, args.processes)
        print("There were {} fasta files created.".format(count),
              file=sys.stderr)
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":

    main()
