#!/usr/bin/env python
"""
Run faidx from samtools on all genomes below the given directory.
Decompress gzip comressed fasta files.  Create an index
of genomic and taxonomic information.
This program requires SAMTOOLS and UCSC tools.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import os
import subprocess
import pickle
import datetime
import sys
import json
from collections import defaultdict
from get_paths import get_paths
from multiprocessing import Pool

path_dict = get_paths()

SAMTOOLS = path_dict['SAMTOOLS']
UCSC = path_dict['UCSC']

FASTA_ENDINGS = ['fasta', 'fa', 'fna', 'faa']

def process_fai(fai_file_fd):
    """
    Use an fai file to update the index.

    Parameters
    ----------
    fai_file_fd: file object
        An open file object for the fai file.

    Returns
    -------
    contigs, genome_length
        The number of contigs followed by the number of nucleotide bases.

    """
    contigs = 0
    genome_length = 0
    for line in fai_file_fd:
        contigs += 1
        genome_length += int(line.split('\t')[1])
    return contigs, genome_length


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


def make_fai_individual(root_directory, path, files, twoBitToFa=False,
                        leave_compressed=False, verbose=0, index_only=False):
    """
    Make the fai files and decompress the fasta files.

    Parameters
    ----------
    root_directory: str
        A string giving the location of the root directory where genomes
        are found.
    path: str
        The path to a directory containing files
    files: iterable
        An iterable of files under the path.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.

    Returns
    -------
    count, contig, genome_length
        A dictionary that gives which types of files were processed.
        The number of contigs in the fasta file that was processed.
        The number of bases in the fasta file.

    """
    count = defaultdict(int)
    contigs = None
    genome_length = None
    fasta_exists = False
    for f in files:
        if (name_ends(f, FASTA_ENDINGS, addition=".gz") or
                name_ends(f, FASTA_ENDINGS, addition="")):
            fasta_exists = True
            break
    if not fasta_exists:
        return None, None, None
    # if twoBitToFa:
    #     skip = False
    #     for name in files:
    #         if (name_ends(name, FASTA_ENDINGS)):
    #             skip = True
    #             break
    #     if skip:
    #         for name in files:
    #             if name.endswith('.2bit'):
    #                 os.remove(os.path.join(path, name))
    #         return count, contigs, genome_length
    if not index_only:
        for name in files:
            if (name.endswith('fai')):
                os.remove(os.path.join(path, name))
    for name in files:
        if index_only:
            if name.endswith('.fai'):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                fai_fd = open(os.path.join(path, name))
                contigs, genome_length = process_fai(fai_fd)
                count["fai"] += 1
            else:
                continue
        else:
            if (not leave_compressed and
                    name_ends(name, FASTA_ENDINGS, addition='.gz')):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                subprocess.run(["gunzip", os.path.join(path, name)])
                name = name[0:len(name) - 3]
                count["gzip"] += 1
            if name.endswith('.2bit'):
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                new_name = name[0:len(name)-5]
                fa_file = os.path.join(path, new_name)
                twobit_process = subprocess.run([UCSC + "twoBitToFa",
                                                 os.path.join(path, name),
                                                 fa_file])
                if twobit_process.returncode == 0:
                    os.remove(os.path.join(path, name))
                name = new_name
                count["2bit"] += 1
            if (name_ends(name, FASTA_ENDINGS) and not twoBitToFa):
                fa_file = os.path.join(path, name)
                if verbose >= 2:
                    sys.stderr.write(name + "\n")
                subprocess.run([SAMTOOLS + "samtools", "faidx", fa_file])
                count["fai"] += 1
                contigs, genome_length = process_fai(
                    open(os.path.join(path, name + '.fai')))
    return count, contigs, genome_length


def make_fai_multi(index, root_directory, twoBitToFa=False,
                   leave_compressed=False, verbose=0,
                   index_only=False, workers=16):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.
    workers: int
        The number of workers to use for multiprocessing.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    Examples
    --------
        make_fai({}, '/home/jsporter/Genomes/')

    """
    count = {"gzip": 0, "2bit": 0, "fai": 0}
    counter = 0
    pool = Pool(processes=workers)
    pd_list = []
    for path, _, files in os.walk(root_directory):
        # print(path, files)
        pd = pool.apply_async(make_fai_individual, (root_directory,
                                                    path,
                                                    files,
                                                    twoBitToFa,
                                                    leave_compressed,
                                                    verbose,
                                                    index_only))
        pd_list.append(pd)
    pool.close()
    pool.join()
    for pd in pd_list:
        count_local, contigs, genome_length = pd.get()
        accession = os.path.basename(path)
        if contigs is not None and genome_length is not None:
            index['genomes'][accession]['contig_count'] = contigs
            index['genomes'][accession]['base_length'] = genome_length
        else:
            print(path, files)
        for key in count_local:
            count[key] += count_local[key]
        counter += 1
        if counter >= 5000 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
    sys.stderr.flush()
    return count, index


def remove_accessions(index, accessions_to_remove):
    """
    Remove genomes from the index.

    Parameters
    ----------
    index: dict
        The genomes index.
    accessions_to_remove: list
        A list of accessions to remove from the index.

    Returns
    -------
    dict
        The genomes index is returned.

    """
    remove_count = 0
    for taxid in index['taxids']:
        my_accessions_to_remove = []
        for accession in index['taxids'][taxid]:
            if accession in accessions_to_remove:
                my_accessions_to_remove.append(accession)
        for accession in my_accessions_to_remove:
            del index['taxids'][taxid][accession]
            remove_count += 1
    return index, remove_count


def make_fai(index, root_directory, twoBitToFa=False,
             leave_compressed=False, verbose=0,
             index_only=False):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    output_fd: writable (such as a file)
        A writable object such as an open file to periodically save the index.
    twoBitToFasta: bool
        If true, convert 2bit files to fasta.  Does not create a fai index.
    leave_compressed: bool
        Does not gunzip anything.
    index_only: bool
        Updates the index in index only.  Assumes that the fai files have
        already been created.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    Examples
    --------
        make_fai({}, '/home/jsporter/Genomes/')

    """
    count = {"gzip": 0, "2bit": 0, "fai": 0}
    counter = 0
    accessions_to_remove = []
    for path, _, files in os.walk(root_directory):
        counter += 1
        if counter % 5000 == 0 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
        count_local, contigs, genome_length = make_fai_individual(
            root_directory, path, files, twoBitToFa,
            leave_compressed, verbose, index_only)
        accession = os.path.basename(path)
        if count_local is None:
            accessions_to_remove.append(accession)
        else:
            if contigs is not None and genome_length is not None:
                index['genomes'][accession]['contig_count'] = contigs
                index['genomes'][accession]['base_length'] = genome_length
            for key in count_local:
                count[key] += count_local[key]
    index, remove_count = remove_accessions(index, accessions_to_remove)
    count['remove_count'] = remove_count
    sys.stderr.flush()
    return count, index


def main():
    """Parse arguments and dispatches fai file creation."""
    parser = argparse.ArgumentParser(description=('Make fasta index files '
                                                  'below the given directory. '
                                                  'Decompresses gzip fasta '
                                                  'files and converts them to '
                                                  '2bit files.'
                                                  'Add genome length and '
                                                  'contig count to the '
                                                  'index.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("genomes_directory",
                        type=str, help=("The directory storing genomes. "
                                        "All fasta "
                                        "files below here will be processed."))
    parser.add_argument("index", type=str,
                        help=("The location of the index to add to."),)
    parser.add_argument("--output", "-o", type=str,
                        help=("The file location to write the index to."),
                        default='./index.pck')
    parser.add_argument("--twoBitToFa", "-t", action='store_true',
                        help=("Convert any 2bit files to fasta.  "
                              "Does not build fai indexes."))
    parser.add_argument("--leave_compressed", "-l",
                        action='store_true',
                        help="Leave gzip files compressed.")
    parser.add_argument("--index_only", "-i",
                        action='store_true',
                        help=("Updates the index with the fai files "
                              "assuming that they are already created."))
    parser.add_argument("--workers", "-w", type=int,
                        help=("The number of workers to use.  "
                              "If there are more than one, then "
                              "multiprocessing will be used."),
                        default=1)
    parser.add_argument("--verbose", "-v", type=int,
                        help=("When 0 print basic info only. "
                              "When 1, print a counter of the files "
                              "processed.  When 2, the files that are being "
                              "worked on will be printed."),
                        default=0)
    args = parser.parse_args()
    now = datetime.datetime.now()
    if not args.output:
        args.output = args.index
    sys.stderr.write("Starting to go down {} to create fai files. "
                     "Adding the results to {} and writing the index "
                     "to {}.  twoBitToFa: {}, leave_compressed: {}, "
                     "index_only: {}, verbose: {}.\n"
                     .format(args.genomes_directory, args.index,
                             args.output, args.twoBitToFa,
                             args.leave_compressed, args.index_only,
                             args.verbose))
    index = pickle.load(open(args.index, 'rb'))
    if args.workers == 1:
        count, index2 = make_fai(index, args.genomes_directory,
                                 args.twoBitToFa, args.leave_compressed,
                                 args.verbose, args.index_only)
    elif args.workers > 1:
        raise NotImplementedError
        count, index2 = make_fai_multi(index, args.genomes_directory,
                                       args.twoBitToFa, args.leave_compressed,
                                       args.verbose, args.index_only,
                                       args.workers)
    else:
        parser.error("The number of workers must be larger than "
                     "or equal to 1.  Check the input.")
    pickle.dump(index2, open(args.output, 'wb'))
    # pprint(index2, open(args.output + ".txt", 'w'))
    json.dump(index2, open(args.output + '.json', 'w'), indent=4)
    later = datetime.datetime.now()
    sys.stderr.write("There were {} files processed.\n".format(count))
    sys.stderr.write("The process took time {}.\n".format(later - now))


if __name__ == '__main__':
    main()
