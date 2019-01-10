"""
Run faidx from samtools on all genomes below the given directory.
Decompress gzip comressed fasta files.  Create an index
of genomic and taxonomic information.
This program requires SAMTOOLS and UCSC tools.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import os
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from config import SAMTOOLS

FASTA_ENDINGS = ['fasta', 'fa', 'fna', 'faa']

VERBOSE_COUNTER = 5000

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


def make_fai_individual(root_directory, path, files,
                        leave_compressed=False, verbose=0):
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
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.

    Returns
    -------
    count
        A dictionary that gives which types of files were processed.

    """
    count = defaultdict(int)
    fasta_exists = False
    for f in files:
        if (name_ends(f, FASTA_ENDINGS, addition=".gz") or
                name_ends(f, FASTA_ENDINGS, addition="")):
            fasta_exists = True
            break
    if not fasta_exists:
        return None
    for name in files:  # Should we always delete FAI files.
        if (name.endswith('fai')):
            os.remove(os.path.join(path, name))
    for name in files:
        if (not leave_compressed and
                name_ends(name, FASTA_ENDINGS, addition='.gz')):
            if verbose >= 2:
                sys.stderr.write(name + "\n")
            subprocess.run(["gunzip", os.path.join(path, name)])
            name = name[0:len(name) - 3]
            count["gzip"] += 1
        if (name_ends(name, FASTA_ENDINGS)):
            fa_file = os.path.join(path, name)
            if verbose >= 2:
                sys.stderr.write(name + "\n")
            subprocess.run([SAMTOOLS + "samtools", "faidx", fa_file])
            count["fai"] += 1
    return count


# def make_fai_multi(index, root_directory, twoBitToFa=False,
#                    leave_compressed=False, verbose=0,
#                    index_only=False, workers=16):
#     """
#     Make fai files for all fasta files under the root_directory.  Add contig
#     counts and fasta base length to the index.
# 
#     Parameters
#     ----------
#     index: dict
#         A python dictionary representing indexed genome information.
#     root_directory: str
#         A string representing the location of the root directory to search.
#     twoBitToFasta: bool
#         If true, convert 2bit files to fasta.  Does not create a fai index.
#     leave_compressed: bool
#         Does not gunzip anything.
#     verbose: int
#         If larger than one, print additional output.
#     index_only: bool
#         Updates the index in index only.  Assumes that the fai files have
#         already been created.
#     workers: int
#         The number of workers to use for multiprocessing.
# 
#     Returns
#     -------
#     int, index
#         A count of the fasta files processed and the index.
# 
#     Examples
#     --------
#         make_fai({}, '/home/jsporter/Genomes/')
# 
#     """
#     count = {"gzip": 0, "2bit": 0, "fai": 0}
#     counter = 0
#     pool = Pool(processes=workers)
#     pd_list = []
#     for path, _, files in os.walk(root_directory):
#         # print(path, files)
#         pd = pool.apply_async(make_fai_individual, (root_directory,
#                                                     path,
#                                                     files,
#                                                     twoBitToFa,
#                                                     leave_compressed,
#                                                     verbose,
#                                                     index_only))
#         pd_list.append(pd)
#     pool.close()
#     pool.join()
#     for pd in pd_list:
#         count_local, contigs, genome_length = pd.get()
#         accession = os.path.basename(path)
#         if contigs is not None and genome_length is not None:
#             index['genomes'][accession]['contig_count'] = contigs
#             index['genomes'][accession]['contig_sum'] = genome_length
#         else:
#             print(path, files)
#         for key in count_local:
#             count[key] += count_local[key]
#         counter += 1
#         if counter >= 5000 and verbose >= 1:
#             sys.stderr.write("Processed: {}\n".format(counter))
#             sys.stderr.flush()
#     sys.stderr.flush()
#     return count, index


def make_fai(root_directory, leave_compressed=False, verbose=0):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        Controls verbosity of output.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    Examples
    --------
        make_fai({}, '/home/jsporter/Genomes/')

    """
    count = {"gzip": 0, "fai": 0}
    counter = 0
    for path, _, files in os.walk(root_directory):
        counter += 1
        if counter % VERBOSE_COUNTER == 0 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
        count_local = make_fai_individual(
            root_directory, path, files,
            leave_compressed, verbose)
        if count_local:
            for key in count_local:
                count[key] += count_local[key]
    sys.stderr.flush()
    return count