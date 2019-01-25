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


def make_fai_individual(path, files, leave_compressed=False, verbose=0):
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
    for name in files:  # Should we always delete FAI files.  Yes.
        if (name.endswith('fai')):
            os.remove(os.path.join(path, name))
    returncode = 0
    for name in files:
        if (not leave_compressed and
                name_ends(name, FASTA_ENDINGS, addition='.gz')):
            if verbose >= 2:
                sys.stderr.write(name + "\n")
            complete_p = subprocess.run(["gunzip", os.path.join(path, name)])
            returncode = returncode ^ complete_p.returncode
            name = name[0:len(name) - 3]
            if not complete_p.returncode:
                count["gzip"] += 1
        if (name_ends(name, FASTA_ENDINGS)):
            fa_file = os.path.join(path, name)
            if verbose >= 2:
                sys.stderr.write(name + "\n")
            complete_p = subprocess.run([
                SAMTOOLS + "samtools", "faidx", fa_file])
            returncode = returncode ^ complete_p.returncode
            if not complete_p.returncode:
                count["fai"] += 1
    return count, returncode, path


def fai_fail_message(returncode, path):
    """
    Prints out an error message for when faidx or when gzip fails.

    Parameters
    ----------
    returncode: int
        A 0 if the process existed normally.  Something else otherwise.
    path: str
        A string representing the location of the genome accession's files.

    Returns
    -------
    None

    """
    if returncode:
        print("Making the fai file for {} failed with return code {}.".format(
            returncode, path), file=sys.stderr)


def make_fai(root_directory, processes=1, leave_compressed=False, verbose=0):
    """
    Make fai files for all fasta files under the root_directory.  Add contig
    counts and fasta base length to the index.

    Parameters
    ----------
    root_directory: str
        A string representing the location of the root directory to search.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        Controls verbosity of output.

    Returns
    -------
    int
        A count of the fasta files processed.

    """
    def _count_add(count_local, count):
        if count_local:
            for key in count_local:
                count[key] += count_local[key]
    count = {"gzip": 0, "fai": 0}
    counter = 0
    if processes > 1:
        pool = Pool(processes=processes)
        pd_list = []
    for path, _, files in os.walk(root_directory):
        if processes > 1:
            pd = pool.apply_async(make_fai_individual,
                                  args=(path,
                                        files,
                                        leave_compressed,
                                        verbose))
            pd_list.append(pd)
        else:
            counter += 1
            if counter % VERBOSE_COUNTER == 0 and verbose >= 1:
                sys.stderr.write("Processed: {}\n".format(counter))
                sys.stderr.flush()
            count_local, returncode, _ = make_fai_individual(path,
                                                             files,
                                                             leave_compressed,
                                                             verbose)
            fai_fail_message(returncode, path)
            _count_add(count_local, count)
    if processes > 1:
        for pd in pd_list:
            count_local, returncode, path = pd.get()
            fai_fail_message(returncode, path)
            _count_add(count_local, count)
            counter += 1
            if counter % VERBOSE_COUNTER == 0 and verbose >= 1:
                sys.stderr.write("Processed: {}\n".format(counter))
                sys.stderr.flush()
    sys.stderr.flush()
    return count
