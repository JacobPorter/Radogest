"""
Run faidx from samtools on all genomes below the given directory.
Decompress gzip comressed fasta files.
This module requires SAMTOOLS.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import os
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool

try:
    from config import SAMTOOLS
except ImportError:
    SAMTOOLS = ""
from library.util import which

FASTA_ENDINGS = ['fasta', 'fa', 'fna', 'faa', 'jsonhfasta']

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
            return len(end + addition)
    return 0


def make_fai_individual(path,
                        files,
                        samtools_path,
                        leave_compressed=False,
                        verbose=0):
    """
    Make the fai files and decompress the fasta files.

    Parameters
    ----------
    path: str
        The path to a directory containing files.
    files: iterable
        An iterable of files under the path.
    samtools_path: str
        The location of the samtools executable.
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
        if (name_ends(f, FASTA_ENDINGS, addition=".gz")
                or name_ends(f, FASTA_ENDINGS, addition="")):
            fasta_exists = True
            break
    if not fasta_exists:
        return count, 1, path
    for name in files:  # Should we always delete FAI files.  Yes.
        remove_fai(path, name)
    returncode = 0
    for name in files:
        count_temp, returncode, _ = make_fai_core(path, name, samtools_path,
                                                  returncode, leave_compressed,
                                                  verbose)
        for key in count_temp:
            count[key] += count_temp[key]
    return count, returncode, path


def remove_fai(path, name):
    """
    Remove the file if it is an fai file.

    Parameters
    ----------
    path: str
        The path to a directory containing files.
    name: str
        The name of a file in path

    Returns
    -------
    None

    """
    if (name.endswith('fai')):
        os.remove(os.path.join(path, name))


def make_fai_core(path,
                  name,
                  samtools_path,
                  returncode=0,
                  leave_compressed=False,
                  verbose=0):
    """
    Remove the file if it is an fai file.

    Parameters
    ----------
    path: str
        The path to a directory containing files.
    name: str
        The name of a file in path
    samtools_path: str
        The location of the samtools executable.
    returncode: int
        A 0 if exiting normally, non-zero otherwise.
    leave_compressed: bool
        Does not gunzip anything.
    verbose: int
        If larger than one, print additional output.

    Returns
    -------
    count, returncode, path
        A dictionary of counts, a returncode, and the path

    """
    count = {"gzip": 0, "fai": 0}
    if (not leave_compressed
            and name_ends(name, FASTA_ENDINGS, addition='.gz')):
        if verbose >= 2:
            sys.stderr.write(name + "\n")
        complete_p = subprocess.run(["gunzip", "-f", os.path.join(path, name)])
        returncode = returncode ^ complete_p.returncode
        name = name[0:len(name) - 3]
        if not complete_p.returncode:
            count["gzip"] += 1
    if (name_ends(name, FASTA_ENDINGS)):
        fa_file = os.path.join(path, name)
        if verbose >= 2:
            sys.stderr.write(name + "\n")
        complete_p = subprocess.run([samtools_path, "faidx", fa_file])
        returncode = returncode ^ complete_p.returncode
        if not complete_p.returncode:
            count["fai"] += 1
    return count, returncode, path


def fai_fail_message(returncode, path):
    """
    Print out an error message for when faidx or when gzip fails.

    Parameters
    ----------
    returncode: int
        A 0 if the process exited normally.  Something else otherwise.
    path: str
        A string representing the location of the genome accession's files.

    Returns
    -------
    None

    """
    if returncode:
        print("Making the fai file for {} failed with return code {}.".format(
            path, returncode),
              file=sys.stderr)


def make_fai(root_directory, processes=1, leave_compressed=False, verbose=0):
    """
    Make fai files for all fasta files under the root_directory.
    Decompress fasta files.

    Parameters
    ----------
    root_directory: str
        A string representing the location of the root directory to search.
    processes: int
        The number of processes to use.  This allows for parallelism.
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
    samtools_path_in = SAMTOOLS + "samtools"
    samtools_path = which(samtools_path_in)
    if not samtools_path:
        raise FileNotFoundError(samtools_path_in)
    if processes > 1:
        pool = Pool(processes=processes)
        pd_list = []
        for path, _, files in os.walk(root_directory):
            for name in files:
                pd_list.append(pool.apply_async(remove_fai, args=(path, name)))
        for pd in pd_list:
            pd.get()
        pd_list = []
        for path, _, files in os.walk(root_directory):
            returncode = 0
            for name in files:
                pd_list.append(
                    pool.apply_async(make_fai_core,
                                     args=(path, name, samtools_path,
                                           returncode, leave_compressed,
                                           verbose)))
        for pd in pd_list:
            count_local, returncode, path = pd.get()
            counter += 1
            if returncode is not None:
                fai_fail_message(returncode, path)
                _count_add(count_local, count)
            if counter % VERBOSE_COUNTER == 0 and verbose >= 1:
                sys.stderr.write("Processed: {}\n".format(counter))
                sys.stderr.flush()
    else:
        for path, _, files in os.walk(root_directory):
            count_local, returncode, _ = make_fai_individual(
                path, files, samtools_path, leave_compressed, verbose)
            counter += 1
            if counter % VERBOSE_COUNTER == 0 and verbose >= 1:
                sys.stderr.write("Processed: {}\n".format(counter))
                sys.stderr.flush()
            if returncode is not None:
                fai_fail_message(returncode, path)
                _count_add(count_local, count)
    sys.stderr.flush()
    return count
