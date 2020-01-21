"""
Compute Mash sketches for genomes in a directory.
This can be used to create distances and cluster genomes.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""

import os
import subprocess
from multiprocessing import Pool

from library.faidx import FASTA_ENDINGS, name_ends

MASH_LOC = "mash"


def sketch(genome_file, output, k, s, sketch_prog=MASH_LOC):
    """
    Produce a sketch of a genome.

    Parameters
    ----------
    genome_file: str
        The location of a fasta genome file.
    output: str
        The prefix for the output file.
    k: int
        The k-mer size for the sketch.
    s: int
        The sketch size for the sketch.
    sketch_prog: str
        The location of the sketch program executable. e.g. mash

    Returns
    -------
    Return the return output from subprocess.run

    Examples
    --------
    Example of running mash:
    mash sketch -k <k> -s <s> -o <outputfile_prefix>
    """
    cp = subprocess.run([
        sketch_prog, "sketch", "-k",
        str(k), "-s",
        str(s), "-o", output, genome_file
    ], capture_output=True)
    output = cp.stdout
    for line in output:
        if line.startswith("Writing"):
            output = line.strip()
            break
    else:
        output = ""
    return cp.returncode, output


def sketch_dir(path, files, k, s):
    """
    For a directory, produce sketches for all fasta files in that directory.

    Parameters
    ----------
    path: str
        The location of a directory.
    files: list
        A list of all of the files in the directory path
    k: int
        The number of k-mers in a sketch (For mash, [1-32])
    s: int
        The sketch size.

    Parameters
    ----------
    int
        Return code.  1 if abnormal.  0 if normal.

    """
    ret_codes = []
    for f in files:
        len_end = name_ends(f, FASTA_ENDINGS, addition="")
        len_end = len_end if len_end else name_ends(
            f, FASTA_ENDINGS, addition=".gz")
        if len_end:
            genome_file = os.path.join(path, f)
            ret_codes.append(sketch(genome_file, genome_file, k, s))
    return 1 if any(ret_codes) else 0


def sketch_root(root_directory, k, s, processes=1):
    """
    Produce sketches for all fasta files under a directory and its subdirectories.

    Parameters
    ----------
    root_directory: str
        The root directory where fasta files are located.
    k: int
        The number of k-mers in a sketch (For mash, [1-32])
    s: int
        The sketch size.
    processes: int
        The number of processes to use for parallelism.
        Each process takes a directory.

    Returns:
    -------
    int
        A return code.  If everything went well, then 0, else non-0.

    """
    pool = Pool(processes=processes)
    pd_list = []
    retcode = 0
    for path, _, files in os.walk(root_directory):
        pd_list.append(pool.apply_async(sketch_dir, args=(path, files, k, s)))
    for pd in pd_list:
        out = pd.get()
        retcode = retcode ^ out[0]
        print(out[1])
    return retcode
