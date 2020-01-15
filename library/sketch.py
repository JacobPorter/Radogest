"""
Compute Mash sketches for genomes in a directory.
This can be used to create distances and cluster genomes.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""

import subprocess
import os

from library.faidx import FASTA_ENDINGS, name_ends

from multiprocessing import Manager, Pool, Process

MASH_LOC = "mash"

def sketch(genome_file, output, k, s, sketch_prog=MASH_LOC):
    """
    Produce a sketch of a genome.
    
    Parameters
    ----------
    genome_file: str
        The location of a fasta genome file.
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
    cp = subprocess.run([sketch_prog, "sketch", "-k", str(k), "-s",
                           str(s), "-o", output, genome_file])
    return cp.returncode


def sketch_dir(path, files, k, s):
    """
    
    """
    ret_codes = []
    for f in files:
        len_end = name_ends(f, FASTA_ENDINGS, addition="")
        len_end = len_end if len_end else name_ends(f, FASTA_ENDINGS, addition=".gz")
        if len_end:
            genome_file = os.path.join(path, f)
            ret_codes.append(sketch(genome_file, 
                                    genome_file[0:len(genome_file)-len_end-1],
                                    k,
                                    s))
    return 1 if any(ret_codes) else 0
    
    
def sketch_root(root_directory, k, s, processes=1):
    """
    
    """
    pool = Pool(processes=processes)
    pd_list = []
    retcode = 0
    for path, _, files in os.walk(root_directory):
        pd_list.append(pool.apply_async(sketch_dir, args=(path, files, k, s)))
    for pd in pd_list:
        retcode = retcode ^ pd.get()
    return retcode
        