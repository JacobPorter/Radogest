"""
Produce a hierarchical tree (e.g. UPGMA tree)
from a distance matrix of genomes.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""
import subprocess
from multiprocessing import Manager, Pool, Process

from scipy.cluster.hierarchy import linkage, fcluster 
from scipy.spatial.distance import pdist

