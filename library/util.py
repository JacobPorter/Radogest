"""
Some utility files.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import math

def one_loop_stats(iterable, num_func, prefix=""):
    """
    Calcualte the mean, standard deviation, 
    count, sum, max, and min in one loop.
    
    Parameters
    ----------
    iterable: iterable
        An object that can be iterated through such as a list or an open file.
    name_func: function
        A function that when given an input line from iterable 
        will return a number.  This allows this function to be more generic.
    prefix: str
        A string prefix to add to each statistic.
        
    Returns
    -------
    dict
        A dictionary of the statistics.
    
    """
    s = 0.0
    s_2 = 0.0
    n = 0
    mx = -float("inf")
    mn = float("inf")
    for line in iterable:
        ai = num_func(line)
        if ai > mx:
            mx = ai
        if ai < mn:
            mn = ai
        s += ai
        s_2 += ai * ai
        n += 1
    mean = s / n
    var = s_2 / n - mean * mean
    return {prefix + "mean": mean, prefix + "std": math.sqrt(var),
            prefix + "count": n, prefix + "sum": s, 
            prefix + "max": mx, prefix + "min": mn}


def faidx_length(line):
    """
    Get the base length for a contig from a fai file.
    
    Parameters
    ----------
    line: str
        A fai file line.
        
    Returns
    -------
    int
        The number of bases in the contig.
    
    """
    return float(line.split('\t')[1])