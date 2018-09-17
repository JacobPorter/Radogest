###
Calculates one to all mash distance on selected genomes under a species
###

import pickle
import argparse
import json
from math import sqrt
import time
from joblib import Parallel, delayed
import subprocess
import random
import os
import csv

next_taxon_level = {"species": "genus",
                    "genus": "family",
                    "family": "order",
                    "order": "class",
                    "class": "phylum",
                    "phylum": "kingdom",
                    "kingdom": "domain",
                    "domain": "root"}

def main(map_args):
    genome_file = map_args.genome_file
    genomes = pickle.load(open(genome_file, "rb"))
    #print(json.dumps(genomes, indent=4))
    mash_path = map_args.mash
    N = map_args.N
    for keys, vals in genomes['species'].items():
        query_gen = random.sample(list(genomes['species'][keys]), 1)[0]
        gens = []
        gens = list(genomes['species'][keys])
        for items in gens:
            subprocess.check_call([mash_path, "dist",
                                os.path.join(str(genomes['species'][keys][query_gen]), str(query_gen)),
                                os.path.join(str(genomes['species'][keys][items]), str(items))])

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_file', help='Path to sampled\
                        genomes pickle file')
    parser.add_argument('--mash', help='Path to mash')
    parser.add_argument('--N', help='number of genomes at species')
    map_args = parser.parse_args()
    main(map_args)
