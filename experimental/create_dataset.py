import pickle
import argparse
import json
from math import sqrt
import time
from joblib import Parallel, delayed

def main(map_args):
    genome_file = map_args.genome_file
    genomes = pickle.load(open(genome_file, "rb"))
    print(json.dumps(genomes, indent=4))
    



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_file', help='Path to sampled\
                        genomes pickle file')
    map_args = parser.parse_args()
    main(map_args)
