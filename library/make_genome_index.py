#!/usr/bin/env python
"""
Makes genome location index dictionaries from JSON records
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import datetime
import pickle
import os
import sys
import json
from collections import defaultdict
from tqdm import tqdm
from ete3 import NCBITaxa

# COUNT_THRESHOLD = 1000


def add_index(index, genome_info, verbose=True):
    """
    Add the contents of JSON genome information to an index.

    Parameters
    ----------
    index: dictionary
        a dictionary representing an index of taxonomic ids to genomes
    genome_info: iterable of dictionaries
        each record represents a genome
    verbose: boolean
        if set to True, prints out periodic counts of records processed.

    Returns
    -------
    int
        A count of the records processed.

    """
    ncbi = NCBITaxa()
    count = 0
    for genome_record in tqdm(genome_info):
        group = genome_record[1]
        genome_record = genome_record[0]
        accession = genome_record["assembly_accession"]
        index['genomes'][accession][
            'location'] = '/{section}/{group}/{accession}/'.format(
                section=genome_record['section'],
                group=group, accession=accession)
        for key in genome_record:
            index['genomes'][accession][key] = genome_record[key]
        for taxid in ncbi.get_lineage(int(genome_record["taxid"])):
            if accession not in index['taxids'][taxid]:
                index['taxids'][taxid][accession] = True
        # if verbose and count % COUNT_THRESHOLD == 0:
        #    sys.stderr.write("Processed {} records.\n".format(count))
        count += 1
    return count


def main():
    """Parse arguments and dispatches index creation."""
    parser = argparse.ArgumentParser(description=('Make genome location and '
                                                  'information indexes from '
                                                  'json files.  The index '
                                                  'associates taxonomic ids '
                                                  'with genomes, and genome '
                                                  'accessions are associated '
                                                  'with information, which '
                                                  'includes the location of '
                                                  'the genome in the file '
                                                  'system.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("json_file_dir", help=("The location of the json "
                                               "file or a directory "
                                               "containing json files."))
    parser.add_argument("--index", "-i", type=str,
                        help=("An existing index location. "
                              "If empty, a new index is created."),
                        default='')
    parser.add_argument("--output", "-o",
                        help="The location to write the index file.",
                        default="index_initial.pck")
    parser.add_argument("-v", "--verbose", action='store_true', default=False)
    args = parser.parse_args()
    now = datetime.datetime.now()
    index = ({'taxids':defaultdict(dict), 'genomes':defaultdict(dict)} if
             not args.index else pickle.load(open(args.index, 'rb')))
    if os.path.isfile(args.json_file_dir):
        count = add_index(index, json.load(open(args.json_file)), args.verbose)
    elif os.path.isdir(args.json_file_dir):
        onlyfiles = [f for f in os.listdir(args.json_file_dir)
                     if os.path.isfile(os.path.join(args.json_file_dir, f))]
        onlyfiles = [f for f in onlyfiles if f.endswith('.json')]
        count = 0
        for f_f in onlyfiles:
            sys.stderr.write("Processing file {}.\n".format(f_f))
            my_path = os.path.join(args.json_file_dir, f_f)
            count += add_index(index, json.load(open(my_path)), args.verbose)
    else:
        sys.stderr.write("The json parameter was not a directory and was not "
                         "a file.  Please to check.\n")
    sys.stderr.write("There were {} records processed.\n".format(count))
    sys.stderr.write("Dumping the index to: {}\n".format(args.output))
    pickle.dump(index, open(args.output, 'wb'))
    json.dump(index, open(args.output + ".json", 'w'), indent=4)
    sys.stderr.write("The process took {} time.\n".format(
        datetime.datetime.now() - now))

if __name__ == '__main__':
    main()
