#!/usr/bin/env python
"""
A utility to add the genome sampling strategy represented in the tree selection
approach into the random sampler.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""
import argparse
import datetime
import pickle
import sys
import pprint
import json
import re


def set_all_genomes(random_index, boolean=False):
    """
    Set all genome inclusions to the given boolean value
    in the random sampler index.

    Parameters
    ----------
    random_index: dict
        The random sampler index dictionary.
    boolean: bool
        A boolean to set all genome inclusions to.

    Returns
    -------
    random_index: dict
        The same random sampler index dictionary with all genome inclusions
        set to the given boolean value.

    """
    taxids_dict = random_index['taxids']
    for taxid in taxids_dict:
        for accession in taxids_dict[taxid]:
            taxids_dict[taxid][accession] = boolean
    return random_index


def integrate_into_sampler(random_index, tree_index):
    """
    Take the genome selection information from the tree index and update the
    random index to use that for its genome selection.

    Parameters
    ----------
    random_index: dict
        The random sampler index dictionary.
    tree_index: dict
        The tree selection index dictionary.

    Returns
    -------
    (random_index, tree_index): (dict, dict)
        The random sampler index and the tree sampler index.

    """
    random_index = set_all_genomes(random_index, False)
    accession_matcher = re.compile('[A-Z]{3}_[0-9]+.[0-9]')
    for rank in tree_index:
        for taxid in tree_index[rank]:
            for filename in tree_index[rank][taxid]:
                accession = accession_matcher.match(filename)
                if accession:
                    try:
                        random_index['taxids'][taxid][accession] = True
                    except KeyError:
                        print("The accession {} for taxid {} was not found in "
                              "the random index.".format(accession, taxid),
                              file=sys.stderr)
                else:
                    print("An accession for taxid {} at rank {} was not "
                          "found in the genome file name {}"
                          .format(taxid, rank, filename), file=sys.stderr)


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('A utility to integrate the '
                                                  'tree genome sampler into '
                                                  'the random sampler.  '
                                                  'This takes a random '
                                                  'sampler index and a tree '
                                                  'selection index.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument('random_index', type=str,
                        help=('The random sampler index object.'),
                        default="./index.pck")
    parser.add_argument('tree_index', type=str,
                        help=('The tree selection index object.'),
                        default='./taxid2genomes.pck')
    parser.add_argument('--output', '-o', type=str,
                        help=('The location and file name for the output '
                              'random sampler index'),
                        default='./random_tree_index')
    args = parser.parse_args()
    print("Modifying the random sampler index with the selection in "
          "the tree index.", file=sys.stderr)
    random_index, _ = integrate_into_sampler(pickle.load(
        open(args.random_index, 'rb')),
        pickle.load(open(args.tree_index, 'rb')))
    print("Finished updating the random index.", file=sys.stderr)
    pickle.dump(random_index, open(args.output + 'pck', 'wb'))
    pprint.pprint(random_index, open(args.output + '.txt', 'w'))
    json.dump(random_index, open(args.output + '.json', 'w'))
    tock = datetime.datetime.now()
    print('The process took time: {}'.format(tock - tick), file=sys.stderr)


if __name__ == '__main__':
    main()
