#!/usr/bin/env python
"""
Select genomes based on a selection strategy and a depth-first post order
traversal of the taxonomic tree.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import datetime
import sys
import argparse
import pickle
import json
from library.genome_selection.strategy import ProportionalRandom, QualitySortingTree, QualitySortingLeaf
from library.genome_selection.strategy import AllGenomes, MinHashTree
from library.genome_selection.strategy import EXCLUDED_GENOMES


class StrategyNotFound(Exception):
    """An exception to raise when a genome selection strategy does not exit."""

    def __init__(self, *args, **kwargs):
        """Pass to the super class."""
        Exception.__init__(self, *args, **kwargs)


class TaxTreeTraversal:
    """
    Traverse the taxonomic tree and do genome selection with a
    depth-first post order traversal.
    """

    def __init__(self, tree, selection_strategy):
        """Save the input arguments."""
        self.tree = tree
        self.selection_strategy = selection_strategy

    def select_genomes(self, taxid):
        """
        Perform depth-first post order traversal of the taxonomic tree that
        selects genomes at leach level.
        """
        children = self.tree[taxid]
        levels_visited = 0
        for child in children:
            levels_visited += self.select_genomes(child)
        self.selection_strategy.sample(taxid, children)
        return 1 + levels_visited


def main():
    """Parse arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Downselect genomes '
                                                  'according to some '
                                                  'strategy.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("taxid", type=int, help=('The taxonomic id for '
                                                 'the root.'))
    parser.add_argument("--index", '-i', type=str,
                        help=('The pickle file of the genomes index.'),
                        default='./index.pck')
    parser.add_argument("--tree", '-t', type=str,
                        help=('The pickle file containing the '
                              'taxonomic tree.'),
                        default='./tree.pck')
    parser.add_argument("--strategy", '-s', action='store',
                        choices=['PR', 'QST', 'QSL', 'MHT', 'AG'],
                        help=('Choose the genome selection strategy.  '
                              'The choices are: ProportionalRandom (PR), '
                              'QualitySortingTree (QST), QualitySortingLeaf '
                              '(QSL), MinHashTree (MHT), AllGenomes (AG)'),
                        default='PR')
    parser.add_argument('--sample_amount', '-n', type=int,
                        help=('The number of genomes to sample at each level. '
                              'Does not apply with AllGenomes.'),
                        default=10)
    parser.add_argument('--output', '-o', type=str,
                        help=('The location to store the index with '
                              'the down selected genomes.'),
                        default='./index_down.pck')
    args = parser.parse_args()
    strategy_string = args.strategy.upper()
    index = pickle.load(open(args.index, 'rb'))
    tree = pickle.load(open(args.tree, 'rb'))
    sample_amount = args.sample_amount
    if sample_amount <= 0:
        parser.error('The sample amount needs to be a positive non-zero '
                     'integer.')
    if strategy_string == 'PR':
        strategy = ProportionalRandom(index, sample_amount)
    elif strategy_string == 'QST':
        strategy = QualitySortingTree(index, sample_amount)
    elif strategy_string == 'QSL':
        strategy = QualitySortingLeaf(index, sample_amount)
    elif strategy_string == 'MHT':
        strategy = MinHashTree(index, sample_amount)
    elif strategy_string == 'AG':
        strategy = AllGenomes(index)
    else:
        raise StrategyNotFound()
    traversal = TaxTreeTraversal(tree, strategy)
    levels_visited = traversal.select_genomes(args.taxid)
    for accession in EXCLUDED_GENOMES:
        print("WARNING: {} excluded because {}.".format(
            accession, EXCLUDED_GENOMES[accession]), file=sys.stderr)
    print("Levels visited: {}".format(levels_visited), file=sys.stderr)
    print("Creating a pickled index.", file=sys.stderr)
    pickle.dump(index, open(args.output, 'wb'))
    print("Creating a json index", file=sys.stderr)
    json.dump(index, open(args.output + '.json', 'w'), indent=4)
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
