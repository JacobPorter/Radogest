#!/usr/bin/env python
"""
Select genomes based on a selection strategy and a depth-first post order
traversal of the taxonomic tree.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""


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
        selects genomes at each level.
        """
        children = self.tree[taxid]
        levels_visited = 0
        for child in children:
            levels_visited += self.select_genomes(child)
        self.selection_strategy.select(taxid, children)
        return 1 + levels_visited


class TaxTreeListTraversal:
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
        selects genomes at each level.
        """
        children = self.tree[taxid]
        genome_lists = []
        levels_visited = 0
        for child in children:
            levels_child, genomes_child = self.select_genomes(child)
            genome_lists.append(genomes_child)
            levels_visited += levels_child
        return (1 + levels_visited, self.selection_strategy.select(taxid, children, genome_lists))