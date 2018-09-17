#!/usr/bin/env python
"""
Construct a taxonomic rank tree hierarchy from the genomes index.
:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import datetime
import pickle
import sys
import json
from tqdm import tqdm
from collections import defaultdict
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def count_levels(tree, roots):
    """
    Count the number of taxid that are not leaves under the taxid under roots.

    Parameters
    ----------
    tree: dict
        The tree hierarchy.
    roots: list
        A list of taxonomic ids to search under.

    Returns
    -------
    dict
        A dictionary indexed on a taxonomic id root from roots.
        The first number is the number of taxonomic ids under root.
        The second number is the number of taxonomic ids under root that
        have only a single child.

    """
    counter_dict = {root: [0, 0] for root in roots}
    taxid_list = []
    for taxid in tree:
        lineage = ncbi.get_lineage(taxid)
        for root in roots:
            if root in lineage:
                taxid_list.append(taxid)
                counter_dict[root][0] = counter_dict[root][0] + 1
                if len(tree[taxid]) == 1:
                    counter_dict[root][1] = counter_dict[root][1] + 1
                break
    return counter_dict, taxid_list


def remove_only_children(tree, roots):
    """
    Remove nodes with only one child.  Reconnect parent nodes to their
    new children.

    Parameters
    ----------
    tree: dict
        The taxonomic tree.
    roots: iterable
        The root nodes in the tree to traverse down.

    Returns
    -------
    int
        The number of nodes removed.

    """
    nodes_deleted = 0
    for root in roots:
        deletion_list = []
        remapping = defaultdict(list)
        remove_only_children_traversal(tree, root, deletion_list, remapping)
        for taxid in remapping:
            for old_child, new_child in remapping[taxid]:
                tree[taxid].remove(old_child)
                tree[taxid].append(new_child)
        for taxid in deletion_list:
            del tree[taxid]
        nodes_deleted += len(deletion_list)
    return nodes_deleted


def remove_only_children_traversal(tree, node, deletion_list, remapping):
    """
    Traverse the taxonomic tree (depth-first search) to identify nodes
    with only one child.  Gives a remapping from parents to children.

    Parameters
    ----------
    tree: dict
        The taxonomic tree.
    node: int
        The rooth node to start searching.
    deletion_list: list
        A list to store taxonomic ids to delete.
    remapping: defaultdict(list)
        A mapping between parents and their [(old_child, new_child), ..., ]

    Returns
    -------
    (int, int)
        A taxonomic node and its number of children.

    """
    children = tree[node]
    if len(children) == 0:
        return (node, 0)
    for child in children:
        ret_node, amount = remove_only_children_traversal(tree, child,
                                                          deletion_list,
                                                          remapping)
        if ret_node != child and len(children) != 1:
            remapping[node].append((child, ret_node))
    if len(children) == 1:
        deletion_list.append(node)
        return ret_node, amount
    else:
        return (node, len(children))


def check_leaves_species(taxid, tree):
    """
    Validate that the leaves of the tree are all species.

    Parameters
    ----------
    taxid: int
        An integer representing a taxonomic id.
    tree: dict
        A dictionary representing the taxonomic tree.

    Returns
    -------
    list
        A list of taxonomic ids that are leaves and that are not species.

    """
    if not tree[taxid] and 'species' not in ncbi.get_rank([taxid])[taxid]:
        return [taxid]
    elif not tree[taxid]:
        return []
    list_of_taxid = []
    for child in tree[taxid]:
        list_of_taxid + check_leaves_species(child, tree)
    return list_of_taxid


def make_tree(index, verbose=False):
    """
    Construct a tree hierarchy using conventional ranks from the genomes index.

    Parameters
    ----------
    index: dict
        The all genomes index
    verbose: boolean
        Print additional output.

    Returns
    -------
    (dict, dict)
        1. A dictionary object that represents the subranks below each rank.
        2. A dictionary that gives a list of missing ranks for species.

    """
    rank_list = ["superkingdom", "kingdom", "phylum", "class", "order",
                 "family", "genus", "species"]
    missing_ranks = defaultdict(list)
    looked_at = {}
    # ranks = {"superkingdom":defaultdict(set), "kingdom":defaultdict(set),
    #          "phylum":defaultdict(set), "class":defaultdict(set),
    #          "order":defaultdict(set), "family":defaultdict(set),
    #          "genus":defaultdict(set)}
    ranks = defaultdict(set)
    print("Creating the initial tree.", file=sys.stderr)
    for accession in tqdm(index['genomes']):
        taxid = int(index['genomes'][accession]['taxid'])
        if taxid in looked_at:
            continue
        lineage_rank = ncbi.get_rank(ncbi.get_lineage(taxid))
        lineage_rank = {lineage_rank[key]: key for key in lineage_rank if
                        lineage_rank[key] != 'no rank'}
        looked_at[taxid] = True
        text_translator = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
        exclude = False
        for text in text_translator.values():
            if (len(lineage_rank) <= 2 or
                    # 'unclassified' in text or  # Include unclassified?
                    'uncultured' in text or
                    'RNA virus' in text or
                    'satellite RNA' in text or
                    'environmental samples' in text):
                exclude = True
                break
        if exclude:
            continue
        for i, rank in enumerate(rank_list):
            if i == len(rank_list) - 1:
                break
            if rank not in lineage_rank:
                continue
            for j in range(i+1, len(rank_list)):
                if rank_list[j] in lineage_rank:
                    ranks[lineage_rank[rank]].add(
                        lineage_rank[rank_list[j]])
                    break
    for taxid in ranks:
        ranks[taxid] = list(ranks[taxid])
    unclassified = []
    print("Removing species with 'unclassified' in their lineage.",
          file=sys.stderr)
    for parent in tqdm(ranks):
        for taxid in ranks[parent]:
            if 'species' in ncbi.get_rank([taxid])[taxid]:
                text_translator = ncbi.get_taxid_translator(
                    ncbi.get_lineage(taxid))
                exclude = False
                for text in text_translator.values():
                    if 'unclassified' in text.lower():
                        if verbose:
                            print(text_translator, file=sys.stderr)
                        exclude = True
                        break
                if exclude:
                    unclassified.append(taxid)
    unclassified_to_remove = len(unclassified)
    for taxid in ranks:
        for species in unclassified:
            if species in ranks[taxid]:
                ranks[taxid].remove(species)
    # Add the root node for eukaryotes, bacteria, viruses, archaea
    ranks[1] = [2759, 2, 10239, 2157]
    print("Removing nodes with only one child.", file=sys.stderr)
    children_removed = remove_only_children(ranks, [1, 12884])
    print("Removing empties.", file=sys.stderr)
    empties_remove = 0
    to_delete = []
    for taxid in ranks:
        if not ranks[taxid]:
            to_delete.append(taxid)
            empties_remove += 1
    for taxid in to_delete:
        del ranks[taxid]
    print("Counting the levels and producing the taxonomic id list.",
          file=sys.stderr)
    counter_dict, taxid_list = count_levels(ranks, ranks[1])
    taxid_list = [1] + taxid_list
    return (ranks, missing_ranks,
            unclassified_to_remove,
            empties_remove, children_removed, counter_dict, taxid_list)


def main():
    """Parse arguments and prints output."""
    parser = argparse.ArgumentParser(description=('Construct an object '
                                                  'representing a tree of '
                                                  'taxonomic ranks.  '
                                                  'Print a list of '
                                                  'non leaf non viroid '
                                                  'taxids to stdout.'),
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("index", help=("The location of the pickled "
                                       "index file."))
    parser.add_argument("--output", "-o", type=str,
                        help=("The file to write the ranks tree to."),
                        default="./tree")
    parser.add_argument("--output_missing", "-m", type=str,
                        help=("The file location to write missing rank "
                              "information too."),
                        default="./missing.pck")
    parser.add_argument("-v", "--verbose", action='store_true', default=False)
    args = parser.parse_args()
    now = datetime.datetime.now()
    print("Creating the taxonomic tree.", file=sys.stderr)
    (ranks, missing_ranks, unclassified,
     empties, children_removed,
     counter_dict, taxid_list) = make_tree(pickle.load(open(args.index, 'rb')))
    pickle.dump(ranks, open(args.output + '.pck', 'wb'))
    json.dump(ranks, open(args.output + '.json', 'w'), indent=4)
    pickle.dump(missing_ranks, open(args.output_missing, 'wb'))
    print("Unclassified: {}, empties: {}, children removed: {}".
          format(unclassified, empties, children_removed),
          file=sys.stderr)
    print("Counts: {}".format(counter_dict), file=sys.stderr)
    total_taxid = 0
    singleton_taxid = 0
    for taxid in counter_dict:
        total_taxid += counter_dict[taxid][0]
        singleton_taxid += counter_dict[taxid][1]
    print("Total taxid: {}, Singletons: {}, Multiclass:{}".
          format(total_taxid, singleton_taxid, total_taxid - singleton_taxid),
          file=sys.stderr)
    leaves_not_species = check_leaves_species(1, ranks)
    if leaves_not_species:
        print("There were leaves that were not species:\n{}".format(
            leaves_not_species), file=sys.stderr)
    for taxid in taxid_list:
        print(taxid)
    sys.stderr.write("The process took {} time.\n".format(
        datetime.datetime.now() - now))


if __name__ == '__main__':
    main()
