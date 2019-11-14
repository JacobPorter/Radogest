"""
Construct a taxonomic rank tree hierarchy from the genomes index.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""
import sys
from collections import defaultdict
from random import shuffle

from ete3 import NCBITaxa
from tqdm import tqdm

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
        ret_node, amount = remove_only_children_traversal(
            tree, child, deletion_list, remapping)
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
    tuple: (dict, int, int, int, dict, list<str>)
        1. ranks: The tree dictionary
        2. unclassified_to_remove: A count of the unclassified species removed.
        3. empties_remove: Empty ranks removed.
        4. children_removed: The amount of children removed.
        5. counter_dict: A dictionary of counts.
        6. taxid_list: The non-leaf taxonomic id in a random order.

    """
    rank_list = [
        "superkingdom", "kingdom", "phylum", "class", "order", "family",
        "genus", "species"
    ]
    # What does this do?
    # missing_ranks = defaultdict(list)
    looked_at = {}
    ranks = defaultdict(set)
    print("Creating the initial tree.", file=sys.stderr)
    for accession in tqdm(index['genomes']):
        try:
            taxid = int(index['genomes'][accession]['taxid'])
        except KeyError:
            print("Could not find a key for: {}".format(accession),
                  file=sys.stderr)
            raise KeyError
        if taxid in looked_at:
            continue
        if not index['taxids'][taxid]:
            continue
        lineage_rank = ncbi.get_rank(ncbi.get_lineage(taxid))
        lineage_rank = {
            lineage_rank[key]: key
            for key in lineage_rank if lineage_rank[key] != 'no rank'
        }
        looked_at[taxid] = True
        text_translator = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
        exclude = False
        for text in text_translator.values():
            if (len(lineage_rank) <= 2 or
                    # 'unclassified' in text or  # Include unclassified?
                    'uncultured' in text or 'RNA virus' in text
                    or 'satellite RNA' in text
                    or 'environmental samples' in text):
                exclude = True
                break
        if exclude:
            continue
        for i, rank in enumerate(rank_list):
            if i == len(rank_list) - 1:
                break
            if rank not in lineage_rank:
                continue
            for j in range(i + 1, len(rank_list)):
                if rank_list[j] in lineage_rank:
                    ranks[lineage_rank[rank]].add(lineage_rank[rank_list[j]])
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
    remove_domains = [domain for domain in ranks[1] if not ranks[domain]]
    for domain in remove_domains:
        ranks[1].remove(domain)
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
    shuffle(taxid_list)
    return (ranks, unclassified_to_remove, empties_remove, children_removed,
            counter_dict, taxid_list)
