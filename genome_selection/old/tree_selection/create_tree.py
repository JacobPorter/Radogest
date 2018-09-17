import ete3
from ete3 import NCBITaxa
import argparse
import pickle
import os
import json

ncbi = NCBITaxa()

next_taxon_level = {"root": "domain",
                    "domain": "kingdom",
                    "kingdom": "phylum",
                    "phylum": "class",
                    "class": "order",
                    "order": "family",
                    "family": "genus",
                    "genus": "species"}

def make_levels(start_rank, stop_rank):
    levels_list = []
    level = start_rank
    if start_rank == stop_rank:
        return levels_list.append(level)
    else:
        while next_taxon_level[level]!=stop_rank:
            levels_list.append(level)
            level = next_taxon_level[level]
        levels_list.extend((level, stop_rank))
    return levels_list


taxonomy_tree = dict()

def child_exists(root_node, root_rank, final_rank, ranks):
    child_list = []
    parent2child = dict()
    species_list = []
    next_level = next_taxon_level[root_rank]
    next_level_taxa = ranks[str(root_rank)][root_node]
    if next_level == 'species':
        return list(next_level_taxa)
    else:
        for taxa in next_level_taxa:
            if len(ranks[next_level][taxa]) == 0:
                continue
            child_list.append(taxa)
    #parent2child[root_node] = child_list
    #for cl in child_list:
        #print(cl)
    return child_list

def main(map_args):
    child_dict = dict()
    query_id = map_args.taxon_id
    queried_rank = map_args.root_taxon_level
    final_rank = map_args.final_taxon_level
    RANK_DIR = map_args.rank_dir
    ranks = pickle.load(open(os.path.join(RANK_DIR, 'ranks.pck'), 'rb'))

    levels = make_levels(queried_rank, final_rank)
    for current_level in levels[:-1]:
        taxid_list = []
        if current_level == queried_rank:
            if current_level not in taxonomy_tree.keys():
                taxonomy_tree[current_level] = {query_id[0]: child_exists(query_id[0], queried_rank, final_rank, ranks)}

        else:
            parent_level = [k for k, v in next_taxon_level.items() if v == current_level]
            for k in taxonomy_tree[parent_level[0]].keys():

                taxid_list = taxonomy_tree[parent_level[0]][k]
                for items in taxid_list:
                    #import pdb; pdb.set_trace()
                    if current_level not in taxonomy_tree.keys():
                        taxonomy_tree[current_level] = {items: child_exists(items, current_level, final_rank, ranks)}
                    else:
                        taxonomy_tree[current_level][items] = child_exists(items, current_level, final_rank, ranks)

    print(json.dumps(taxonomy_tree, indent=4))
    with open('taxonomy_tree.pickle', 'wb') as handle:
        pickle.dump(taxonomy_tree, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxon_id', type=int, nargs=1,
                        help='Taxon ID of the root node of the tree')
    parser.add_argument('--root_taxon_level', type=str,
                        help='Taxon level of the root node')
    parser.add_argument('--final_taxon_level', type=str,
                        help='Taxon level of the leaves')
    parser.add_argument('--rank_dir',
                        help='Path of the ranks pickle datastructure')

    map_args = parser.parse_args()
    main(map_args)
