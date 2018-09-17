
"""Tree selection"""

from ete3 import NCBITaxa
import argparse
import pickle
import os
import json
import math
from tqdm import tqdm
import random
import time

ncbi = NCBITaxa()
hierarchy = dict()

next_taxon_level = {"species": "genus",
                    "genus": "family",
                    "family": "order",
                    "order": "class",
                    "class": "phylum",
                    "phylum": "kingdom",
                    "kingdom": "domain",
                    "domain": "root"}

hierarchy = dict()

def get_species2genomes(index, taxonomy_tree, data_dir,
                        subselect_species, n):
    if len(subselect_species) is not 0:
        slst = []
        slst = subselect_species[:]
    else:
        slst = []
        for k, v in taxonomy_tree['genus'].items():
            for items in v:
                slst.append(items)
    print('\n Creating dataset at taxonomy - species', '\n')
    taxid2genomes = dict()
    for species in tqdm(slst):
        genomes = dict()
        try:
            gens = random.sample(index['taxids'][species].keys(), n)
        except:
            gens = random.sample(index['taxids'][species].keys(),
                        len(index['taxids'][species].keys()))
        for items in gens:
            gen_dir = os.path.join(data_dir,
                                   index['genomes'][items]['location'][1:])
            for file in os.listdir(gen_dir):
                if file.endswith('.fna'):
                    genomes[file] = (str(species), gen_dir)
        taxid2genomes[species] = genomes
    hierarchy['species'] = taxid2genomes

    # print(json.dumps(hierarchy, indent=4))
    # with open('taxid2genomes.pickle', 'wb') as handle:
    #     pickle.dump(hierarchy, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Remove root_rank and use ncbi
# Remove ranks from the tree and use ncbi if needed
# Update the index data structure instead of creating a new hierarchy
def sampleGenomes(root_rank, n, taxonomy_tree, leaf_rank, species_select,
                 data_dir, index, hierarchy):
    if root_rank == leaf_rank:
        return
    next_level_taxa = []
    taxid2genomes = dict()
    for k, v in hierarchy[leaf_rank].items():
        next_level_taxa.append(k)
    print('\n Creating dataset at taxonomy -', hierarchy[leaf_rank], '\n')
    for keys, vals in tqdm(taxonomy_tree[next_taxon_level[leaf_rank]].items()):
        genomes = dict()
        child = []
        if species_select == 0:
            child = taxonomy_tree[next_taxon_level[leaf_rank]][keys]
        else:
            child = []
            for items in tqdm(vals):
                if items in next_level_taxa:
                    child.append(items)
            if len(child) > 0:
                prob = 1.0/len(child)
                for samples in child:
                    genomes_from_child = math.ceil(prob*n)
                    total = len(hierarchy[leaf_rank][samples])
                    if genomes_from_child > total:
                        for k, v in hierarchy[leaf_rank][samples].items():
                            genomes[k] = (str(samples), v[1])
                    else:
                        indices_to_choose = random.sample(list(
                            hierarchy[leaf_rank][samples]),
                            genomes_from_child)
                        for items in indices_to_choose:
                            genomes[items] = (str(samples),
                                              hierarchy[
                                                  leaf_rank][
                                                      samples][items][1])
        if len(genomes) > n:
            for i in random.sample(genomes.keys(), len(genomes)-n):
                del genomes[i]
        if len(genomes) > 0:
            taxid2genomes[keys] = genomes
    hierarchy[next_taxon_level[leaf_rank]] = taxid2genomes
    sampleGenomes(root_rank, n, taxonomy_tree, next_taxon_level[leaf_rank],
                  species_select, data_dir, index, hierarchy)


def count_kmers(seq, k, id):
    set_kmers = set([])
    num_kmers = len(seq)-k+1
    for i in range(num_kmers):
        kmer = seq[i:i+k]
        set_kmers.add(kmer)
    return set_kmers

def main(map_args):
    start_time = time.time()
    query_id = map_args.taxon_id[0]
    query_rank = map_args.root_taxon_level
    final_rank = map_args.final_taxon_level
    data_dir = map_args.data_dir
    tree_file = map_args.tree
    n = map_args.genomes
    index = pickle.load(open(os.path.join(data_dir, 'index.pck'), 'rb'))
    # output_directory = map_args.output_directory
    # species_select = map_args.species_select
    # train_set = map_args.train_set
    # fragment_len = map_args.fragment_len
    # train_fraction = map_args.train_fraction
    subselect_species = map_args.subselect_species
    taxonomy_tree = pickle.load(open(tree_file, 'rb'))

    # print(subselect_species)
    get_species2genomes(index, taxonomy_tree, data_dir,
                        subselect_species, n)
    # print("--- %s seconds ---" % (time.time() - start_time))
    leaf_rank = 'species'
    sampleGenomes(query_rank, n, taxonomy_tree, leaf_rank,
                  43, data_dir, index, hierarchy)
    print(json.dumps(hierarchy, indent=4))
    with open('taxid2genomes_tuple.pickle', 'wb') as handle:
        pickle.dump(hierarchy, handle, protocol=pickle.HIGHEST_PROTOCOL)
    end_time = time.time()
    print(end_time-start_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxon_id', type=int, nargs=1,
                        help='Taxon ID of the root node of the tree')
    parser.add_argument('--root_taxon_level', type=str,
                        help='Taxon level of the root node')
    parser.add_argument('--final_taxon_level', type=str,
                        help='Taxon level of the leaves')
    parser.add_argument('--data_dir',
                        help='Path to Genomes directory')
    parser.add_argument('--tree',
                        help='Path of the taxonomy tree pickle file')
    parser.add_argument('--genomes', type=int,
                        help='number of genomes at each node')
    #parser.add_argument('--output_directory',
    #                    help='Path to storing genomes')
    #parser.add_argument('--train_set', type=str,
    #                    help='directory for storing original genomes')
    #parser.add_argument('--fragment_len', type=int,
    #                    help='length of fragments')
    #parser.add_argument('--train_fraction', type=float,
    #                    help='fraction of training samples')
    parser.add_argument('--subselect_species', nargs="*",\
                        type=int, default=[])

    map_args = parser.parse_args()
    main(map_args)
