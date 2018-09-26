import ete3
from ete3 import NCBITaxa
import argparse
import pickle
import os
import json
import math
from tqdm import tqdm
from shutil import copyfile
import random
import Bio
from Bio import SeqIO
import time

ncbi = NCBITaxa()

next_taxon_level = {"species": "genus",
                    "genus": "family",
                    "family": "order",
                    "order": "class",
                    "class": "phylum",
                    "phylum": "kingdom",
                    "kingdom": "domain",
                    "domain": "root"}


genomes = dict()

def get_species2genomes(index, taxonomy_tree, output_directory, data_dir,
                        subselect_species, train_set, train_fraction, fragment_len):
    os.makedirs(os.path.join(output_directory, 'species'))
    if len(subselect_species) is not 0:
        slst = []
        slst = subselect_species[:]
    else:
        slst = []
        for k, v in taxonomy_tree['genus'].items():
            for items in v:
                slst.append(items)

    directory_for_ml_training = 'ml_training'
    directory_for_db_creation = 'db_creation'
    directory_for_ml_testing  = 'ml_test'

    print('\n Creating dataset at taxonomy - species', '\n')
    for species in tqdm(slst):
        new_output_directory = os.path.join(output_directory, 'species')
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_training, str(species)))
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_testing, str(species)))
        skip_length = fragment_len
        start_pos = 0

        if train_set == 'True':
            os.makedirs(os.path.join(new_output_directory, directory_for_db_creation, str(species)))

        with open(os.path.join(new_output_directory, directory_for_ml_training,
                str(species), str(species) + ".fasta"), 'a') as fasta_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_training,
                str(species), str(species) + ".taxid"), 'a') as taxid_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                str(species), str(species) + ".fasta"), 'a') as fasta_handle_test,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                str(species), str(species) + ".taxid"), 'a') as taxid_handle_test:

            if train_set == 'True':
                fasta_db_handle = open(os.path.join(new_output_directory, directory_for_db_creation,
                                        str(species), str(species) + ".fasta"), 'a')
                taxid_db_handle = open(os.path.join(new_output_directory, directory_for_db_creation,
                                        str(species), str(species) + ".taxid"), 'a')


            for keys, vals in tqdm(index['taxids'][species].items()):
                gen_dir = os.path.join(data_dir, index['genomes'][keys]['location'][1:])
                for file in os.listdir(gen_dir):
                    if file.endswith('.fna'):
                        kmer_lst = []
                        for fasta_record in SeqIO.parse(os.path.join(gen_dir, file), 'fasta'):
                            kmer_lst = list(count_kmers(fasta_record.seq, 100, fasta_record.id))
                            test_indices = random.sample(range(0, len(kmer_lst)),
                                                        len(kmer_lst)-int(train_fraction*len(kmer_lst)))
                            test_list = [kmer_lst[i] for i in test_indices]
                            # train_list = [item for item in kmer_lst if item not in test_list]
                            train_list = list(set(kmer_lst) - set(test_list))
                            for train, test in zip(train_list, test_list):
                                    fasta_handle_train.write(">" + fasta_record.id + "\n")
                                    fasta_handle_train.write(str(train) + "\n")
                                    taxid_handle_train.write(str(species) + "\n")
                                    fasta_handle_test.write(">" + fasta_record.id + "\n")
                                    fasta_handle_test.write(str(test) + "\n")
                                    taxid_handle_test.write(str(species) + "\n")

                            if train_set == "True":
                                fasta_db_handle.write(">" + fasta_record.id + "\n")
                                fasta_db_handle.write(str(fasta_record.seq) + "\n")
                                taxid_db_handle.write(str(species) + "\n")

            if train_set == "True":
                fasta_db_handle.close()
                taxid_db_handle.close()

def sampleGenomes(root_rank, n, taxonomy_tree, leaf_rank, species_select, output_directory):
    if root_rank == leaf_rank:
        return
    os.makedirs(os.path.join(output_directory, next_taxon_level[leaf_rank]))
    for keys, vals in taxonomy_tree[next_taxon_level[leaf_rank]].items():
        child = []
        if species_select == 0:
            child = taxonomy_tree[next_taxon_level[leaf_rank]][keys]
        else:
            next_dir = os.path.join(output_directory, leaf_rank)
            for items in vals:
                if os.path.exists(os.path.join(next_dir, str(items))):
                    child.append(items)
            if len(child) > 0:
                prob = 1.0/len(child)
                gen_dir = os.path.join(output_directory,
                        next_taxon_level[leaf_rank], str(keys))
                os.makedirs(gen_dir)
                for samples in child:
                    genomes_from_child = math.ceil(prob*n)
                    DIR = os.path.join(output_directory, leaf_rank, str(samples))
                    total = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])
                    if genomes_from_child > total:
                        for files in os.listdir(DIR):
                            if files.endswith('.fna'):
                                copyfile(os.path.join(DIR, files), os.path.join(gen_dir, files))
                    else:
                        # import pdb; pdb.set_trace()
                        filelist = []
                        for files in os.listdir(DIR):
                            if files.endswith('.fna'):
                                filelist.append(files)
                        indices_to_choose = random.sample(range(0, len(filelist)), genomes_from_child)
                        for i, file in enumerate(filelist):
                            if i in indices_to_choose:
                                copyfile(os.path.join(DIR, files), os.path.join(gen_dir, files))

    sampleGenomes(root_rank, n, taxonomy_tree, next_taxon_level[leaf_rank],
                  species_select, output_directory)

def count_kmers(seq, k, id):
    counts = {}
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
    output_directory = map_args.output_directory
    #species_select = map_args.species_select
    train_set = map_args.train_set
    fragment_len = map_args.fragment_len
    train_fraction = map_args.train_fraction
    subselect_species = map_args.subselect_species
    taxonomy_tree = pickle.load(open(tree_file, 'rb'))

    # print(subselect_species)
    lst = get_species2genomes(index, taxonomy_tree, output_directory, data_dir,
                            subselect_species, train_set, train_fraction, fragment_len)

    print("--- %s seconds ---" % (time.time() - start_time))
    leaf_rank = 'species'
    # sampleGenomes(query_rank, n, taxonomy_tree, leaf_rank,
    #               species_select, output_directory)

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
    parser.add_argument('--output_directory',
                        help='Path to storing genomes')
    # parser.add_argument('--species_select', type=str,
    #                     help='number of species to be sub-selected,\
    #                     True if sub_selection is required,\
    #                     False if all species need to be used')
    parser.add_argument('--train_set', type=str,
                        help='directory for storing original genomes')
    parser.add_argument('--fragment_len', type=int,
                        help='length of fragments')
    parser.add_argument('--train_fraction', type=float,
                        help='fraction of training samples')
    parser.add_argument('--subselect_species', nargs="*",\
                        type=int, default=[])

    map_args = parser.parse_args()
    main(map_args)
