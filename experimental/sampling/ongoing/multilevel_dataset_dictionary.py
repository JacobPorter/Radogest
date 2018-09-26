import os
from Bio import SeqIO
import random
import pdb
from ete2 import NCBITaxa
from tqdm import tqdm
import math
import time

"""
script to return tree hierarchy with file names and genomes under it
not based on genomes
"""

# ----------- #
#  Functions  #
# ----------- #

hierarchy = dict()
output_directory = "/home/anamikas/work/FunGNet/data/multilevel_dataset_dictionary/"
train_set = True
fragment_len = 100
train_fraction = 0.8
genomes_per_node_train = 5000
genomes_per_node_test = (1-train_fraction)*genomes_per_node_train

next_taxon_level = {"species": "genus",
                    "genus": "family",
                    "family": "order",
                    "order": "class",
                    "class": "phylum",
                    "phylum": "kingdom",
                    "kingdom": "domain",
                    "domain": "root"}

def count_total_records(filename):
    count = 0
    for record in SeqIO.parse(filename, "fasta"):
        count = count + 1
    return count

def get_desired_ranks(taxid, desired_ranks):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in
                         lineage2ranks.items())

    return [ranks2lineage.get(rank, '0') for rank in desired_ranks]

def get_parent2children(all_fasta_files, parent_level):

    parent2child = dict()
    for fasta_file in os.listdir(all_fasta_files):
        key = get_desired_ranks(int(fasta_file.strip().split('.')[0]),
                                [parent_level])[0]
        if key not in parent2child:
            parent2child[key] = []
        parent2child[key].append(int(fasta_file.strip().split('.')[0]))
    return parent2child

def get_parent2children_next_level(child_dataset_path, parent_level):
    
    parent2child = dict()
    for key, value in child_dataset_path.items():
        rank = get_desired_ranks(key, [parent_level])[0]
        if rank not in parent2child:
            parent2child[rank] = []
        parent2child[rank].append(key)
    return parent2child


def create_from_species(all_fasta_files, current_taxon_level,
                        final_taxon_level):
    
    try:
        next_taxon_level[final_taxon_level]
    except:
        return
    
    if current_taxon_level == next_taxon_level[final_taxon_level]:
        return
    
    directory_for_ml_training = 'ml_training'
    directory_for_db_creation = 'db_creation'
    directory_for_ml_testing = 'ml_test'
   
    print "\n Creating dataset at taxonomy level: {} \n".format(current_taxon_level)

    species2subspecies = get_parent2children(all_fasta_files, current_taxon_level)
    
    new_output_directory = os.path.join(output_directory, current_taxon_level)

    skip_length = 100
    start_pos = 0

    for key, value in tqdm(species2subspecies.items()):
	
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_training,
                    str(key)))
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_testing,
                    str(key)))

        if train_set is True:
            os.makedirs(os.path.join(new_output_directory,
                        directory_for_db_creation, str(key)))

        with open(os.path.join(new_output_directory, directory_for_ml_training,
                  str(key), str(key) + ".fasta"), 'a') as fasta_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_training,
                 str(key), str(key) + ".taxid"), 'a') as taxid_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                 str(key), str(key) + ".fasta"), 'a') as fasta_handle_test,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                 str(key), str(key) + ".taxid"), 'a') as taxid_handle_test:
        
            # if train_set is True:
            #     fasta_db_handle = open(os.path.join(new_output_directory,
            #                            directory_for_db_creation, str(key),
            #                            str(key) + ".fasta"), 'a')
            #     taxid_db_handle = open(os.path.join(new_output_directory,
            #                            directory_for_db_creation, str(key),
            #                            str(key) + ".taxid"), 'a')
        
            for fasta_file in tqdm(value):
                for fasta_record in SeqIO.parse(os.path.join(all_fasta_files, str(fasta_file) + ".fna"),  'fasta'):
                    
                    rec_range = range(len(xrange(start_pos,len(fasta_record.seq) - fragment_len + 1,skip_length)))
                    random.shuffle(rec_range)
                    test_indexes = rec_range[int(train_fraction * len(rec_range)) + 1:]
                                        
                    for cur_position in xrange(start_pos, len(fasta_record.seq) - fragment_len + 1, skip_length):
                        text = str(fasta_record.seq)[cur_position:cur_position + fragment_len]
                        if len(text) != fragment_len:
                            continue
                        
                        if cur_position/skip_length in test_indexes:
                            fasta_handle_test.write(">" + fasta_record.id + "\n")
                            fasta_handle_test.write(text + "\n")
                            taxid_handle_test.write(str(fasta_file) + "\n")

                        else:
                            fasta_handle_train.write(">" + fasta_record.id + "\n")
                            fasta_handle_train.write(text + "\n")
                            taxid_handle_train.write(str(fasta_file) + "\n")

                        
            #         if train_set is True:
            #             fasta_db_handle.write(">" + fasta_record.id + "\n")
            #             fasta_db_handle.write(str(fasta_record.seq) + "\n")
            #             taxid_db_handle.write(str(fasta_file) + "\n")

            # if train_set is True:
            #     fasta_db_handle.close()
            #     taxid_db_handle.close()
   
    hierarchy[current_taxon_level] = species2subspecies

    
    return create_next_dataset(hierarchy[current_taxon_level], next_taxon_level[current_taxon_level],
                               final_taxon_level)

def create_next_dataset(child_dataset, current_taxon_level,
                        final_taxon_level):
    
    try:
        next_taxon_level[final_taxon_level]
    except:
        return

    if current_taxon_level == next_taxon_level[final_taxon_level]:
        return
    
    directory_for_ml_training = 'ml_training'
    directory_for_ml_testing = 'ml_test'

    print "\n \n Creating dataset at taxonomy level: {} \n".format(current_taxon_level)

    parent2child = get_parent2children_next_level(child_dataset,
                                                  current_taxon_level)
    hierarchy[current_taxon_level] = parent2child
    new_output_directory = os.path.join(output_directory, current_taxon_level)
    start_pos = 0   
    skip_length = 100
    
    for key, value in tqdm(hierarchy[current_taxon_level].items()):
        prob = 1.0/len(value)
        genomes_from_child_train = math.ceil(prob*genomes_per_node_train)
        genomes_from_child_test = math.ceil(prob*genomes_per_node_test)
        
        
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_training,
                                 str(key)))
        os.makedirs(os.path.join(new_output_directory, directory_for_ml_testing,
                                 str(key)))

        
        child = [k for k, v in next_taxon_level.items() if v==current_taxon_level]
        child_dataset_path = os.path.join(output_directory, child[0]) 

        
        with open(os.path.join(new_output_directory, directory_for_ml_training,
                  str(key), str(key) + ".fasta"), 'a') as fasta_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_training,
                 str(key), str(key) + ".taxid"), 'a') as taxid_handle_train,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                 str(key), str(key) + ".fasta"), 'a') as fasta_handle_test,\
            open(os.path.join(new_output_directory, directory_for_ml_testing,
                 str(key), str(key) + ".taxid"), 'a') as taxid_handle_test: 
  
              
            for directory in tqdm(value):
                
                current_train_directory = os.path.join(child_dataset_path,
                                            directory_for_ml_training,
                                            str(directory))
            
                for fasta_file in os.listdir(current_train_directory):
                    if fasta_file.endswith('.taxid'):
                        continue
                    
                    total_records = count_total_records(os.path.join(current_train_directory, fasta_file))
                    indices_to_choose = []
                    try:
                        indices_to_choose = random.sample(range(0, total_records), int(genomes_from_child_train))
                    except:
                        indices_to_choose = random.sample(range(0, total_records), int(total_records))
                    count = 0                    
                    for fasta_record in SeqIO.parse(os.path.join(current_train_directory, fasta_file), 'fasta'):
                        count += 1
                        if count in indices_to_choose:
                            fasta_handle_train.write(">" + fasta_record.id + "\n")
                            fasta_handle_train.write(str(fasta_record.seq) + "\n")
                            taxid_handle_train.write(str(directory) + "\n")
                
                
                current_test_directory = os.path.join(child_dataset_path,
                                            directory_for_ml_testing,
                                            str(directory))
            
                for fasta_file in os.listdir(current_test_directory):
                    if fasta_file.endswith('.taxid'):
                        continue
                    
                    total_records = count_total_records(os.path.join(current_test_directory, fasta_file))
                    indices_to_choose = []
                    try:
                        indices_to_choose = random.sample(range(0, total_records), int(genomes_from_child_test))
                    except:
                        indices_to_choose = random.sample(range(0, total_records), int(total_records))
                    count = 0
                    
                    for fasta_record in SeqIO.parse(os.path.join(current_test_directory, fasta_file), 'fasta'):
                        count += 1
                        if count in indices_to_choose:
                            fasta_handle_test.write(">" + fasta_record.id + "\n")
                            fasta_handle_test.write(str(fasta_record.seq) + "\n")
                            taxid_handle_test.write(str(directory) + "\n")

                            
    return create_next_dataset(hierarchy[current_taxon_level], next_taxon_level[current_taxon_level], final_taxon_level)


if __name__=="__main__":
    start_time = time.time()
    create_from_species('/home/anamikas/work/FunGNet/data/all_genomes/', "species", "class")
    end_time = time.time()
    print 'Time taken in buiding dataset', end_time-start_time
