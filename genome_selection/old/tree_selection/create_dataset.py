import argparse
import pickle
import os
from joblib import Parallel, delayed
from multiprocessing import Pool
from Bio import SeqIO
import time
from tqdm import tqdm
import random

taxid2genomes = pickle.load(open('taxid2genomes_random10_tuple.pickle', 'rb'))

directory_for_ml_training = 'ml_training'
directory_for_db_creation = 'db_creation'
directory_for_ml_testing = 'ml_test'
directory = [directory_for_ml_training, directory_for_ml_testing, directory_for_db_creation]

def create_directory(level, output_dir, k, train_fraction, train_set):
    for items in taxid2genomes[level].keys():
        create_dataset(items, level, output_dir, k, train_fraction, train_set)
    #Parallel(n_jobs=len(taxid2genomes[level].keys()))(delayed(create_dataset)(items, level, output_dir,
    #                                        k, train_fraction, train_set) for items in taxid2genomes[level].keys())
    for items in taxid2genomes[level].keys():
        write_files(items, level, output_dir, k, train_fraction, train_set)
    #Parallel(n_jobs=len(taxid2genomes[level].keys()))(delayed(write_files)(items, level, output_dir,
    #                                        k, train_fraction, train_set) for items in taxid2genomes[level].keys())

def create_dataset(keys, level, output_dir, k, train_fraction, train_set):
    for items in directory:
        os.makedirs(os.path.join(output_dir, level, str(items), str(keys)))
    #Parallel(n_jobs=len(taxid2genomes[level].keys()))(delayed(os.makedirs)(os.path.join(output_dir, level, str(items), str(keys))) for items in directory)


def write_files(keys, level, output_dir, k, train_fraction, train_set):
    with open(os.path.join(output_dir, level, directory_for_ml_training,
            str(keys), str(keys) + ".fasta"), 'a') as fasta_handle_train,\
        open(os.path.join(output_dir, level, directory_for_ml_training,
            str(keys), str(keys) + ".taxid"), 'a') as taxid_handle_train,\
        open(os.path.join(output_dir, level, directory_for_ml_testing,
            str(keys), str(keys) + ".fasta"), 'a') as fasta_handle_test,\
        open(os.path.join(output_dir, level, directory_for_ml_testing,
            str(keys), str(keys) + ".taxid"), 'a') as taxid_handle_test:


        if train_set == 'True':
            with open(os.path.join(output_dir, level, directory_for_db_creation,
                str(keys), str(keys) + ".fasta"), 'a') as fasta_db_handle,\
                open(os.path.join(output_dir, level, directory_for_db_creation,
                str(keys), str(keys) + ".taxid"), 'a') as taxid_db_handle:

                for gens, locs in tqdm(taxid2genomes[level][keys].items()):
                    for fasta_record in SeqIO.parse(os.path.join(locs[1], gens), 'fasta'):
                        lst_kmers=[]
                        test_list = []
                        train_list = []
                        kmer_lst = list(count_kmers(fasta_record.seq, int(k)))
                        test_indices = random.sample(range(0, len(kmer_lst)),
                                                    len(kmer_lst)-int(train_fraction*len(kmer_lst)))
                        test_list = [kmer_lst[i] for i in test_indices]
                        train_list = list(set(kmer_lst) - set(test_list))

                        # for train, test in zip(train_list, test_list):
                        #         fasta_handle_train.write(">" + fasta_record.id + "\n")
                        #         fasta_handle_train.write(str(train) + "\n")
                        #         taxid_handle_train.write(str(locs[0]) + "\n")
                        #         fasta_handle_test.write(">" + fasta_record.id + "\n")
                        #         fasta_handle_test.write(str(test) + "\n")
                        #         taxid_handle_test.write(str(locs[0]) + "\n")

                        for kmer in train_list:
                            fasta_handle_train.write(">" + fasta_record.id + "\n")
                            fasta_handle_train.write(str(kmer) + "\n")
                            taxid_handle_train.write(str(locs[0]) + "\n")

                        for kmer in test_list:
                            fasta_handle_test.write(">" + fasta_record.id + "\n")
                            fasta_handle_test.write(str(kmer) + "\n")
                            taxid_handle_test.write(str(locs[0]) + "\n")


                        if train_set == "True":
                            fasta_db_handle.write(">" + fasta_record.id + "\n")
                            fasta_db_handle.write(str(fasta_record.seq) + "\n")
                            taxid_db_handle.write(str(locs[0]) + "\n")

        if train_set == "True":
            fasta_db_handle.close()
            taxid_db_handle.close()

def count_kmers(seq, k):
    set_kmers = set([])
    num_kmers = len(seq)-k+1
    for i in range(num_kmers):
        kmer = seq[i:i+k]
        set_kmers.add(kmer)
    return set_kmers

def main(map_args):
    start_time = time.time()
    output_dir = map_args.output_dir
    k = map_args.fragment_len
    train_fraction = map_args.train_fraction
    train_set = map_args.train_set
    levels = []
    for keys, _ in taxid2genomes.items():
        levels.append(keys)
    for level in levels:
        create_directory(level, output_dir, k, train_fraction, train_set)
    #Parallel(n_jobs=len(levels))(delayed(create_directory)(level, output_dir,
    #                                    k, train_fraction, train_set) for level in levels)
    end_time = time.time()
    print(end_time-start_time)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', help='Path to hold\
			genomes')
    parser.add_argument('--fragment_len', type=int,
                        help='Length of kmers')
    parser.add_argument('--train_fraction', type=float,
                        help='Train and test split')
    parser.add_argument('--train_set', help='True if whole sequences to be kept')
    map_args = parser.parse_args()
    main(map_args)
