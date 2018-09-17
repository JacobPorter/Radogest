#!/usr/bin/env python
"""
Sort and downselect genomes.
"""

import pickle
from ete3 import NCBITaxa
from pprint import pprint

index = pickle.load(open('genomes.withcounts.index.pck', 'rb'))
ncbi = NCBITaxa()
NUMBER_OF_GENOMES_TO_KEEP = 5

accessions_to_keep = {}
assembly_level = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]

def append_by_assembly(my_genomes, genomes_list):
    for genome in genomes_list:
        for assembly in assembly_level:
            if genome['assembly_level'] == assembly:
                my_genomes.append(genome['assembly_accession'])

def contig_key(genome):
    return float(genome['base_length'] / genome['contig_count'])

def add_genomes_to_include(my_genomes):
    for accession in my_genomes:
        accessions_to_keep[accession] = True

#printMe = True
for taxid in index['taxids']:
    try:
        rank_type = ncbi.get_rank([taxid])[taxid]
    except KeyError:
        continue
    if rank_type == 'species':
        my_genomes = []
        genomes = index['taxids'][taxid]
        if len(genomes) <= NUMBER_OF_GENOMES_TO_KEEP:
            add_genomes_to_include(genomes.keys())
            continue
        ref_genomes = [index['genomes'][accession] for
                       accession in genomes if
                       index['genomes'][accession]['refseq_category'] ==
                       'reference genome']
        repr_genomes = [index['genomes'][accession] for
                        accession in genomes if
                        index['genomes'][accession]['refseq_category'] ==
                        'representative genome']
        other_genomes = [index['genomes'][accession] for
                         accession in genomes if
                         index['genomes'][accession]['refseq_category'] !=
                         'reference genome' and
                         index['genomes'][accession]['refseq_category'] !=
                         'representative genome']
        append_by_assembly(my_genomes, repr_genomes)
        append_by_assembly(my_genomes, ref_genomes)
        if len(my_genomes) > NUMBER_OF_GENOMES_TO_KEEP:
            add_genomes_to_include(my_genomes[0:NUMBER_OF_GENOMES_TO_KEEP])
            continue
        other_genomes.sort(key=contig_key, reverse=True)
        # if printMe:
        #     printMe = False
        #     for genome in other_genomes:
        #         print(contig_key(genome))
        #     pprint(other_genomes)
        number_to_add = NUMBER_OF_GENOMES_TO_KEEP - len(my_genomes)
        for genome in other_genomes:
            if number_to_add <= 0:
                break
            my_genomes.append(genome['assembly_accession'])
            number_to_add -= 1
        add_genomes_to_include(my_genomes)

for taxid in index['taxids']:
    genomes_dictionary = index['taxids'][taxid]
    for accession in genomes_dictionary:
        if accession not in accessions_to_keep:
            genomes_dictionary[accession] = False

pickle.dump(index, open('genomes.withcounts.selected.index.pck', 'wb'))
pprint(index, open('genomes.withcounts.selected.index.pck.txt','w'))
