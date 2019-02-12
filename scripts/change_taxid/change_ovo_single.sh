#!/bin/sh
# $1 type 1
# $2 type 2
# $3 data set

DATA_PATH="/project/biocomplexity/fungcat/genomes/Genomes_CD/Data/1_ova/$3/"
OUTPUT_PATH="/project/biocomplexity/fungcat/genomes/Genomes_CD/Data/1_ovo/$3/"
mkdir $OUTPUT_PATH/$1_$2/
/home/jsp4cu/Applications/Root_Classifier/utils/transform_vector_inputs.py extract $DATA_PATH/1.fasta $DATA_PATH/$1/1.taxid $OUTPUT_PATH/$1_$2/$1.fasta $OUTPUT_PATH/$1_$2/$1.taxid
/home/jsp4cu/Applications/Root_Classifier/utils/transform_vector_inputs.py extract $DATA_PATH/1.fasta $DATA_PATH/$2/1.taxid $OUTPUT_PATH/$1_$2/$2.fasta $OUTPUT_PATH/$1_$2/$2.taxid
cat $OUTPUT_PATH/$1_$2/$1.fasta $OUTPUT_PATH/$1_$2/$2.fasta &> $OUTPUT_PATH/$1_$2/1.fasta
cat $OUTPUT_PATH/$1_$2/$1.taxid $OUTPUT_PATH/$1_$2/$2.taxid &> $OUTPUT_PATH/$1_$2/1.taxid
/home/jsp4cu/Applications/Radogest/library/old/permute_split_fasta_taxid.py $OUTPUT_PATH/$1_$2/1.fasta $OUTPUT_PATH/$1_$2/1.taxid
rm $OUTPUT_PATH/$1_$2/$1.fasta
rm $OUTPUT_PATH/$1_$2/$2.fasta
rm $OUTPUT_PATH/$1_$2/$1.taxid
rm $OUTPUT_PATH/$1_$2/$2.taxid


