#!/bin/sh
# Extract all labels into files.
# $1 records file
# $2 labels file
# $3 output file
./intercds_fasta.py extract --label 0  $1 $2 1> $3.0.records
./intercds_fasta.py extract --label 1  $1 $2 1> $3.1.records
./intercds_fasta.py extract --label 2  $1 $2 1> $3.2.records
./intercds_fasta.py extract --label 3  $1 $2 1> $3.3.records
