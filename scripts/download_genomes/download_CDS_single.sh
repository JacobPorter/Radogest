#!/bin/sh
nohup srun --partition standard --exclusive -t 6-23 ~/Applications/Radogest/library/ncbi-genome-download/ncbi-genome-download-runner.py -s $1 -F cds-fasta -o /scratch/jsp4cu/Data/GenomesCD/ -p 20 -r 5 -j /scratch/jsp4cu/Data/GenomesCD/json/$2.$1.json $2 &
