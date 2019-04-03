#!/bin/sh
# $1 genbank or refseq
# $2 section
nohup srun -A fungcat --partition bii -t 6-23 --exclusive ~/Applications/Radogest/utils/intercds_fasta.py all_cds --processes 20 --intercds /project/biocomplexity/fungcat/genomes/Genomes_iCD/$1/$2/ /project/biocomplexity/fungcat/genomes/Genomes_CD/$1/$2/ /project/biocomplexity/fungcat/genomes/Genomes_NT/$1/$2/ &> /project/biocomplexity/fungcat/genomes/Genomes_iCD/out/$1.$2.intercds.out &
