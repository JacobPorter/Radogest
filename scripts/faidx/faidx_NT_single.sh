#!/bin/sh
# $1 section
# $2 kingdom
nohup srun -A fungcat -t 6-23 --exclusive --partition bii ~/Applications/Radogest/radogest.py faidx -p 10  --genomes  /project/biocomplexity/fungcat/genomes/Genomes_NT/$1/$2/ &> /project/biocomplexity/fungcat/genomes/Genomes_NT/out/faidx.$1.$2.out &
