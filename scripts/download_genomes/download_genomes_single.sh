#!/bin/sh
# $1 section (refseq, genbank)
# $2 taxonomic group
# $3 file format
# $4 directory
# $5 time
nohup srun -A fungcat --partition bii --exclusive -t $5 ~/Applications/Radogest/radogest.py download -s $1 -F $3 -o /project/biocomplexity/fungcat/genomes/$4/ -p 40 -r 10 /project/biocomplexity/fungcat/genomes/$4/json/$2.$1.json $2 &> /project/biocomplexity/fungcat/genomes/$4/out/$2.$1.$3.out &
