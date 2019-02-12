#!/bin/sh
# $1 file format
# $2  directory
nohup ./download_genomes_single.sh refseq archaea $1 $2 3-1  &
nohup ./download_genomes_single.sh genbank archaea $1 $2 3-1 &
nohup ./download_genomes_single.sh refseq bacteria $1 $2 6-23 &
nohup ./download_genomes_single.sh genbank bacteria $1 $2 6-23 &
nohup ./download_genomes_single.sh refseq fungi $1 $2 3-1 &
nohup ./download_genomes_single.sh genbank fungi $1 $2 3-1 &
nohup ./download_genomes_single.sh refseq invertebrate $1 $2 2-1 &
nohup ./download_genomes_single.sh genbank invertebrate $1 $2 2-1 &
nohup ./download_genomes_single.sh refseq plant $1 $2 3-1 &
nohup ./download_genomes_single.sh genbank plant $1 $2 3-1 &
nohup ./download_genomes_single.sh refseq protozoa $1 $2 3-1 & 
nohup ./download_genomes_single.sh genbank protozoa $1 $2 3-1 &
nohup ./download_genomes_single.sh refseq vertebrate_mammalian $1 $2 1-1 &
nohup ./download_genomes_single.sh genbank vertebrate_mammalian $1 $2 2-1 &
nohup ./download_genomes_single.sh refseq vertebrate_other $1 $2 1-1 &
nohup ./download_genomes_single.sh genbank vertebrate_other $1 $2 2-1 &
nohup ./download_genomes_single.sh refseq viral $1 $2 3-1 &
nohup ./download_genomes_single.sh genbank viral $1 $2 3-1 &
