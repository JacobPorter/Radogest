#!/bin/sh
# $1 file format
# $2  directory
nohup ./download_genomes_single.sh refseq bacteria $1 $2 6-23:55:55 &
nohup ./download_genomes_single.sh genbank bacteria $1 $2 6-23:55:55 &
