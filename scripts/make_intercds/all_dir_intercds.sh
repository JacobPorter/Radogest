#!/bin/sh
./single_dir_intercds.sh refseq archaea &
./single_dir_intercds.sh refseq bacteria &
./single_dir_intercds.sh refseq fungi &
./single_dir_intercds.sh refseq invertebrate &
./single_dir_intercds.sh refseq plant &
./single_dir_intercds.sh refseq protozoa &
./single_dir_intercds.sh refseq vertebrate_mammalian &
./single_dir_intercds.sh refseq vertebrate_other &
./single_dir_intercds.sh refseq viral &

./single_dir_intercds.sh genbank archaea &
./single_dir_intercds.sh genbank bacteria &
./single_dir_intercds.sh genbank fungi &
./single_dir_intercds.sh genbank invertebrate &
./single_dir_intercds.sh genbank plant &
./single_dir_intercds.sh genbank protozoa &
./single_dir_intercds.sh genbank vertebrate_mammalian &
./single_dir_intercds.sh genbank vertebrate_other &
./single_dir_intercds.sh genbank viral &
