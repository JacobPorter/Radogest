#!/bin/sh
./change_ovo_single.sh archaea bacteria test
./change_ovo_single.sh archaea bacteria validate
./change_ovo_single.sh archaea bacteria train

./change_ovo_single.sh archaea eukaryotes test
./change_ovo_single.sh archaea eukaryotes validate
./change_ovo_single.sh archaea eukaryotes train

./change_ovo_single.sh archaea virus test
./change_ovo_single.sh archaea virus validate
./change_ovo_single.sh archaea virus train

./change_ovo_single.sh bacteria eukaryotes test
./change_ovo_single.sh bacteria eukaryotes validate
./change_ovo_single.sh bacteria eukaryotes train

./change_ovo_single.sh bacteria virus test
./change_ovo_single.sh bacteria virus validate
./change_ovo_single.sh bacteria virus train

./change_ovo_single.sh virus eukaryotes test
./change_ovo_single.sh virus eukaryotes validate
./change_ovo_single.sh virus eukaryotes train
