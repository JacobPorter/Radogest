#!/bin/sh
# section (genbank, refseq)
./faidx_NT_single.sh $1 archaea &
./faidx_NT_single.sh $1 bacteria &
./faidx_NT_single.sh $1 fungi &
./faidx_NT_single.sh $1 invertebrate &
./faidx_NT_single.sh $1 plant &
./faidx_NT_single.sh $1 protozoa &
./faidx_NT_single.sh $1 vertebrate_mammalian &
./faidx_NT_single.sh $1 vertebrate_other &
./faidx_NT_single.sh $1 viral &
