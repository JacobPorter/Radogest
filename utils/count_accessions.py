#!/usr/bin/env python
from collections import defaultdict

file_name = "1.fasta"

accession_set = defaultdict(set)

fd = open(file_name)

for line in fd:
    if line.startswith(">"):
        line_list = line.split(":")
        accession_set[int(line_list[1])].add(line_list[0][1:])

for taxid in accession_set:
    print(taxid, len(accession_set[taxid]))
