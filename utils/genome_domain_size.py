#!/usr/bin/env python

import pickle
from collections import defaultdict

index_loc = "/project/biocomplexity/fungcat/genomes/Genomes_NT/index_initial.pck"

index = pickle.load(open(index_loc, "rb"))
sizes = defaultdict(int)
counts = defaultdict(int)
maxes = defaultdict(int)
mins = defaultdict(int)
mins_acc = {}
maxs_acc = {}
for accession in index["genomes"]:
    domain = index["genomes"][accession]['location'].split("/")[2]
    try:
        contig_sum = float(index["genomes"][accession]['contig_sum'])
    except KeyError:
        continue
    counts[domain] += 1
    sizes[domain] += contig_sum
    if contig_sum > maxes[domain]:
        maxes[domain] = contig_sum
        maxs_acc[domain] = accession
    if contig_sum < mins[domain] or mins[domain] == 0:
        mins[domain] = contig_sum
        mins_acc[domain] = accession

avg_list = [(domain, sizes[domain]/counts[domain]) for domain in counts]
avg_list.sort(key=lambda x: x[1])
# print(maxs_acc)
# print(mins_acc)
div_amt = 1000000
print("Numbers are in megabases.")
for domain, avg in avg_list:
    output_str = "Domain: {}, avg: {}, max: {} ({}), min: {} ({})".format(domain, avg/div_amt, maxes[domain]/div_amt, maxs_acc[domain], mins[domain]/div_amt, mins_acc[domain])
    print(output_str)
