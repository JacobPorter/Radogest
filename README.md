# Radogest: random genome sampler for trees
This is a random sampler for selecting fixed length nucleotide substrings from a taxonomic tree of genomes.  

It requires Python 3 and ete3.  This requires downloading samtools and BED tools.

## Install Radogest

```
git clone ....
git submodule init
git submodule update
```

## Important Files

### radogest.py
The main program.  Usage: `radogest.py <command> <options>`.  For help, type `radogest.py -h` or `radogest.py <command> -h`.

### config.py
This file stores the locations of samtools and BED tools.  This file may need to be modified with the locations of the programs to make radogest.py work correctly.

## Commands

### download
Used to download genomes from NCBI.  Modified by Jacob Porter from [https://github.com/kblin/ncbi-genome-download].  In addition to downloading genomes, this downloads genome information in JSON so that a genomes index can be created.

### faidx
Runs FAIDX from samtools on the genomes in the genomes directory and updates the genomes index.  The FAI file is used to perform random sampling. The genomes index is updated with the number of nucleotide bases and the number of contigs in the genome.

### index
Makes a python dictionary representing genome locations and genome information.  This is called the genomes index.  Creates an index that associates taxonomic ids to genomes.

### tree
Makes a python dictionary that represents the taxonomic tree.  This is called the tree. Requires the index from make_genome_index.py.  This presently uses only conventional taxonomic ranks: superkingdom, kingdom, phylum, class, order, genus, species.  It excludes environmental samples, unclassified samples, and RNA viruses.  Viroids are present but not at the root node 1.

### select
Implements strategies to down select genomes at each taxonomic level.  This includes ProportionalRandom, which samples genomes proportionaly at random and QualitySortTree, which sorts genomes by quality (first by reference type, then by assembly level, and then by number of bases divided by contigs) and chooses the best ones.

Performs a post-order depth first search of the taxonomic tree to down select genomes at each taxonomic level.  Requires the tree and the genomes index data structures.  This uses the strategies found in genome_selection/strategy.py.

### sample
Randomly samples DNA substrings for a given taxonomic id.  The taxonomic id must be present in the tree and in the genomes index.
This creates a taxid file labeling each DNA substring drawn with a taxonomic id from the children of the taxonomic id given to the sampler.
This creates a fasta file with each sampled substring.

Given a taxonomic id file and a fasta file, this creates up to three sets of randomly permuted data.  This is useful for training, validating, and testing a classifier.

Given a list of taxonomic ids, this can randomly sample DNA substrings for each taxonomic id in parallel using multiprocessing.

## Other Files

### SeqIterator.py and Constants.py
A utility file for iterating through fasta files.  Used by other components.

## Workflow

1. download
2. faidx
3. index
4. tree
5. select
6. sample
 
## Examples

TODO
