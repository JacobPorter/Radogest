# Random Sampler
This is a random sampler for selecting fixed length nucleotide substrings from a taxonomic tree of genomes.  It requires Python 3.

This requires samtools, UCSC tools, and BED tools.

## Components

### NCBI-Genome-Download
Used to download genomes from NCBI.  Modified by Jacob Porter from [https://github.com/kblin/ncbi-genome-download].  In addition to downloading genomes, this downloads genome information in JSON so that a genomes index can be created.

### get_paths.py
This file stores the locations of samtools, UCSC tools, and BED tools.  It stores the location of where the genomes are stored.

### make_genome_index.py
Makes a python dictionary representing genome locations and genome information.  This is called the genomes index.  Creates an index that associates taxonomic ids to genomes.

### run_faidx_genomes.py
Runs FAIDX from samtools on the genomes in the genomes directory and updates the genomes index.  The FAI file is used to perform random sampling. The genomes index is updated with the number of nucleotide bases and the number of contigs in the genome.

### make_tree.py
Makes a python dictionary that represents the taxonomic tree.  This is called the tree. Requires the index from make_genome_index.py.  This presently uses only conventional taxonomic ranks: superkingdom, kingdom, phylum, class, order, genus, species.  It excludes environmental samples, unclassified samples, and RNA viruses.  Viroids are present but not at the root node 1.

### genome_selection/strategy.py
Implements strategies to down select genomes at each taxonomic level.  This includes ProportionalRandom, which samples genomes proportionaly at random and QualitySortTree, which sorts genomes by quality (first by reference type, then by assembly level, and then by number of bases divided by contigs) and chooses the best ones.

### genome_selection/traversal.py
Performs a post-order depth first search of the taxonomic tree to down select genomes at each taxonomic level.  Requires the tree and the genomes index data structures.  This uses the strategies found in genome_selection/strategy.py.

### random_genome_sampler.py
Randomly samples DNA substrings for a given taxonomic id.  The taxonomic id must be present in the tree and in the genomes index.
This creates a taxid file labeling each DNA substring drawn with a taxonomic id from the children of the taxonomic id given to the sampler.
This creates a fasta file with each sampled substring.

### permute_split_fasta_taxid.py
Given a taxonomic id file and a fasta file, this creates up to three sets of randomly permuted data.  This is useful for training, validating, and testing a classifier.

### parallel_sampler.py
Given a list of taxonomic ids, this can randomly sample DNA substrings for each taxonomic id in parallel using multiprocessing.

### SeqIterator.py and Constants.py
A utility file for iterating through fasta files.  Used by other components.

## Workflow

1. Download genomes with NCBI-Genome-Download with saved JSON files.
2. Update get_paths.py with the locations of the genomes folder, and the locations of samtools, BED tools, and UCSC tools.  This is done by editing the file.
3. Make the genomes index with make_genome_index.py
4. Create fasta index files (fai files) for all of the genomes.  This updates the genomes index from the previous step.
5. Make the taxonomic tree with make_tree.py
6. Genome down selection can be done with genome/selection/traversal.py.  There are several strategies.  This will modify the genomes index data structure, and it requires the tree data structure.
7. random_genome_sampler.py and permute_split_fasta_taxid.py can be used to sample from a single taxonomic id and split the data.  parallel_sampler.py can sample from multiple taxonomic ids in the tree in parallel and will split the data.
 
## Example

TODO
