# Radogest: __ra__n__do__m __ge__nome __s__ampler for __t__rees
Radogest randomly samples fixed length nucleotide substrings (kmers) for given taxonomic ids from a data store of genomes downloaded from the National Center for Biotechnology Information (NCBI).  There is support for whole genomes, coding domain nucleotide data, and amino acid data.  (However, each data type needs to be stored in its own directory.)  Radogest is useful for generating kmers to train and analyze metagenomic classifiers, and it labels each kmer sampled with the taxonomic id that it represents.  Radogest can generate data for any NCBI taxonomic id that is present in the NCBI data store that is of a conventional rank.  Conventional ranks include superkingdom, kingdom, phylum, class, order, genus, and species.


## Install Radogest

Radogest source code can be found on github.  It requires the SeqIterator submodule, which can be downloaded at the same time with recursion.


```
git --recurse-submodules clone ...
```

Radogest is written in Python 3, and it requires the Python packages ete3, tqdm, appdirs, and requests.  Radogest requires installing SAMtools (https://github.com/samtools/samtools) and bedtools (https://bedtools.readthedocs.io/en/latest/).  SAMtools is used to generate fasta files, and bedtools is used to generate the locations of randomly drawn kmers.  The location where SAMtools and bedtools are stored can be documented in config.py.  Radogest tries the locations in the config.py file, and if those are not valid executables, then Radogest searches for the exectuables for SAMtools and bedtools in the operating system PATH variable.

## Important Files

### radogest.py
The main program.  Usage: `radogest.py <command> <options>`.  For help, type `radogest.py -h` or `radogest.py <command> -h`.

### config.py
This file stores the locations of samtools and BED tools.  This file may need to be modified with the locations of the programs to make radogest.py work correctly.

## Commands

### download
Used to download genomes from NCBI.  Modified by Jacob Porter from [https://github.com/kblin/ncbi-genome-download].  In addition to downloading genomes, this downloads genome information in JSON so that a genomes index can be created.  

*All file types need to be downloaded into separate directories because Radogest only recognizes one fasta file in a directory.*  For examples, all NCBI whole genomes could be downloaded into GenomesNT, and all coding domain fasta files could be downloaded into GenomesCD, and all amino acid fasta files could be downloaded into GenomesAA. 

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
Randomly samples DNA substrings for a given taxonomic id.  The taxonomic id must be present in the tree and in the genomes index.  This creates a taxid file labeling each DNA substring drawn with a taxonomic id from the children of the taxonomic id given to the sampler.  This creates a fasta file with each sampled substring.

Given a taxonomic id file and a fasta file, this creates up to three sets of randomly permuted data.  This is useful for training, validating, and testing a classifier.

Radogest can sample from multiple taxonomic ids by giving it a list of taxonomic ids where each id is on a single line in a file.  The sampling can use multiprocessing where a process is used to sample from each  taxonomic id.  When a single taxonomic id is given, Radogest can use multiple processes to generate data for that single taxonomic id.

### util_permute
Randomly permutes and optionally splits a fasta file and a taxid file.

### util_chop
Cuts up a genome into kmers based on a sliding window.  Outputs a fasta file.

### util_rc
Randomly applies the reverse complement to DNA sequences in a fasta file.  Outputs a new fasta file.


## Workflow

The commands must be done in the following order to generate fasta files of sampled kmers.  Once step 5, select, is done, any number of fasta files may be generated in step 6.  Step 5 may need to be done several times if many genome selection strategies are desired.

1. download
2. faidx
3. index
4. tree
5. select
6. sample
 
## Examples

Download both refseq and genbank data into the /project/GenomesNT directory.  Data can be downloaded for each domain, but this example downloads all NCBI genomes at once.  This example uses 20 processes on a 20 core machine to accelerate downloading.  There are 5 attempts to download a file until the process times out.

```
radogest.py download --format fasta -p 20 -r 5 -s refseq,genbank -o /project/GenomesNT /project/Genomes/NT/json/all.json all
```

The next step is to build a samtools fai index for each fasta file.  This step requires samtools.  This can be done in parallel, so we will use 10 processes.

```
radogest.py faidx --genomes /project/GenomesNT -p 10
```


Radogest's genomes index must be built.  This requires the json output from download.  It should only take 10-60 minutes to build.  The index is a pickled Python dictionary that can be opened and queried.

```
radogest.py index --index /project/GenomesNT/index.pck --genomes /project/GenomesNT /project/Genomes/NT/json/
```

A taxonomic tree is built so that Radogest can know which conventional ranks exist.  It requires Radogest's genomes index.  It produces a list of non-species taxonomic ids at conventional ranks.  This is useful if a data set for all taxonomic ids is desired.

```
radogest.py tree --index /project/GenomesNT/index.pck --tree /project/GenomesNT/tree.pck --taxid /project/GenomesNT/taxid_list.txt
```

Strategies exist for Radogest to down select genomes for each taxonomy in the taxonomic tree.  This example uses the method Quality Sorting Leaf (QSL) with at most 10 genomes selected for each species.  This outputs a new Radogest index that can be used in sampling.

```
radogest.py select --index /project/GenomesNT/index.pck --tree /project/GenomesNT/tree.pck --strategy QSL --sample_amount 10 --output /project/GenomesNT/index_QSL.pck
```

Suppose that one million 100-mers representing superkingdom data is desired.  The following command will generate that data given the Radogest index and tree.  It will put the data in a Data directory, and the /tmp/ directory will be used to store temporary files that will be deleted when Radogest is finished.  The taxonomic id 1 represents the superkingdom.  This command creates samples kmers from four kingdoms: virus, archaea, eukaryota, bacteria.  

In this example, four processes will be used to generate the data, and there should be at least four processing cores.  When using multiprocessing on a single taxonomic id, two more processes than specified with the processes parameter will be used.  One process is the main process that coordinates communication, and one process writes the files.  The other processes are worker processes that generate kmer samples.

```
radogest.py sample --index /project/GenomesNT/index_QSL.pck --tree /project/GenomesNT/tree.pck --genomes /project/GenomesNT --kmer_size 100 --number 1000000 -p 4 --data_dir /project/GenomesNT/Data/ --temp_dir /tmp/ 1
```

The kmers can be split into train, validation, and test sets.  The percentage of kmers allocated to these sets can be given with the `split_amount` option.  The following example will produce a training set of 75% of the kmers and a testing set of 25% of the kmers since there are two percentages given for the `split_amount.`  If a validation data set is desired, then a third percentage must be given.

```
radogest.py sample --index /project/GenomesNT/index_QSL.pck --tree /project/GenomesNT/tree.pck --genomes /project/GenomesNT --split --split_amount 0.75,0.25 --data_dir /project/GenomesNT/Data/ --temp_dir /tmp/ 1
```

A file listing taxonomic ids can be used instead of single taxonomic ids.  This option allows for parallel processing where one process is used to sample one taxonomic id.  Each taxonomic id will produce a fasta file for that id.  The range option gives a range of taxonomic ids to sample from similar to a Python slice.  For example, `--range 5 10` will generate samples for the fifth through the tenth taxonomic id in the taxonomic id file.  The taxonomic id file generated by the tree command is useful here.  Use `--amino_acid` to indicate that amino acid data will be sampled from.  Amino acid data must be downloaded using `download`, and the index and tree must be generated. 

```
radogest.py sample --index /project/GenomesAA/index_QSL.pck --tree /project/GenomesAA/tree.pck --genomes /project/GenomesAA --range 5 10 -p 5 --data_dir /project/GenomesAA/Data/ --amino_acid taxid_list.txt
```

Two options, `thresholding` and `chop`, can be used to cut up whole genomes rather than draw kmers uniformly at random.  Thresholding cuts up the genome into kmers if the amount of genomic data requested from a genome is larger than the size of the genome.  Chop cuts the genome into kmers no matter what.  Both options work with the `window_length` option to determine the length of the sliding window from which kmers are produced.  When using the genome holdout strategy, the chop option may be desired to cut up whole genomes for training and testing data.  Genome holdout only produces training and testying data, and the `split` option has no effect.

The option `prob` gives the probability that the reverse complelent will be applied.  This is only valid for nucleotide data as amino acids do not have reverse complements.

The following example includes chopping with the genome holdout strategy where no reverse complements are taken.

```
radogest.py sample --index /project/GenomesNT/index_GHT.pck --tree /project/GenomesNT/tree.pck --genomes /project/GenomesNT --prob 0.0 --data_dir /project/GenomesNT/Data/ 1
```


## Down selection strategies
Down selection uses a post-order depth-first search of the taxonomic tree to chooses genomes for each taxonomic id in the taxonomic tree.  There are several genome selection strategies.  

The genome holdout strategies hold out entire genomes.  The sampling strategy will chop up whole genomes if this strategy is selected.  Genome holdout strategies may not be appropriate for very high level taxonomic ids since the number of kmers sampled will be unbalanced since eukaryotic genomes could be much larger than virus genomes, for example.  For sampling kmers at the genus level, genome holdout makes no sense since this strategy holds out entire species for training or testing.

### AllGenomes (AG)
Uses all non-redundant genomes.  This set could be very large for high level taxonomic ids.

### ProportionalRandom (PR)
At each level, this selection strategy chooses a fixed number of random genomes.  Only genomes selected by children nodes are chosen.

### QualitySortingTree (QST)
At each level, the genomes are sorted by quality, and a fixed number of genomes is chosen.

### QualitySortingLeaf (QSL)
At the species level, the genomes are sorted by quality and a fixed number of genomes is chosen.  These genomes are propagated up the tree so that each parent node selects all of the genomes that each species node under the parent contains.

### GenomeHoldoutLeaf (GHL)
Alternately label each species as part of the training or the testing data and choose a fixed number of genomes sorted by quality.  Propagate this labeling up the taxonomic tree. 

### GenomeHoldoutTree (GHT)
This is the recommended holdout strategy since it will be better balanced at intermediate ranks.  It works like GenomeHoldoutLeaf except that the number of genomes in the testing and the training data at each taxonomic level is fixed.

## Index format

The Radogest genomes index gives Radogest information about genomes and about which genomes are included for which taxonomic ids.  The index is marked with a genome down selection strategy.

While in a Python iterpreter, the following code can be used to load Radogest's genomes index.

```
import pickle
index = pickle.load(open("index.pck", "rb"))
```

There are three parts of the index.

The genomes section, `index['genomes']`, contains a dictionary indexed on NCBI genome accession id.  For each accession, assorted information about the genome or fasta file is included.  For example, `index['genomes']['GCF_000913175.1']` contains information on a virus genome.

The taxid section, `index['taxids']`, gives the genome accession ids included for each taxonomic id.  For example, `index['taxids'][9604]` gives a dictionary keyed on genome accessions with boolean or integer values that indicate whether that genome accession should be included for that taxonomic id while sampling.  The id `9604` is for the family Hominidae.

The selection strategy section, `index['select']`, gives the genome selection strategy used during the `select` command to create the index.  Typing `index['select']` gives a genome strategy id and the sample amount used during the select command.

## Tree format

The Radogest taxonomic tree is a Python dictionary that gives the parent and child relationship between taxonomic ids.  The key is a taxonomic id, and the value is a set of taxonomic ids that are children of the key.  Only conventional ranks are included in this tree.  If there exists an unbranching simple path between nodes, it is removed.

While in a Python interpreter, the tree can be loaded with the following code.

```
import pickle
tree = pickle.load(open("tree.pck", "rb"))
```

Typing `tree[1]`, for example, could give the set `{2, 2157, 10239, 2579}`.  The taxonomic id `1` may be a root node, and the taxonomic ids `2`, `2157`, `10239`, and `2579` are the ids for bacteria, archaea, viruses, and eukaryotes.

## Fasta record format

Radogest sample generates fasta records with a fasta id in the following format:

`>[Accession]:[TaxonomicID]:[ContigID]:[StartPosition-EndPosition]:[+|-]`

The last character is a `+` if the kmer is the forward strand, and it is a `-` for the reverse complement.  The start and end positions are relative to the forward strand.  This means that the start position is always numerically less than the end position.