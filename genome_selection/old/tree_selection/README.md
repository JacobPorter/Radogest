### Sampling for Plinko
Author: Anamika Sen <anamikas@vt.edu>

This folder contains scripts used for sampling genomes in Plinko. The approach used in sampling can be read about in here:
<https://docs.google.com/document/d/1tCVdyc4faQwFC7jhHiqdhO6i1pqQa2oQvy-rIGJCXaw/edit?usp=sharing>

#### The scripts used are:

    1) create_tree.py
    - Given a root taxon id and the final taxon level, it will create a taxonomy tree of all the children in the tree.
    - The parameters are:
        --taxon_id: Taxon ID of the root node of the tree
        --root_taxon_level: Taxon level of the root node
        --final_taxon_level: Taxon level of the leaves
        --rank_dir: Path of the ranks pickle datastructure
    - Output is a pickled file. Example: Refer the file taxonomy_tree.pickle
##### Example:
 ```  
 python -W ignore create_tree.py --taxon_id 91061 --root_taxon_level class --final_taxon_level species --rank_dir '/groups/fungcat/dataset/current/fasta/Genomes'
```
    2) multilevel_sampleGenomes.py
    - Given a taxonomy tree and nummber of genomes to be sampled, N genomes are randomnly sampled from leaves until the root of      the tree.
    - The parameters are:
        --taxon_id: Taxon ID of the root node of the tree
        --root_taxon_level: Taxon level of the root node
        --final_taxon_level: Taxon level of the leaves
        --data_dir: Path to Genomes directory
        --tree Taxonomy tree pickled file
        --genomes Number of genomes to be sampled at each node
        --subselect_species list of species TaxIDs to be considered instead of all the species        
    - Output is a pickled file. Example: Refer the file taxid2genomes_random10_tuple.pickle
    Output is of the form:
    file:
        {
            'species' :
                {   
                    '1304': {
                            'Genome1': (TaxID, 'path/to/genome/Genome1',
                            'Genome2': (TaxID, 'path/to/genome/Genome2',
                            'Genome3': (TaxID, 'path/to/genome/Genome3',
                            .
                            .
                            
                            }
                    '1308': {
                            'Genome1': (TaxID, 'path/to/genome/Genome1',
                            'Genome2': (TaxID, 'path/to/genome/Genome2',
                            'Genome3': (TaxID, 'path/to/genome/Genome3',
                            .
                            .
                            }
                       .
                       .
                 }
             'genus' :
                {   
                    '1279': {
                            'Genome1': (TaxID of parent, 'path/to/genome/Genome1',
                            'Genome2': (TaxID of parent, 'path/to/genome/Genome2',
                            'Genome3': (TaxID of parent, 'path/to/genome/Genome3',
                            .
                            .
                            }
                    '1386': {
                            'Genome1': (TaxID of parent, 'path/to/genome/Genome1',
                            'Genome2': (TaxID of parent, 'path/to/genome/Genome2',
                            'Genome3': (TaxID of parent, 'path/to/genome/Genome3',
                            .
                            .
                            }
                       .
                       .
                 }
                 .
                 .
               'root_taxon_level' :
                {
                    '91061': {
                             'Genome1': (TaxID of parent, 'path/to/genome/Genome1',
                             'Genome2': (TaxID of parent, 'path/to/genome/Genome2',
                             .
                             .
                             }                                  
                }               
         }
(Genomes from species are randomly selected as N in the current implementation. Come up with a better implementation.)
##### Example:
 ```  
 python -W ignore multilevel_sampleGenomes.py --taxon_id 91061 --root_taxon_level class --final_taxon_level species --data_dir '/groups/fungcat/datasets/current/fasta/Genomes/' --tree taxonomy_tree.pickle --genomes 10 --subselect_species 1304 1308 1311 1314 1328 1336 1351 1280 1390 1396 1402 1406 1428 1580 1582 1587 1598 1613 28035 315405 47715 76860 1282 1307 1309 1313 1318 1334 1338 1352 1392 1398 1404 1423 1579 1581 1584 1590 1604 1624 283734 33959 61624
```
    3) create_dataset.py
    - Create a dataset from sampled genomes pickled file.
    - The parameters are:
        --output_dir: Path to hold genomes
        --fragment_len: Length of kmers
        --train_fraction: Train and test split
        --train_set: True if whole sequences are to be kept in dataset
     - Output is the dataset for training and testing. 
##### Example:
 ```  
 python -W ignore create_dataset.py --output_dir /localscratch/anamikas/output/ --fragment_len 100 --train_frac 0.8 --train_set True
```
