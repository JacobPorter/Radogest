"""
Genome selection strategies.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import sys
from collections import defaultdict
from itertools import accumulate, chain
from random import random, shuffle

EXCLUDED_GENOMES = {}

# The id to use with Genome holdout genome strategies
# when a taxid has a single genome.
SINGLETON = -1


def filter_genomes(accessions, index):
    """
    Filter unwanted genomes.

    Parameters
    ----------
    accessions: iterable
        An iterable of genome accession strings.
    index: dict
        The genomes index.

    Returns
    -------
    include, exclude
        Two lists of accessions.  One to include and one to exclude.

    """
    def genbank_duplicate(accession_info, index):
        return (accession_info['section'] == 'genbank'
                and accession_info['gbrs_paired_asm'] != ''
                and accession_info['paired_asm_comp'] == 'identical'
                and index['genomes'][accession_info['gbrs_paired_asm']]
                and (index['genomes'][accession_info['gbrs_paired_asm']]
                     ['species_taxid'] == accession_info['species_taxid']))

    include = []
    exclude = []
    for accession in accessions:
        if 'contig_sum' not in index['genomes'][accession]:
            exclude.append(accession)
            EXCLUDED_GENOMES[accession] = 'the contig_sum does not exist'
        elif len(accessions) == 1:
            include.append(accession)
        elif genbank_duplicate(index['genomes'][accession], index):
            exclude.append(accession)
            EXCLUDED_GENOMES[accession] = 'duplicate genome found in refseq'
        else:
            include.append(accession)
    return include, exclude


def genome_sort(genomes, index):
    """
    Sort the genomes based on the quality of the assembly, etc.
    First, representative genomes are chosen.  This list is sorted by
    assembly level and then in reverse order by the number of bases
    divided by the number of contigs.  Second, the same is done for
    reference genomes.  Third, the same is done for all other genomes.

    Parameters
    ----------
    genomes: list
        A list of genome accession ids.
    index: dict
        The genomes index.  This will be used to lookup information about each
        genome to do sorting.

    Returns
    -------
    A sorted list of genome accession numbers.

    """
    assembly_level = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]

    def contig_key(genome):
        try:
            return float(genome['contig_sum'] / genome['contig_count'])
        except KeyError:
            print(genome, file=sys.stderr)
            raise KeyError

    def sort_by_assembly_and_contig(genomes_data_list):
        assembly_level_dict = {level: [] for level in assembly_level}
        for genome_data in genomes_data_list:
            assembly_level_dict[genome_data['assembly_level']].append(
                genome_data)
        for level in assembly_level_dict:
            assembly_level_dict[level].sort(key=contig_key, reverse=True)
        sorted_list = []
        for level in assembly_level:
            sorted_list += assembly_level_dict[level]
        return sorted_list

    ref_genomes = [
        index['genomes'][accession] for accession in genomes
        if index['genomes'][accession]['refseq_category'] == 'reference genome'
    ]
    repr_genomes = [
        index['genomes'][accession] for accession in genomes
        if index['genomes'][accession]['refseq_category'] ==
        'representative genome'
    ]
    other_genomes = [
        index['genomes'][accession] for accession in genomes
        if index['genomes'][accession]['refseq_category'] != 'reference genome'
        and index['genomes'][accession]['refseq_category'] !=
        'representative genome'
    ]
    my_genomes = sort_by_assembly_and_contig(ref_genomes)
    my_genomes += sort_by_assembly_and_contig(repr_genomes)
    my_genomes += sort_by_assembly_and_contig(other_genomes)
    return [genome_data['assembly_accession'] for genome_data in my_genomes]


def select_equal(list_lists, select_number):
    """
    Select equal elements from a list of lists.

    Parameters
    ----------
    list_lists: list
        A list of lists of elements.  The lists could be unequal in size.
    select_number: int
        An integer of samples to take.

    Examples
    --------
    list_lists is of form [[a, b, c] [1, 2, 3] [x, y].  If select_number is 5,
    then the function will return [a, 1, x, b, 2].

    Returns
    -------
    list
        A list of sampled elements

    """
    pos = [0] * len(list_lists)
    empty = [0] * len(list_lists)
    output_list = []
    iterations = 0
    while len(output_list) < select_number and sum(empty) < len(list_lists):
        i = iterations % len(list_lists)
        j = pos[i]
        if j < len(list_lists[i]):
            output_list.append(list_lists[i][j])
            pos[i] += 1
        else:
            empty[i] = 1
        iterations = i + 1
    return output_list


class GenomeSelection:
    """The base class for genome selection strategies."""
    def __init__(self, index):
        """
        Initialize the GenomeSelection class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.

        """
        self.index = index

    def set_all_genomes(self, boolean=False):
        """
        Set all genome inclusions to the given value
        in the index.

        Parameters
        ----------
        boolean: primitive
            A value to set all genome inclusions to.

        Returns
        -------
        None

        """
        taxids_dict = self.index['taxids']
        for taxid in taxids_dict:
            for accession in taxids_dict[taxid]:
                taxids_dict[taxid][accession] = boolean


class ProportionalRandom(GenomeSelection):
    """
    For a given level in the taxonomy tree, select a random fixed number
    of genomes at each level.
    """
    def __init__(self, index, select_number):
        """
        Initialize the GenomeSelection class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_number: int
            The number of genomes to select at each level.

        """
        super().__init__(index)
        self.select_number = select_number
        self.set_all_genomes(boolean=False)

    def select(self, parent, children):
        """
        Choose genomes from the children randomly to include in the parent
        taxonomic id.  Genomes are selected so that the children genomes are
        uniformly represented.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        if children:
            children_genomes = []
            for child in children:
                child_genomes = []
                child_dict = self.index['taxids'][child]
                include, _ = filter_genomes(child_dict.keys(), self.index)
                for accession in include:
                    if child_dict[accession]:
                        child_genomes.append(accession)
                shuffle(child_genomes)
                children_genomes.append(child_genomes)
            my_genomes = select_equal(children_genomes, self.select_number)
            samples = len(my_genomes)
            for accession in my_genomes:
                self.index['taxids'][parent][accession] = True
            return samples
        else:
            my_genomes, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                           self.index)
            shuffle(my_genomes)
            i = 0
            for accession in my_genomes:
                if i >= self.select_number:
                    break
                self.index['taxids'][parent][accession] = True
                i += 1
            return i


class QualitySortingTree(GenomeSelection):
    """
    For a given level in the taxonomy tree, select a fixed number of genomes
    at each level based on the sorted quality of the genomes.
    """
    def __init__(self, index, select_number):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_number: int
            The number of genomes to select at each level.

        """
        super().__init__(index)
        self.select_number = select_number
        self.set_all_genomes(boolean=False)

    def select(self, parent, children):
        """
        For each child, sort the genomes by quality.  For the parent,
        take equal portions of genomes from each child.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        i: int
            The number of genomes selected.

        """
        if children:
            children_genomes = []
            for child in children:
                include, _ = filter_genomes(self.index['taxids'][child].keys(),
                                            self.index)
                child_genomes = genome_sort(include, self.index)
                children_genomes.append(child_genomes)
            my_genomes = select_equal(children_genomes, self.select_number)
            for accession in my_genomes:
                self.index['taxids'][parent][accession] = True
            return len(my_genomes)
        else:
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = genome_sort(include, self.index)
            i = 0
            for accession in my_genomes:
                if i >= self.select_number:
                    break
                self.index['taxids'][parent][accession] = True
                i += 1
            return i


class QualitySortingLeaf(GenomeSelection):
    """
    For a leaf node, sort the genomes by quality and select a certain number
    of them.  Propagate all selected genomes up the tree.
    """
    def __init__(self, index, select_number):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_number: int
            The number of genomes to select at each level.

        """
        super().__init__(index)
        self.select_number = select_number
        self.set_all_genomes(boolean=False)

    def select(self, parent, children):
        """
        Sort the genomes at the leaf nodes by quality.  Choose a fixed
        number of them.  At every higher level, include all the genomes
        from the lower levels.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        # if ncbi.get_rank([parent])[parent] == self.leaf_rank:
        if not children:
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = genome_sort(include, self.index)
            samples = 0
            for accession in my_genomes:
                if samples >= self.select_number:
                    break
                samples += 1
                self.index['taxids'][parent][accession] = True
        else:
            samples = 0
            for child in children:
                include, _ = filter_genomes(self.index['taxids'][child].keys(),
                                            self.index)
                for accession in include:
                    if self.index['taxids'][child][accession]:
                        self.index['taxids'][parent][accession] = True
                        samples += 1
        return samples


class AllGenomes(GenomeSelection):
    """Choose all of the genomes."""
    def __init__(self, index):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.

        """
        self.index = index
        self.set_all_genomes(boolean=True)

    def select(self, parent, children):
        """
        Set all genomes in the index to true.  Filter out some genomes.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf node
            is given when children is [].

        Returns
        -------
        int
            The number of genomes selected.

        """
        include, exclude = filter_genomes(self.index['taxids'][parent].keys(),
                                          self.index)
        for accession in exclude:
            self.index['taxids'][parent][accession] = False
        return len(include)


class GenomeHoldout(GenomeSelection):
    """Include whole genomes into separate train anid test data sets."""
    def __init__(self, index, select_number, random=False):
        """
        Initialize the GenomeSelection class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_number: dict
            The number of genomes to include in test and train
            to select at each level.

        """
        super().__init__(index)
        self.select_number = select_number
        # 1 is train, 2 is test,
        # 3 is validate (if applicable, may not be implemented.).
        self.select_type = 0
        # The number of data categories.  2 means train and test.
        # 3 means train, test, and validate.  (Not implemented.  Probably.)
        self.num_categories = 2
        # If random, select random genomes.  Otherwise sort the genomes.
        self.random = random
        # Set the inclusion to all genomes to False (0).
        self.set_all_genomes(boolean=0)

    def get_genomes(self, genome_list):
        """
        Get a list of genomes in a sorted list or in a random list.

        Parameters
        ----------
        genome_list: a list of genome accessions.

        Returns
        -------
        A list of genome accessions.

        """
        if self.random:
            shuffle(genome_list)
            return genome_list
        else:
            return genome_sort(genome_list, self.index)


class GHLeaf(GenomeHoldout):
    """
    Handle an inner node in genome holdout leaf species strategies.
    All genomes from children are propagated to the parent.
    """
    def inner_node(self, parent, children):
        """
        Propagate all genomes from the children to the parent.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        samples = 0
        for child in children:
            include, _ = filter_genomes(self.index['taxids'][child].keys(),
                                        self.index)
            for accession in include:
                if self.index['taxids'][child][accession]:
                    self.index['taxids'][parent][accession] = self.index[
                        'taxids'][child][accession]
                    samples += 1
        return samples


# class GHLeafS(GenomeHoldout):
#     """
#     Handle an inner node in genome holdout leaf species strategies.
#     All genomes from children are propagated to the parent.
#     """
#
#     def inner_node(self, parent, children):
#         """
#         Propagate all genomes from the children to the parent.
#
#         Parameters
#         ----------
#         parent: int
#             Taxonomic id of the parent.
#         children: iterable
#             An iterable of children taxonomic ids of the parent.  A leaf
#             node is represented by [] or False.
#
#         Returns
#         -------
#         samples: int
#             The number of genomes selected.
#
#         """
#         samples = 0
#         for child in children:
#             include, _ = filter_genomes(self.index['taxids'][child].keys(),
#                                         self.index)
#             for accession in include:
#                 if self.index['taxids'][child][accession]:
#                     self.index['taxids'][parent][accession] = self.index[
#                         'taxids'][child][accession]
#                     samples += 1
#         return samples
#
#
# class GHLeafG(GenomeHoldout):
#     """
#     Handle an inner node in genome holdout leaf species strategies.
#     All genomes from children are propagated to the parent.
#     """
#
#     def inner_node(self, parent, children):
#         """
#         Propagate all genomes from the children to the parent.
#
#         Parameters
#         ----------
#         parent: int
#             Taxonomic id of the parent.
#         children: iterable
#             An iterable of children taxonomic ids of the parent.  A leaf
#             node is represented by [] or False.
#
#         Returns
#         -------
#         samples: int
#             The number of genomes selected.
#
#         """
#         samples = 0
#         for child in children:
#             include, _ = filter_genomes(self.index['taxids'][child].keys(),
#                                         self.index)
#             for accession in include:
#                 if self.index['taxids'][child][accession]:
#                     self.index['taxids'][parent][accession] = self.index[
#                         'taxids'][child][accession]
#                     samples += 1
#         return samples


class GHTree(GenomeHoldout):
    """
    Handle an inner node in genome holdout tree strategies.
    A select number of children genomes will be propagated to the parent.
    """
    def inner_node(self, parent, children):
        """
        Propagate some number of genomes from the children to the parent.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        total_selected = sum(self.select_number)
        genome_label_dict = defaultdict(list)
        # Get all genomes from the child nodes and their labels.
        len_genome_set = 0
        for child in children:
            child_dict = defaultdict(list)
            filt_genomes, _ = filter_genomes(
                self.index['taxids'][child].keys(), self.index)
            for accession in filt_genomes:
                if self.index['taxids'][child][accession]:
                    child_dict[self.index['taxids'][child][accession]].append(
                        accession)
                    len_genome_set += 1
            for key in child_dict:
                genome_label_dict[key].append(self.get_genomes(
                    child_dict[key]))
        # If there are singleton genomes and the genome set is large,
        # randomly assign singleton genomes to another label.
#         if parent == 89373:
#             print(children, len(children))
#             print([len(genome_label_dict[key]) for key in genome_label_dict])
        if len_genome_set >= self.num_categories and genome_label_dict[
                SINGLETON]:
            select_type = 0
            singleton_dict = defaultdict(list)
            samples = [0] * self.num_categories
            single_genomes = self.get_genomes(
                list(chain.from_iterable(genome_label_dict[SINGLETON])))
            freq = list(accumulate(self.select_number))
            freq = [num / total_selected for num in freq]
            for accession in single_genomes:
                if min([
                        samples[i] >= self.select_number[i]
                        for i in range(self.num_categories)
                ]):
                    # Randomly assign remaining genomes to labels.
                    prob = random()
                    for i, f in enumerate(freq):
                        if prob <= f:
                            singleton_dict[i + 1].append(accession)
                            break
                else:  # Fill each category up to its limit.
                    singleton_dict[select_type + 1].append(accession)
                    samples[select_type] += 1
                    for i in range(select_type + 1,
                                   select_type + 1 + self.num_categories):
                        j = i % self.num_categories
                        if samples[j] < self.select_number[j]:
                            select_type = j
                            break
            genome_label_dict[SINGLETON] = []
            for key in singleton_dict:
                genome_label_dict[key].append(
                    self.get_genomes(singleton_dict[key]))
        # Down select genomes from each category
        # and assign a label for the parent node.
#         if parent == 89373:
#             print([len(genome_label_dict[key]) for key in genome_label_dict])
#             c = 0
#             dup = {}
#             for key in genome_label_dict:
#                 for l in genome_label_dict[key]:
#                     for a in l:
#                         if a in dup:
#                             print(a)
#                         else:
#                             dup[a] = 1
#                         c += 1
#             print(len_genome_set, c)
        samples = 0
        for label in genome_label_dict:
            shuffle(genome_label_dict[label])  # Randomize choice of child.
            if label != SINGLETON:
                my_genomes = select_equal(genome_label_dict[label],
                                          self.select_number[label - 1])
            else:
                my_genomes = list(
                    chain.from_iterable(genome_label_dict[SINGLETON]))
            for accession in my_genomes:
                samples += 1
                self.index['taxids'][parent][accession] = label
        return samples


# GHTree with down selection for each child.
# class GHTree(GenomeHoldout):
#     """
#     Handle an inner node in genome holdout tree strategies.
#     A select number of children genomes will be propagated to the parent.
#     """
#
#     def inner_node(self, parent, children):
#         """
#         Propagate some number of genomes from the children to the parent.
#
#         Parameters
#         ----------
#         parent: int
#             Taxonomic id of the parent.
#         children: iterable
#             An iterable of children taxonomic ids of the parent.  A leaf
#             node is represented by [] or False.
#
#         Returns
#         -------
#         samples: int
#             The number of genomes selected.
#
#         """
#         total_selected = sum(self.select_number)
#         genome_label_dict = defaultdict(list)
#         # Get all genomes from the child nodes and their labels.
#         for child in children:
#             len_genome_set = 0
#             child_dict = defaultdict(list)
#             filt_genomes, _ = filter_genomes(self.index['taxids']
#                                              [child].keys(), self.index)
#             for accession in filt_genomes:
#                 if self.index['taxids'][child][accession]:
#                     child_dict[self.index['taxids'][child]
#                                [accession]].append(accession)
#                     len_genome_set += 1
#             # If there are singleton genomes and the genome set is large,
#             # randomly assign singleton genomes to another label
#             # after filling the different categories with genomes.
#             if (len_genome_set > total_selected and
#                 child_dict[SINGLETON]):
#                 samples = [len(child_dict[i + 1])
#                            for i in range(self.num_categories)]
#                 select_type = 0
#                 single_genomes = self.get_genomes(child_dict[SINGLETON])
#                 freq = list(accumulate(self.select_number))
#                 freq = [num / total_selected for num in freq]
#                 for accession in single_genomes:
#                     if min([samples[i] >= self.select_number[i] for
#                             i in range(self.num_categories)]):
#                         # random assignment
#                         prob = random()
#                         for i, f in enumerate(freq):
#                             if prob <= f:
#                                 child_dict[i + 1].append(accession)
#                                 break
#                     child_dict[select_type + 1].append(accession)
#                     samples[select_type] += 1
#                     for i in range(select_type + 1,
#                                    select_type + 1 + self.num_categories):
#                         j = i % self.num_categories
#                         if samples[j] < self.select_number[j]:
#                             select_type = j
#                             break
#
# #
# #
# #                 all_genomes = list(chain.from_iterable(
# #                     [child_dict[key] for key in child_dict]))
# #                 all_genomes = self.get_genomes(all_genomes)
# #
# #
# #                 for accession in all_genomes:
# #                     if min([samples[i] >= self.select_number[i] for
# #                             i in range(self.num_categories)]):
# #                         break
# #                     self.index['taxids'][parent][accession] = select_type + 1
# #                     samples[select_type] += 1
# #                     for i in range(select_type + 1,
# #                                    select_type + 1 + self.num_categories):
# #                         j = i % self.num_categories
# #                         if samples[j] < self.select_number[j]:
# #                             select_type = j
# #                             break
# #
# #
# #
# #
# #
# #                 singleton_dict = defaultdict(list)
# #                 freq = list(accumulate(self.select_number))
# #                 freq = [num / total_selected for num in freq]
# #                 for g_list in child_dict[SINGLETON]:
# #                     for accession in g_list:
# #                         prob = random()
# #                         for i, f in enumerate(freq):
# #                             if prob <= f:
# #                                 singleton_dict[i + 1].append(accession)
#                 child_dict[SINGLETON] = []
# #                 for key in singleton_dict:
# #                     child_dict[key].append(
# #                         self.get_genomes(singleton_dict[key]))
#             for key in child_dict:
#                 genome_label_dict[key].append(
#                     self.get_genomes(child_dict[key]))
#         # Down select genomes from each category
#         # and assign a label for the parent node.
#         samples = 0
#         for label in genome_label_dict:
#             shuffle(genome_label_dict[label])  # Randomize choice of child.
#             if label != SINGLETON:
#                 my_genomes = select_equal(genome_label_dict[label],
#                                           self.select_number[label - 1])
#             else:
#                 my_genomes = list(chain.from_iterable(
#                     genome_label_dict[SINGLETON]))
#             for accession in my_genomes:
#                 samples += 1
#                 self.index['taxids'][parent][accession] = label
#         return samples


class GHSpecies(GenomeHoldout):
    """
    Handle a leaf node in genome holdout species strategies.
    Whole species are marked as in one set or another.
    """
    def leaf_node(self, parent):
        """
        Mark chosen genomes as in a set.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                    self.index)
        my_genomes = self.get_genomes(include)
        samples = 0
        for accession in my_genomes:
            if samples >= self.select_number[self.select_type]:
                break
            self.index['taxids'][parent][accession] = self.select_type + 1
            samples += 1
        self.select_type = (self.select_type + 1) % self.num_categories
        return samples


class GHGenome(GenomeHoldout):
    """
    Handle a leaf node in genome holdout genome strategies.
    Whole genomes are marked as in one set or another.
    """
    def leaf_node(self, parent):
        """
        Mark chosen genomes as in a set.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                    self.index)
        my_genomes = self.get_genomes(include)
        samples = [0] * self.num_categories
        select_type = 0
        if len(my_genomes) == 1:
            for accession in my_genomes:
                self.index['taxids'][parent][accession] = SINGLETON
            return 1
        for accession in my_genomes:
            if min([
                    samples[i] >= self.select_number[i]
                    for i in range(self.num_categories)
            ]):
                break
            self.index['taxids'][parent][accession] = select_type + 1
            samples[select_type] += 1
            for i in range(select_type + 1,
                           select_type + +1 + self.num_categories):
                j = i % self.num_categories
                if samples[j] < self.select_number[j]:
                    select_type = j
                    break
        samples = sum(samples)
        return samples


class GHSpeciesLeaf(GHLeaf, GHSpecies):
    """
    Include whole species into separate train and test data sets.
    Down select only at leaves and propagate the selected genomes up the tree.
    """
    def select(self, parent, children):
        """
        Alternately choose species for the train and test data sets.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        if children:
            samples = self.inner_node(parent, children)
        else:
            samples = self.leaf_node(parent)
        return samples


class GHSpeciesTree(GHTree, GHSpecies):
    """
    Include whole species into separate train and test data sets.
    Select a uniform number of genomes at each level.
    """
    def select(self, parent, children):
        """
        Alternately choose species for the train and test data sets.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        if children:  # Inner node
            samples = self.inner_node(parent, children)
        else:  # Leaf node
            samples = self.leaf_node(parent)
        return samples


class GHGenomeLeaf(GHLeaf, GHGenome):
    """
    Include whole genomes into separate train and test data sets.
    Down select only at leaves and propagate the selected genomes up the tree.
    Some inner nodes will have more genomes than others.
    """
    def select(self, parent, children):
        """
        Choose whole genomes for the train and test data sets.
        Propagate chosen genomes up the tree.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        if children:
            samples = self.inner_node(parent, children)
        else:
            samples = self.leaf_node(parent)
        return samples


class GHGenomeTree(GHTree, GHGenome):
    """
    Include whole genomes in train and test data sets.
    Pass up the genomes in the tree and down select genomes
    at each level to a uniform number.
    """
    def select(self, parent, children):
        """
        Choose genomes for train and test sets at each level.
        Choose genomes up to a fixed number at each taxonomic level.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        if children:  # Inner node
            samples = self.inner_node(parent, children)
        else:  # Leaf node
            samples = self.leaf_node(parent)
        return samples


class MinHashTree(GenomeSelection):
    """
    Use a minhash hierarchical clustering tree to select genomes for a given
    taxonomic id.
    """
    def __init__(self, index, select_number):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_number: int
            The number of genomes to select at each level.

        """
        raise NotImplementedError
        super().__init__(index)
        self.select_number = select_number
        self.set_all_genomes(boolean=False)

    def select(self, parent, children):
        """
        Use a minhash tree to select genomes.

        Parameters
        ----------
        parent: int
            Taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf node
            is represented when children is [].

        Returns
        -------
        samples: int
            The number of genomes selected.

        """
        raise NotImplementedError
