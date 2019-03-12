"""
Genome selection strategies.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import sys
from random import shuffle


EXCLUDED_GENOMES = {}


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
        return (accession_info['section'] == 'genbank' and
                accession_info['gbrs_paired_asm'] != '' and
                accession_info['paired_asm_comp'] == 'identical' and
                index['genomes'][accession_info['gbrs_paired_asm']]
                )
    include = []
    exclude = []
    for accession in accessions:
        if 'contig_sum' not in index['genomes'][accession]:
            exclude.append(accession)
            EXCLUDED_GENOMES[accession] = 'the contig_sum does not exist'
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

    ref_genomes = [index['genomes'][accession] for
                   accession in genomes if
                   index['genomes'][accession]['refseq_category'] ==
                   'reference genome']
    repr_genomes = [index['genomes'][accession] for
                    accession in genomes if
                    index['genomes'][accession]['refseq_category'] ==
                    'representative genome']
    other_genomes = [index['genomes'][accession] for
                     accession in genomes if
                     index['genomes'][accession]['refseq_category'] !=
                     'reference genome' and
                     index['genomes'][accession]['refseq_category'] !=
                     'representative genome']
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
            my_genomes, _ = filter_genomes(
                self.index['taxids'][parent].keys(), self.index)
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

    def __init__(self, index, select_number):
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
        # 1 is train, 2 is test, 3 is validate (if applicable).
        self.select_type = 0
        # The number of data categories.  2 means train and test.
        # 3 means train, test, and validate.  (May not work well.)
        self.num_categories = 2
        self.set_all_genomes(boolean=0)


class GenomeHoldoutLeaf(GenomeHoldout):
    """
    Include whole genomes into separate train and test data sets.
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
            samples = 0
            for child in children:
                include, _ = filter_genomes(self.index['taxids'][child].keys(),
                                            self.index)
                for accession in include:
                    if self.index['taxids'][child][accession]:
                        self.index['taxids'][parent][accession] = self.index[
                            'taxids'][child][accession]
                        samples += 1
        else:
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = genome_sort(include, self.index)
            samples = 0
            for accession in my_genomes:
                if samples >= self.select_number[self.select_type]:
                    break
                self.index['taxids'][parent][accession] = self.select_type + 1
                samples += 1
            self.select_type = (self.select_type + 1) % self.num_categories
        return samples


class GenomeHoldoutTree(GenomeHoldout):
    """
    Include whole genomes into separate train and test data sets.
    Down select genomes at each level of the tree.
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
        samples = 0
        if children:  # Inner node
            for i in range(self.num_categories):
                my_category = i + 1
                children_genomes = []
                for child in children:
                    include, _ = filter_genomes(self.index['taxids'][child].
                                                keys(), self.index)
                    include = [accession for accession in include
                               if self.index['taxids'][child][accession] ==
                               my_category]
                    child_genomes = genome_sort(include, self.index)
                    children_genomes.append(child_genomes)
                my_genomes = select_equal(children_genomes,
                                          self.select_number[i])
                for accession in my_genomes:
                    self.index['taxids'][parent][accession] = my_category
        else:  # Leaf node
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = genome_sort(include, self.index)
            for accession in my_genomes:
                if samples >= self.select_number[self.select_type]:
                    break
                self.index['taxids'][parent][accession] = self.select_type + 1
                samples += 1
            self.select_type = (self.select_type + 1) % self.num_categories
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
