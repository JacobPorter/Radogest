"""
Genome selection strategies.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import os
import pickle
import sys
from collections import Counter, defaultdict
from itertools import accumulate, chain
from random import random, shuffle

import numpy as np
from sklearn.cluster import AgglomerativeClustering

EXCLUDED_GENOMES = {}

# The id to use with Genome holdout genome strategies
# when a taxid has a single genome.
SINGLETON = -1

EPS = 2e-24


class DistRecordError(Exception):
    pass


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
    A sorted list of genome accessions.

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


def select_equal(list_lists, select_amount):
    """
    Select equal elements from a list of lists.

    Parameters
    ----------
    list_lists: list
        A list of lists of elements.  The lists could be unequal in size.
    select_amount: int
        An integer of samples to take.

    Examples
    --------
    list_lists is of form [[a, b, c] [1, 2, 3] [x, y].  If select_amount is 5,
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
    while len(output_list) < select_amount and sum(empty) < len(list_lists):
        i = iterations % len(list_lists)
        j = pos[i]
        if j < len(list_lists[i]):
            output_list.append(list_lists[i][j])
            pos[i] += 1
        else:
            empty[i] = 1
        iterations = i + 1
    return output_list


def select_genomes(genome_list,
                   index,
                   down_select="random",
                   select_amount=None,
                   X=None,
                   mapping=None):
    """
    Get a list of genomes in a sorted list or in a random list.

    Parameters
    ----------
    genome_list: list
        A list of genome accessions.
    index: dict
        The Radogest index structure
    down_select: str
        The type of list ordering to perform.
        'random' puts the list in random order
        'sort' puts the list in sorted order
        'dist' uses genome distances (e.g. Mash)
    select_amount: int
        The number of genomes to select.
        If this is None, then the list will not be cut.
    dist_location: str
        The location of the distance matrices, if any.

    Returns
    -------
    A list of genome accessions.

    """
    def splice(a_list):
        if select_amount is not None and select_amount < len(a_list):
            return a_list[0:select_amount]
        else:
            return a_list

    if down_select.startswith("random"):
        shuffle(genome_list)
        return splice(genome_list)
    elif down_select.startswith("sort"):
        return splice(genome_sort(genome_list, index))
    elif down_select.startswith("dist"):
        return cluster(genome_list, select_amount, X, mapping)[0]
    else:
        raise NotImplementedError


def cluster(genome_list, n_clusters, X, mapping):
    """
    Get genomes representing each cluster from a taxid.
    The genomes are sorted in descending order of cluster size.

    Parameters
    ----------
    genome_list: list
            A list of genome accessions.
    taxid: int
        The taxonomic id to get clusters for.
    n_clusters: int
        The number of clusters to request.
    dist_location:
        The directory where distances are stored.

    Returns
    -------
    genome_list: list
        A list of up to n_clusters genomes.
    labels: numpy array
        A list of cluster labels

    """
    # If there is only one genome, just include it.
    # Include all genomes if the number of clusters is too high.
    if len(genome_list) == 1 or n_clusters > len(genome_list):
        return (genome_list, np.zeros((1, 1)))
    # Cluster the genomes.
    labels = AgglomerativeClustering(n_clusters=n_clusters,
                                     affinity="precomputed",
                                     linkage="average").fit(X).labels_
    # Find a genome that best represents the cluster
    # by finding the genome that has the least distance from
    # all other genomes.
    id_label = defaultdict(list)
    genome_list = []
    for i, label in enumerate(labels):
        id_label[label].append(i)
    labels_count = Counter(labels)
    for label in id_label:
        rows = {}
        shuffle(id_label[label])
        for i in id_label[label]:
            rows[i] = sum([X[i][j] for j in id_label[label]])
        min_value = min(rows.values())
        min_pos = [pos for pos in rows if abs(rows[pos] - min_value) <= EPS][0]
        genome_list.append((mapping[min_pos], label, labels_count[label]))
    # Sort the genomes in descending order of cluster size.
    genome_list.sort(key=lambda x: x[2], reverse=True)
    genome_list = [item[0] for item in genome_list]
    return genome_list, labels


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


class StandardSelect(GenomeSelection):
    """

    """
    def __init__(self, index, select_amount, down_select):
        """
        Initialize the StandardSelect class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_amount: int
            The number of genomes to select at each level.
        down_select: str
            A string indicating which down selection to use.
            e.g. 'random' or 'sort'

        """
        super().__init__(index)
        self.select_amount = select_amount
        self.down_select = down_select
        self.set_all_genomes(boolean=False)


class LeafSelect(StandardSelect):
    """
    Perform down selection only at the leaves of the tree.
    """
    def select(self, parent, children):
        """
        Perform down selection only at the leaves of the tree.

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
        if not children:
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = select_genomes(include,
                                        self.index,
                                        down_select=self.down_select,
                                        select_amount=self.select_amount)
            samples = 0
            for accession in my_genomes:
                if samples >= self.select_amount:
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


class TreeSelect(StandardSelect):
    """
    Perform down selection at each level of the tree.
    """
    def select(self, parent, children):
        """
        Perform down selection at each level of the tree.

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
                child_genomes = select_genomes(include,
                                               self.index,
                                               down_select=self.down_select,
                                               select_amount=None)
                children_genomes.append(child_genomes)
            my_genomes = select_equal(children_genomes, self.select_amount)
            for accession in my_genomes:
                self.index['taxids'][parent][accession] = True
            return len(my_genomes)
        else:
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
            my_genomes = select_genomes(include,
                                        self.index,
                                        down_select=self.down_select,
                                        select_amount=None)
            i = 0
            for accession in my_genomes:
                if i >= self.select_amount:
                    break
                self.index['taxids'][parent][accession] = True
                i += 1
            return i
        
        
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
        super().__init__(index)
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


class TreeDistSuper(GenomeSelection):
    """Choose genomes with maximum distance in the taxonomic tree."""
    def __init__(self, index, select_amount, down_select, dist_location):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_amount: int
            The amount of genomes to select.
        down_select: str
            Down selection method: random, sort

        """
        super().__init__(index)
        self.select_amount = select_amount
        self.down_select = down_select
        self.dist_location = dist_location
        self.set_all_genomes(boolean=False)

    def _merge(self, l_lists):
        """
        Merge a list of lists into one list.  Choose elements so that
        the first element of the first list is chosen first, and then the
        first element of the second list is chosen second, and so on.

        Parameters
        ----------
        l_lists: list
            A list of lists to be merged.

        Returns
        -------
            A list that represents merged elements.

        """
        def cond(self, pos):
            for i, l in enumerate(l_lists):
                if pos[i] < len(l):
                    return True
            return False

        new_list = []
        pos = [0] * len(l_lists)
        while cond(self, pos):
            for i, l in enumerate(l_lists):
                if pos[i] < len(l):
                    new_list.append(l[pos[i]])
                    pos[i] += 1
        return new_list

    def _genomes_from_mapping(self, taxid):
        dist_taxid = os.path.join(self.dist_location, str(taxid))
        if not os.path.isdir(dist_taxid):
            raise DistRecordError("The distance records for {} "
                                  "could not be found".format(taxid))
        dist_taxid += "/"
        mapping = pickle.load(open(dist_taxid + str(taxid) + "_map.pck", "rb"))
        X = np.load(dist_taxid + str(taxid) + "_X.npy")
        return (X, mapping)

    def _get_clustered_genomes(self, parent, the_select_amount):
        """Get some genomes by clustering them."""
        try:
            X, mapping = self._genomes_from_mapping(parent)
            include = list(mapping.values())
        except DistRecordError:
            X, mapping = None, None
            include, _ = filter_genomes(self.index['taxids'][parent].keys(),
                                        self.index)
        my_genomes = select_genomes(include,
                                    self.index,
                                    down_select=self.down_select,
                                    select_amount=the_select_amount,
                                    X=X,
                                    mapping=mapping)
        return my_genomes


class TreeDist(TreeDistSuper):
    """Choose genomes with maximum distance in the taxonomic tree."""
    def __init__(self, index, select_amount, down_select, dist_location):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_amount: int
            The amount of genomes to select.
        down_select: str
            Down selection method: random, sort

        """
        super().__init__(index, select_amount, down_select, dist_location)

    def select(self, parent, children, genome_lists):
        """
        Select genomes to include.

        Parameters
        ----------
        parent: int
            The taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.
        genome_lists: list
            A list of lists of genome accessions.

        Returns
        -------
            A list of genome accessions.

        """
        if children:  # Inner node
            shuffle(genome_lists)
            my_genomes = self._merge(genome_lists)[:self.select_amount]
            for accession in my_genomes:
                self.index['taxids'][parent][accession] = True
            return my_genomes
        else:  # Leaf (species) node
            my_genomes = self._get_clustered_genomes(parent,
                                                     self.select_amount)
            for accession in my_genomes:
                try:
                    self.index['taxids'][parent][accession] = True
                except TypeError:
                    print(parent, accession, my_genomes)
                    raise
            return my_genomes


class GenomeHoldout(GenomeSelection):
    """Include whole genomes into separate train anid test data sets."""
    def __init__(self, index, select_amount, down_select=False):
        """
        Initialize the GenomeSelection class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_amount: dict
            The number of genomes to include in test and train
            to select at each level.

        """
        GenomeSelection.__init__(self, index)
        self.select_amount = select_amount
        # 1 is train, 2 is test,
        # 3 is validate (if applicable, may not be implemented.).
        self.select_type = 0
        # The number of data categories.  2 means train and test.
        # 3 means train, test, and validate.  (Not implemented.  Probably.)
        self.num_categories = 2
        # Set the down selection method.  i.e. 'sort' or 'dist'
        self.down_select = down_select
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
        return select_genomes(genome_list, self.index, self.down_select, 
                              None, None, None)


class GHTreeDist(GenomeHoldout, TreeDistSuper):
    """Apply full genome holdout using the TreeDist traversal method. """
    def __init__(self, index, select_amount, down_select, dist_location):
        """
        Initialize the class.

        Parameters
        ----------
        index: dict
            A dictionary representing the index of genomes.
        select_amount: int
            The amount of genomes to select.
        select_type: str
            Down selection method: random, sort

        """
        GenomeHoldout.__init__(self, index, select_amount,
                               0 if down_select == "sort" else 1)
        TreeDistSuper.__init__(self, index, select_amount, down_select,
                               dist_location)

    def select(self, parent, children, genome_lists):
        """
        Select genomes to include.

        Parameters
        ----------
        parent: int
            The taxonomic id of the parent.
        children: iterable
            An iterable of children taxonomic ids of the parent.  A leaf
            node is represented by [] or False.
        genome_lists: list
            A list of lists of lists of genome accessions.
            Each inner list has three lists for train, test, and singleton
            [[[] [] []] [[] [] []] [[] [] []]]

        Returns
        -------
            A list of genome accessions.

        """
        if children:  # Inner node
            shuffle(genome_lists)
            the_end_list = []
            for i in range(3):
                the_list = [item[i] for item in genome_lists]
                my_genomes = self._merge(the_list)[:self.select_amount[i % 2]]
                for accession in my_genomes:
                    self.index['taxids'][parent][
                        accession] = i + 1 if i < 2 else SINGLETON
                the_end_list.append(my_genomes)
            return the_end_list
        else:  # Leaf (species) node
            my_genomes = self._get_clustered_genomes(parent,
                                                     sum(self.select_amount))
            samples = [0] * self.num_categories
            select_type = 0
            if len(my_genomes) == 1:
                for accession in my_genomes:
                    self.index['taxids'][parent][accession] = SINGLETON
                return [[], [], my_genomes]
            selected_genomes = [[], [], []]
            for accession in my_genomes:
                if min([
                        samples[i] >= self.select_amount[i]
                        for i in range(self.num_categories)
                ]):
                    break
                self.index['taxids'][parent][accession] = select_type + 1
                selected_genomes[select_type].append(accession)
                samples[select_type] += 1
                for i in range(select_type + 1,
                               select_type + +1 + self.num_categories):
                    j = i % self.num_categories
                    if samples[j] < self.select_amount[j]:
                        select_type = j
                        break
            samples = sum(samples)
            return selected_genomes


class GHLeaf(GenomeHoldout):
    """
    Handle an inner node in genome holdout leaf strategies.
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
        total_selected = sum(self.select_amount)
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
        if len_genome_set >= self.num_categories and genome_label_dict[
                SINGLETON]:
            select_type = 0
            singleton_dict = defaultdict(list)
            samples = [0] * self.num_categories
            single_genomes = self.get_genomes(
                list(chain.from_iterable(genome_label_dict[SINGLETON])))
            freq = list(accumulate(self.select_amount))
            freq = [num / total_selected for num in freq]
            for accession in single_genomes:
                if min([
                        samples[i] >= self.select_amount[i]
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
                        if samples[j] < self.select_amount[j]:
                            select_type = j
                            break
            genome_label_dict[SINGLETON] = []
            for key in singleton_dict:
                genome_label_dict[key].append(
                    self.get_genomes(singleton_dict[key]))
        # Down select genomes from each category
        # and assign a label for the parent node.
        samples = 0
        for label in genome_label_dict:
            shuffle(genome_label_dict[label])  # Randomize choice of child.
            if label != SINGLETON:
                my_genomes = select_equal(genome_label_dict[label],
                                          self.select_amount[label - 1])
            else:
                my_genomes = list(
                    chain.from_iterable(genome_label_dict[SINGLETON]))
            for accession in my_genomes:
                samples += 1
                self.index['taxids'][parent][accession] = label
        return samples


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
            if samples >= self.select_amount[self.select_type]:
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
                    samples[i] >= self.select_amount[i]
                    for i in range(self.num_categories)
            ]):
                break
            self.index['taxids'][parent][accession] = select_type + 1
            samples[select_type] += 1
            for i in range(select_type + 1,
                           select_type + +1 + self.num_categories):
                j = i % self.num_categories
                if samples[j] < self.select_amount[j]:
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
