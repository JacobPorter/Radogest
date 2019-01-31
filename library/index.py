"""Make index."""

import os
import sys
from tqdm import tqdm
from ete3 import NCBITaxa
from library.util import one_loop_stats, faidx_length

ncbi = NCBITaxa()

# May need to update the ncbi taxonomy database.
# ncbi.update_taxonomy_database()


def process_fai(fai_file_fd):
    """
    Use an fai file to update the index.

    Parameters
    ----------
    fai_file_fd: file object
        An open file object for the fai file.

    Returns
    -------
        None if the fai file indicates no contigs.
        Otherwise, returns a dictionary of statistics.

    """
    stat_dictionary = one_loop_stats(fai_file_fd,
                                     faidx_length,
                                     prefix="contig_")
    if stat_dictionary["contig_count"] == 0:
        return None
    else:
        return stat_dictionary
#     contig_lengths = [int(line.split('\t')[1]) for line in fai_file_fd]
#     if len(contig_lengths) == 0:
#         return None, None
#     else:
#         return ({"contig_count": len(contig_lengths),
#                  "base_length": sum(contig_lengths),
#                  "contig_max": max(contig_lengths),
#                  "contig_min": min(contig_lengths),
#                  "contig_avg": sum(contig_lengths) / len(contig_lengths)},
#                  contig_lengths)


def create_initial_index(index, genome_info, verbose=True):
    """
    Add the contents of JSON genome information to an index.

    Parameters
    ----------
    index: dictionary
        a dictionary representing an index of taxonomic ids to genomes
    genome_info: iterable of dictionaries
        each record represents a genome
    verbose: boolean
        if set to True, prints out periodic counts of records processed.

    Returns
    -------
    int
        A count of the records processed.

    """
    count = 0
    for genome_record in tqdm(genome_info):
        group = genome_record[1]
        genome_record = genome_record[0]
        accession = genome_record["assembly_accession"]
        index['genomes'][accession][
            'location'] = '/{section}/{group}/{accession}/'.format(
                section=genome_record['section'],
                group=group, accession=accession)
        for key in genome_record:
            index['genomes'][accession][key] = genome_record[key]
        for taxid in ncbi.get_lineage(int(genome_record["taxid"])):
            if accession not in index['taxids'][taxid]:
                index['taxids'][taxid][accession] = True
        count += 1
    return count


def update_index_path(root_directory, path, files, verbose=0):
    """
    Calculate necessary things from a single path for the index.

    Parameters
    ----------
    root_directory: str
        A string giving the location of the root directory where genomes
        are found.
    path: str
        The path to a directory containing files
    files: iterable
        An iterable of files under the path.
    verbose: int
        If larger than one, print additional output.

    Returns
    -------
    success, contig_dict, contig_length_dist
        A 1 if an fai file was successfully found and a 0 otherwise.
        A dictionary of contig statistics.
        The length distribution of contigs.

    """
    for name in files:
        if name.endswith('.fai'):
            if verbose >= 2:
                sys.stderr.write(name + "\n")
            fai_fd = open(os.path.join(path, name))
            contigs_dict = process_fai(fai_fd)
            if contigs_dict:
                return 1, contigs_dict
    return 0, None


def update_index_root(index, root_directory, verbose=0):
    """
    Add contig counts and fasta base length to the index.

    Parameters
    ----------
    index: dict
        A python dictionary representing indexed genome information.
    root_directory: str
        A string representing the location of the root directory to search.
    verbose: int
        Controls verbosity of output.

    Returns
    -------
    int, index
        A count of the fasta files processed and the index.

    """
    count = {"fai": 0}
    counter = 0
    accessions_to_remove = {}
    for path, _, files in os.walk(root_directory):
        counter += 1
        if counter % 5000 == 0 and verbose >= 1:
            sys.stderr.write("Processed: {}\n".format(counter))
            sys.stderr.flush()
        number, contigs_dict = update_index_path(root_directory,
                                                 path,
                                                 files,
                                                 verbose=verbose)
        accession = os.path.basename(path)
        if number == 0:
            accessions_to_remove[accession] = True
        else:
            if contigs_dict:
                for key in contigs_dict:
                    index['genomes'][accession][key] = contigs_dict[key]
            count["fai"] += number
    index, remove_count = remove_accessions(index, accessions_to_remove,
                                            verbose=verbose)
    count['remove_count'] = remove_count
    sys.stderr.flush()
    return index, count


def remove_accessions(index, accessions_to_remove, verbose=0):
    """
    Remove genomes from the index.

    Parameters
    ----------
    index: dict
        The genomes index.
    accessions_to_remove: dict
        A dict keyed on accessions to remove from the index.

    Returns
    -------
    dict
        The genomes index is returned.

    """
    remove_count = 0
    print("Removing empty genome folders from the index.", file=sys.stderr)
    count = 0
    for taxid in index['taxids']:
        my_accessions_to_remove = []
        for accession in index['taxids'][taxid]:
            if accession in accessions_to_remove:
                my_accessions_to_remove.append(accession)
        for accession in my_accessions_to_remove:
            del index['taxids'][taxid][accession]
            remove_count += 1
        count += 1
        if verbose >= 1 and count % 5000 == 0:
            print("Processed: {}".format(count), file=sys.stderr)
            sys.stderr.flush()
    return index, remove_count
