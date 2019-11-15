"""
From a fasta file generated by Radogest, determine the species taxids
and genome accessions for each fasta header.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>

"""
import sys

from SeqIterator.SeqIterator import SeqReader


def get_taxid_genomes(fasta_file, index):
    """
    Get a set of genomes ids, species ids, etc.
    from a fasta file generated by Radogest.

    Parameters
    ----------
    fasta_file: str
        The location of the fasta file.
    index: dict
        The Radogest index data structure

    Returns
    -------
    genomes, samples, species
        A set of genome ids.
        A set of taxids that were used for sampling
        A set of species taxids.

    """
    reader = SeqReader(fasta_file, file_type='fasta')
    genomes = set()
    samples = set()
    species = set()
    for record in reader:
        genome_id = record[0].split(':')[0]
        sample_taxid = int(record[0].split(':')[1])
        species_taxid = int(index['genomes'][genome_id]['species_taxid'])
        genomes.add((genome_id, sample_taxid, species_taxid))
        samples.add(sample_taxid)
        species.add(species_taxid)
        print("{}\t{}\t{}\t{}".format(record[0], genome_id, sample_taxid,
                                      species_taxid),
              file=sys.stdout)
    return genomes, samples, species
