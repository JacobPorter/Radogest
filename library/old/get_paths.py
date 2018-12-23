"""Stores the paths to tools and the data."""

SAMTOOLS = "/home/jsporter/Applications/samtools-1.8/"
UCSC = "/home/jsporter/Applications/UCSC/"
BEDTOOLS = "/home/jsporter/Applications/bedtools2/bin/"
GENOMES = "/groups/fungcat/datasets/current/fasta/Genomes/"
GENOMESAA = "/groups/fungcat/datasets/current/fasta/GenomesAA/"


def get_paths():
    """Give the locations of the paths."""
    return {'SAMTOOLS': SAMTOOLS,
            'UCSC': UCSC,
            'BEDTOOLS': BEDTOOLS,
            'GENOMES': GENOMESAA}
