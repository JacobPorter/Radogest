"""
Obsolete.
"""

from radogest_config import SAMTOOLS, BEDTOOLS
from radogest_config import GENOMES_NT, GENOMES_AA, GENOMES_CD

def get_paths():
    """Give the locations of the paths."""
    return {'SAMTOOLS': SAMTOOLS,
            'BEDTOOLS': BEDTOOLS,
            'GENOMES_NT': GENOMES_NT,
            'GENOMES_AA': GENOMES_AA,
            'GENOMES_CD': GENOMES_CD}

# Tools and genomes locations.
paths = get_paths()