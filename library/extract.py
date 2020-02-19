"""
Extract fasta records with a given label or taxonomic id.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>
"""

from seq.SeqIterator import SeqReader, SeqWriter


def extract(fasta_file, taxid_file, fasta_output, taxid_output, extract_list):
    """
    Extract fasta records with labels from extract_list.
    Write to taxid and fasta files.

    Parameters
    ----------
    fasta_file: str
        The input fasta file.
    taxid_file: str
        The input taxid file.
    fasta_output: str
        The location to write the fasta file output.
    taxid_output: str
        The location to write the taxid file output.
    extract_list: list
        A list of labels (taxonomic ids) to extract.

    Returns
    -------
    int, int
        The number of records written followed
        by the number of records processed.

    """
    fasta_iterator = SeqReader(fasta_file, file_type='fasta')
    taxid_input = open(taxid_file)
    fasta_writer = SeqWriter(open(fasta_output, 'w'), file_type='fasta')
    taxid_out_file = open(taxid_output, 'w')
    counter_write = 0
    counter_processed = 0
    for fasta_record, taxid in zip(fasta_iterator, taxid_input):
        taxid = int(taxid)
        counter_processed += 1
        if taxid in extract_list:
            counter_write += 1
            fasta_writer.write(fasta_record)
            print(str(taxid), file=taxid_out_file)
    taxid_input.close()
    taxid_out_file.close()
    return counter_write, counter_processed
