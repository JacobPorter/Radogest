FILE_LOCATION = "/scratch/fungcat/jsporter/Data/Plinko_Root/Binary_OvA/test_set/1/train/1.fasta"
FILE_OUTPUT = "/scratch/fungcat/jsporter/Data/Plinko_Root/Binary_OvA/test_set/1.noN.fasta"
TAXID_OUTPUT = "/scratch/fungcat/jsporter/Data/Plinko_Root/Binary_OvA/test_set/1.noN.taxid"

from SeqIterator import SeqIterator, SeqWriter

fasta_input = SeqIterator(FILE_LOCATION, file_type='fasta')
fasta_output = SeqWriter(open(FILE_OUTPUT,"w"), file_type='fasta')
taxid_output = open(TAXID_OUTPUT, "w")

for record in fasta_input:
    if 'N' not in record[1] and 'n' not in record[1]:
        fasta_output.write(record)
        taxid = record[0].split(":")[1]
        print(taxid, file=taxid_output)
