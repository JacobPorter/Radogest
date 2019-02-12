import gffutils

gff_in = '/scratch/jsp4cu/Data/GCF_000001405.38_GRCh38.p12_genomic.gff'
gff_out = '/scratch/jsp4cu/Data/GCF_000001405.38_GRCh38.p12_genomic.gff.db'
db = gffutils.create_db(gff_in, gff_out, id_spec={'gene': 'db_xref'})

test_region = db.region(seqid="NC_000001.11", start=11877, end=11974)
for i in test_region:
    print(i.featuretype)
