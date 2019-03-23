import gffutils

gff_in = '/scratch/jsp4cu/Data/GCF_000001405.38_GRCh38.p12_genomic.gff'
gff_out = '/scratch/jsp4cu/Data/GCF_000001405.38_GRCh38.p12_genomic.gff.db'
# db = gffutils.create_db(gff_in, gff_out, id_spec={'gene': 'db_xref'})

db = gffutils.FeatureDB(gff_out)

test_region = db.region(seqid="NC_000001.11", start=11877, end=11974)
for i in test_region:
    print(i, i.featuretupe)

chroms = [i['seqid'] for i in db.execute('SELECT DISTINCT seqid FROM features;')]

def genes_on_chrom(chrom):
    """
    Yield genes on `chrom`, sorted by start position
    """
    for g in db.features_of_type(('gene', 'exon'), order_by='start', limit=(chrom, 0, 1e9)):
        g.strand = '.'
        yield g

def intergenic():
    """
    Yield intergenic features
    """
    for chrom in chroms:
        genes = genes_on_chrom(chrom)
        for intergenic in db.interfeatures(genes):
            yield intergenic

counter = 0
for feature in intergenic():
    counter += 1
    print(feature)
    if counter == 10:
        break
