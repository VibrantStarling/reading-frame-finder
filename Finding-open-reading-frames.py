## Find-open-reading-frames

from Bio import SeqIO
import gffutils
from sqlalchemy import create_engine, select, MetaData, update, and_

# load gff into a database

def get_the_reading_frame(gene_start):
   
    frame = gene_start % 3
   
    return frame

def get_contig_length(genome_fasta, contig_name):

    record_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta,'fasta'))
    seq = record_dict[contig_name]
    length = len(seq)

    return length

# need to make this fit in check_feature()
def reverse_frame_start(liftoff_gff, genome_fasta):
    for gene in database.features_of_type('gene'):
        for i in feature_db.children(gene, featuretype=["CDS"], order_by='start'):
            end = i.stop
            chrom = gene.chrom
            reverse_start = get_contig_length(genome_fasta, chrom) - end
    return reverse_start

# check gene feature is CDSS
# access reading direction
def check_feature(liftoff_gff, genome_fasta):

    # load the gff
    database = gffutils.create_db(liftoff_gff, dbfn = 'test.db', \
                                  merge_strategy="create_unique", \
                                    keep_order=True)
    feature_db = gffutils.FeatureDB("test.db", keep_order=True)
    
    # sqlalchemy time
    engine = create_engine('sqlite:///test.db')
    meta = MetaData()
    meta.reflect(bind=engine)
    conn = engine.connect()
    group_table = meta.tables['features']

    ## replace queries to test (working :) )
    positive_update_statement = update(group_table)\
        .where((group_table.c.featuretype=='CDS') & (group_table.c.strand=='+')).values(frame=group_table.c.start % 3)
    conn.execute(positive_update_statement).fetchall()

    negative_update_statement = update(group_table)\
        .where((group_table.c.featuretype=='CDS') & (group_table.c.strand=='+'))\
            .values(frame=(get_contig_length(genome_fasta, chrom) - group_table.c.end) % 3)
    conn.execute(negative_update_statement).fetchall()

    with engine.begin() as conn:
        result = conn.execute(positive_strand_query)
        result = conn.execute(negative_update_statement)

    ## example queries that work
    positive_strand_query = select(group_table.c.frame).where((group_table.c.featuretype=='CDS') & (group_table.c.start % 3 == 1))
    conn.execute(positive_strand_query).fetchall()
    negative_strand_query = select_statement = select(group_table.c.start)\
        .where((group_table.c.featuretype=='CDS') & ((reverse_start - group_table.c.start) % 3 == 0))



    # need to check if this works for children objects
    # for gene in database.features_of_type('gene'):
    #     for i in feature_db.children(gene, featuretype=["CDS"], order_by='start'):
    #         start = i.start
    #         end = i.stop
    #         orientation = i.strand
    #         chrom = gene.chrom

    #         feature_db.conn
    #         cur = feature_db.execute(
    #             "SELECT frame text FROM features"
    #             "WHERE featuretype='CDS' AND strand='+'"
    #             "REPLACE('.','.','1');"
    #         )
    #         # if +ve then take the start coordinate
    #         if orientation=="+":
    #             reading_frame = i.frame
    #             frame = get_the_reading_frame(start)               
                
                
    #         elif orientation=="-":
    #             reverse_start = get_contig_length(genome_fasta, chrom) - end
    #             get_the_reading_frame(reverse_start)
        
        
        return 





# if -ve then  find contig length, and take away the end coordinate 
# from contig length to give start point on the reverse strand.

# Use % to work out how many time 3 goes into the start coordinate

# designate the remaining 0,1,2 as the reading frame value for the current gene feature


