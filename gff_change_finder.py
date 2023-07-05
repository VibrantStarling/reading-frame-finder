from Bio import pairwise2 as pw2
Bio.Align.PairwiseAligner # pairwise2 will be depreciated in a future biopython version
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(
     description="examine the similarity of an old and new gff for one genome",
     epilog="written by Helen Davison")
parser.add_argument('-g1', '--original_gff', \
                    help="The original gff for the genome",
                    required=True)
parser.add_argument('-g2', '--new_gff', \
                    help="The new gff for the genome to compare against the original one",
                    required=True)
parser.add_argument('-f','--genome_fasta', \
                    help="A genome fasta file that corresponds to the liftoff gff",
                    required=True)
args = parser.parse_args()


genome_fasta = args.genome_fasta

# function to load the gff into a pandas dataframe and unpack the attributes column to access the unique IDs
def get_gff_attributes(gff):
    gff_df = pd.read_csv(gff, sep='\t', comment="#", \
                         names=['seqname', 'source', 'feature', 'start','end','score','strand','frame','attribute'])
    attributes = gff_df['attribute'].str.split(';')
    attr_rows = []
    for i in attributes:
        pairs = [tuple(x.split('=')) for x in i]
        attr_rows.append(dict(pairs))

    attr_df = pd.DataFrame(attr_rows)
    gff_with_attributes = pd.concat([gff_df, attr_df], axis=1)

    return gff_with_attributes


def get_matching_genes(original_gff_attributes,new_gff_attributes):
    # get id, start, end, and sequence from the original gff    
    status = []
    percent_match = []
    record_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta,'fasta'))

    def get_alignment_similarity(og_start,og_end,new_start,new_end,genome_fasta):

        for key, value in record_dict.items():
            original_gene = value[og_start-1:og_end]
            new_gene = value[new_start-1:new_end]
            # compare gene similarity
            global_align = pw2.align.globalxx(original_gene,new_gene)
            seq_length = min(len(original_gene), len(new_gene))
            #  acceses third item of global_align[0],('ATGGATCATTGA', '------CATTGA', 6.0, 0, 12)
            #  The third list element (the first number) is the number of matching residues, \
            #  the other two numbers are the beginning and the end of the alignment.
            matches = global_align[0][2] 
            percent_match = (matches/seq_length)*100
        
        return percent_match
    
    for index, row in original_gff_attributes.iterrows():
        if row['ID'] not in list(new_gff_attributes['ID']):
            status.append('missing')
            percent_match.append(np.nan)
        elif row['ID'] in list(new_gff_attributes['ID']):
            if new_gff_attributes.loc[(row['ID'] == new_gff_attributes['ID']) & (row['start'] == new_gff_attributes['start']) & (row['end'] == new_gff_attributes['end']) & (row['feature'] == new_gff_attributes['feature'])].all().any():
                status.append('true match')
                percent_match.append(1)
            else:
                status.append('partial match')
                for new_row in new_gff_attributes[new_gff_attributes['ID']==row['ID']]:
                    og_start,og_end = row[['start','end']]
                    new_start,new_end = new_row[['start','end']]
                    similarity = get_alignment_similarity(og_start,og_end,new_start,new_end,genome_fasta)
                    percent_match.append(similarity)

    original_gff_attributes['match_status'] = status


    return original_gff_attributes

'''THIS IS THE TEST CODE'''

    # find closest start and end id in the new gff
old_ids = original_gff_attributes['ID'].unique()

for index, row in new_gff_attributes.iterrows():
    if row['ID'] in old_ids:
        print(row['ID'])
        # check which row has all matches
        df = (original_gff_attributes[['ID','start','end','feature','attribute']] == row[['ID','start','end','feature','attribute']])
        df.all(axis=1)==True
        # all match
        new_gff_attributes.loc[df[df.all(axis=1)==True].index,'status'] = 'true match'
        new_gff_attributes.loc[df[df.all(axis=1)==True].index,'percent_match'] = 1

        # start is not equal
        # ROWS ARE GENERALLY NOT UNIQUE
        unequal_starts = list(df.loc[(df['start'] == False) & (df['end'] == True) & (df['feature'] == True) & (df['attribute'] == True)].index)
        unequal_ends = list(df.loc[(df['start'] == True) & (df['end'] == False)  & (df['feature'] == True) & (df['attribute'] == True)].index)
        unequal_start_and_end = list(df.loc[(df['start'] == False) & (df['end'] == False) & (df['feature'] == True) & (df['attribute'] == True)].index)


        for i in unequal_starts:
            start_overhang = original_gff_attributes.loc[i,'start'] - new_gff_attributes.loc[i, 'start']
        
        
            end_overhang = original_gff_attributes.loc[i,'end'] - new_gff_attributes.loc[i, 'end']

            




    elif row['ID'] not in new_ids:  
        original_gff_attributes.loc[df[df].index,'status'] = 'missing'
        original_gff_attributes.loc[df[df].index,'percent_match'] = 0 
        

    # get start, end, and sequence in the new gff


# open gffs
original_gff_attributes = get_gff_attributes(args.original_gff)
new_gff_attributes =  get_gff_attributes(args.new_gff)


