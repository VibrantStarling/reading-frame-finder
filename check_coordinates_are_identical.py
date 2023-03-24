# sequence_similarity
# script to check sequence similarity by coordinate
import numpy as np
import pandas as pd
from Bio import SeqIO

def get_record_by_coord(start,end,contig_fasta):
    """Gets gene feature by coordinate
    
    Args: 
        start (int): A single whole number for the start of gene feature

        end (int): A single whole number for the end of gene feature

        contig_fasta (fasta): The contig to be processed to extract the gene feature from

    Returns:
        List: contains the start and end coordinate as well as the gene feature sequence \
            as a string
    """
    for record in SeqIO.parse(contig_fasta,'fasta'):
        id = record.id
        id_record = record.seq[int(start):int(end)]

    return [start, end, str(id_record).upper()] 


def check_record(coord_list, ncbi_record, tritrypdb_record):
    """Check that gene feratures for two versions of a contig from different assemblies \
        are identical 

    Args:
        coord_list (tsv, txt): a tab delimited file containing two columns for gene \
            feature coordinates from a gff file, the first column are start coordinate \
            and the second columns are end coordinate.

        ncbi_record (fasta): the fasta file containing the first version of a contig

        tritrypdb_record (fasta): the fasta file containing the second version of a contig
         

    Returns:
        DataFrame

    Raises:

    """
    try:
        coords = np.loadtxt(coord_list)

    except FileNotFoundError as fnf_error:
        print(fnf_error)

    assert (coords.dtype == 'int32' == True) & (np.shape(coords)[1]==2),\
        'Please check your file, coords must be a two column tab delimited text file of numbers'
    

    output = []

    for start, end in coords:

        ncbi = get_record_by_coord(start,end,ncbi_record)
        tryp = get_record_by_coord(start,end,tritrypdb_record)
    
        if ncbi[2] == tryp[2]:
            output.append(str(int(start)) + " " + str(int(end)) + " same")

        elif ncbi[2] != tryp[2]:
            output.append(str(int(start)) + " " + str(int(end)) + " different") 
    
    return pd.DataFrame(output)

