import pytest
import sys
import pandas as pd
from Bio import SeqIO
sys.path.append(r'../src')

from alignment import identity_coverage


def test_alignment():
    besthits = pd.read_csv('../test_data/besthits.tsv',sep='\t')
    dna_target = str(besthits[0:1]['dna'].values[0]).strip('\n')
    protein_target = str(besthits[0:1]['protein'].values[0]).strip('\n')
    for seq_record in SeqIO.parse('../test_data/genes_dna.fna', "fasta"):
        dna_query = str(seq_record.seq).strip('\n')
        break

    for seq_record in SeqIO.parse('../test_data/genes_protein.faa', "fasta"):
        protein_query = str(seq_record.seq).strip('\n')
        break

    print(identity_coverage(dna_query,protein_query,dna_target,protein_target))



