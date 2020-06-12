import skbio.alignment
from skbio.sequence import DNA,Protein
from .BLOSUM import blosum62,blosum50

def num_alignment(query,target):
    num = 0
    for nucl_q , nucl_t in zip(query ,target):
            if nucl_q == nucl_t:
                num += 1
    return num

def extract_sw(sw):
    query_align = str(sw[0][0])
    target_aligh = str(sw[0][1])
    align = num_alignment(query_align, target_aligh)
    identity = align / len(target_aligh)
    target_start , target_end = sw[2][1]
    align_length = (target_end-target_start+1)
    return identity,align_length

def identity_coverage(dna_query,protein_query,dna_target,protein_target):
    """
        def category(query, dna_seq, protein_seq):
       if identity_coverage(query, dna_seq) >= (0.95, 0.95): return "EXACT"
       if identity_coverage(query, protein_seq) >= (0.8, 0.8): return "SIMILAR"
       if identity_coverage(query, protein_seq) >= (0.5, 0.5): return "MATCH"
       return "NO MATCH"
    """
    if dna_query != '':
        try:
            sw_dna = skbio.alignment.local_pairwise_align_ssw(DNA(dna_query),DNA(dna_target))
        except:
            sw_dna = skbio.alignment.local_pairwise_align_nucleotide(DNA(dna_query), DNA(dna_target))
        dna_identity,align_length = extract_sw(sw_dna)
        dna_coverage = align_length / min(len(dna_query),len(dna_target))
        if dna_identity >= 0.95 and dna_coverage >= 0.95:

            return 'EXACT'

        else:
            try:
                sw_protein = skbio.alignment.local_pairwise_align_ssw(Protein(protein_query), Protein(protein_target),
                                                                      substitution_matrix=blosum62, gap_open_penalty=11,
                                                                      gap_extend_penalty=1)
            except:
                sw_protein = skbio.alignment.local_pairwise_align_protein(Protein(protein_query),
                                                                          Protein(protein_target),
                                                                          substitution_matrix=blosum62,
                                                                          gap_open_penalty=11, gap_extend_penalty=1)
            protein_identity, align_length = extract_sw(sw_protein)
            protein_coverage = align_length / min(len(protein_query), len(protein_target))
            if protein_identity >= 0.8 and protein_coverage >= 0.8:
                return 'SIMILAR'

            if protein_identity >= 0.5 and protein_coverage >= 0.5:
                return 'MATCH'

            return 'NO MATCH'

    else:
        try:
            sw_protein = skbio.alignment.local_pairwise_align_ssw(Protein(protein_query),Protein(protein_target),substitution_matrix = blosum62,gap_open_penalty=11,gap_extend_penalty=1)
        except:
            sw_protein = skbio.alignment.local_pairwise_align_protein(Protein(protein_query),Protein(protein_target),substitution_matrix = blosum62,gap_open_penalty=11,gap_extend_penalty=1)
        protein_identity,align_length = extract_sw(sw_protein)
        protein_coverage = align_length / min(len(protein_query),len(protein_target))
        if protein_identity >= 0.8 and protein_coverage >= 0.8:
            return 'SIMILAR'

        if protein_identity >= 0.5 and protein_coverage >= 0.5:
            return 'MATCH'

        return 'NO MATCH'


