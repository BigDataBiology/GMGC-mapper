import skbio.alignment
from skbio.sequence import DNA
import numpy as np


def num_alignment(query,target):
    num = 0
    for nucl_q , nucl_t in zip(query ,target):
            if nucl_q == nucl_t:
                num += 1
    return num

sw = skbio.alignment.local_pairwise_align_nucleotide\
    (DNA('AAAATTTTCCCCAAAAA'),
     DNA('AAAAA'))
print(sw)
a = str(sw[0][0])
b = str(sw[0][1])
temp = sw[2][1]
target_start , target_end = temp
print(target_start,target_end)


