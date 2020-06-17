from gmgc_mapper.main import split_file
from os import path
import pytest

def test_split(tmpdir):
    faa_file = path.join(tmpdir, 'input.faa')
    fna_file = path.join(tmpdir, 'input.faa')

    # Both faa_file & fna_file will have 6 elements:
    with open(faa_file, 'wt') as out_a:
        with open(fna_file, 'wt') as out_n:
            for i in range(6):
                out_a.write(f'>g{i}\nMEPTAP\n')
                out_n.write(f'>g{i}\nATGTTC\n')

    odir = path.join(tmpdir, 'output')
    assert split_file(faa_file, odir, is_dna=False ,max_size=1) == 6
    assert split_file(fna_file, odir, is_dna=True, max_size=1) == 6
    assert split_file(faa_file, odir, is_dna=False, max_size=2) == 3
    assert split_file(fna_file, odir, is_dna=True, max_size=2) == 3
    assert split_file(faa_file, odir, is_dna=False, max_size=3) == 2
    assert split_file(fna_file, odir, is_dna=True, max_size=3) == 2
    assert split_file(faa_file, odir, is_dna=False, max_size=6) == 1
    assert split_file(fna_file, odir, is_dna=True, max_size=6) == 1


def test_split(tmpdir):
    faa_file = path.join(tmpdir, 'input.faa')
    fna_file = path.join(tmpdir, 'input.faa')

    total = 20

    # Both faa_file & fna_file will have `total` elements:
    with open(faa_file, 'wt') as out_a:
        with open(fna_file, 'wt') as out_n:
            for i in range(total):
                out_a.write(f'>g{i}\nMEPTAP\n')
                out_n.write(f'>g{i}\nATGTTC\n')


    odir = path.join(tmpdir, 'output')
    n = split_file(faa_file, odir, is_dna= False, max_size=11)
    n = split_file(fna_file, odir, is_dna= True, max_size=11)
    total_nt = 0
    total_aa =0
    for it in range(n):
        with open(path.join(odir, 'split_{}.faa'.format(it+1))) as ifile:
            for line in ifile:
                total_aa += (line[0] == '>')
        with open(path.join(odir, 'split_{}.fna'.format(it+1))) as ifile:
            for line in ifile:
                total_nt += (line[0] == '>')
    assert total_nt == total
    assert total_aa == total
