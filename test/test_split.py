from genome2gmgc.main import split_file
from os import path

def test_split_file(tmpdir):
    faa_file = path.join(tmpdir, 'input.faa')
    fna_file = path.join(tmpdir, 'input.faa')

    # Both faa_file & fna_file will have 6 elements:
    with open(faa_file, 'wt') as out_a:
        with open(fna_file, 'wt') as out_n:
            for i in range(6):
                out_a.write(f'>g{i}\nMEPTAP\n')
                out_n.write(f'>g{i}\nATGTTC\n')

    odir = path.join(tmpdir, 'output')
    assert split_file(faa_file, fna_file, odir, 1) == 6
    assert split_file(faa_file, fna_file, odir, 2) == 3
    assert split_file(faa_file, fna_file, odir, 3) == 2
    assert split_file(faa_file, fna_file, odir, 4) == 2
