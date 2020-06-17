# GMGC-mapper

![gmgc_mapper_test](https://github.com/BigDataBiology/GMGC-Finder/workflows/gmgc_mapper_test/badge.svg)


Command line tool to query the Global Microbial Gene Catalog (GMGC).

## Install

GMGC-mapper runs on Python 3.6-3.8 and requires
[prodigal](https://github.com/hyattpd/Prodigal) to be available for genome
mode.

### Install from source

```bash
python setup.py install
```


## Examples

1. Input is a genome sequence.

```bash
gmgc-mapper -i input.fasta -o output
```

2. Input is DNA/protein gene sequences

```bash
gmgc-mapper --nt-genes genes.fna --aa-genes genes.faa -o output
```

The nucleotide input is optional (but should be used if available so that the
quality of the hits can be refined):

```bash
gmgc-mapper --aa-genes genes.faa -o output
```

If yout input is a metagenome, you can use
[NGLess](https://github.com/ngless-toolkit/ngless) for assembly and gene
prediction. For more details, [read the
docs](https://gmgc-mapper.readthedocs.io/en/latest/usage/).

## Output

The output folder will contain

1. Outputs of gene prediction (prodigal).
2. Complete data table, listing all the hits in GMGC, per gene.
3. Complete table, listing all the genome bins (MAGs) that are found in the results.
4. Human readable summary.

For more details, [read the
docs](https://gmgc-mapper.readthedocs.io/en/latest/output/). A description of
the outputs is also written to output folder for convenience.

## Parameters

* `-i/--input`: path to the input genome file(.fasta/.gz/.bz2).

* `-o/--output`: Output directory (will be created if non-existent).

* `--nt-genes`: path to the input DNA gene file(.fasta/.gz/.bz2).

* `--aa-genes`: path to the input Protein gene file(.fasta/.gz/.bz2).

