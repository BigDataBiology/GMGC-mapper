# Usage

## Commands

* `-i/--input` : path to the input genome file(.fasta/.gz/.bz2).

* `-o/--output` : Output directory (will be created if non-existent).

* `--nt-genes` : path to the input DNA gene file(.fasta/.gz/.bz2).

* `--aa-genes` : path to the input Protein gene file(.fasta/.gz/.bz2).

The input must contain a genome file or DNA and Protein gene file or just Protein gene file.

## Examples

1. Input is a genome sequence (`input.fasta`).

```bash
gmgc-mapper -i input.fasta -o output
```

GMGC-mapper will call `prodigal` to predict genes and then process each gene.

2. Input is DNA/protein gene sequences (`genes.fna` and `genes.faa`,
   respectfully).

```bash
gmgc-mapper --nt-genes genes.fna --aa-genes genes.faa -o output
```
```bash
gmgc-mapper --aa-genes genes.faa -o output
```
# Processing metagenomes using NGLess

If your input is metagenome, you can use
[NGLess](https://github.com/ngless-toolkit/ngless) for assembly and gene
prediction and, then, pass the results to GMGC-mapper.


## Install

The recommended way to install NGLess is through [bioconda](https://bioconda.github.io/):

```
conda install -c bioconda ngless 
```

## Assembly and gene prediction

```bash
ngless "1.0"

sample = 'SAMEA2621155'
input = load_mocat_sample(sample)

preprocess(input, keep_singles=False) using |read|:
    read = substrim(read, min_quality=25)
    if len(read) < 45:
        discard

contigs = assemble(input)
write(contigs, ofile='contigs.fna')

orfs = orf_find(contigs)
write(contigs, ofile='orfs.fna')
```

