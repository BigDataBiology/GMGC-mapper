# Usage

## Commands

* `-i/--input` : path to the input genome file(.fasta/.gz/.bz2).

* `-o/--output` : Output directory (will be created if non-existent).

* `-nt_input` : path to the input DNA gene file(.fasta/.gz/.bz2).

* `-aa_input` : path to the input Protein gene file(.fasta/.gz/.bz2).

The input must contain a genome file or both DNA and Protein gene file.

## Examples

Input is genome sequence.

```bash
genome2gmgc -i input.fasta -o output
```

Input is DNA/protein gene sequence.

```bash
genome2gmgc -nt_input genes.fna -aa_input genes.faa -o output
```

If input is metagenome , you can use [NGLess](https://github.com/ngless-toolkit/ngless) for assemble and gene prediction.

# NGLess

## Install

The recommended way to install NGLess is through [bioconda](http://bioconda.github.io/):

```
conda install -c bioconda ngless 
```

## Assembly and gene prediction

```bash
ngless "0.6"


sample = 'SAMEA2621155.sampled'
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

