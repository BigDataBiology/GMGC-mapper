# GMGC-Finder

Command line tool to query input genome to GMGC project. 

## Install

GMGC-Finder requires [prodigal](https://github.com/hyattpd/Prodigal)

Install from source

```bash
python setup.py install
```

## Parameters

* `-i/--input`: path to the input genome file(.fasta/.gz/.bz2).

* `-o/--output`: Output directory (will be created if non-existent).

* `-nt_input`: path to the input DNA gene file(.fasta/.gz/.bz2).

* `-aa_input`: path to the input Protein gene file(.fasta/.gz/.bz2).

The input must contain a genome file or both DNA and Protein gene files.

## Examples

Input is genome sequence.

```bash
gmgc-finder -i input.fasta -o output
```

Input is DNA/protein gene sequence(You can just input the protein gene sequences).

```bash
gmgc-finder -nt_input genes.fna -aa_input genes.faa -o output
```
```bash
gmgc-finder -aa_input genes.faa -o output
```

If yout input is a metagenome, you can use
[NGLess](https://github.com/ngless-toolkit/ngless) for assembly and gene
prediction. For more details, [read the
docs](https://gmgc-finder.readthedocs.io/en/latest/usage/).

## Output

The output folder will contain

1. Outputs of gene prediction (prodigal).
2. Complete data table, listing all the hits in GMGC, per gene.
3. Complete table, listing all the genome bins (MAGs) that are found in the results.
4. Human readable summary.

For more details, [read the
docs](https://genome2gmgc.readthedocs.io/en/latest/output/). A description of
the outputs is also written to output folder for convenience.
