# Genome2gmgc

Command line tool to query input genome to GMGC project. 



## Install

Genome2gmgc requires [prodigal](https://github.com/hyattpd/Prodigal)

Install from source

```bash
python setup.py install
```



## Parameters

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

If input is metagenome , you can use [NGLess](https://github.com/ngless-toolkit/ngless) for assemble and gene prediction. For more details , you can [read the docs](https://genome2gmgc.readthedocs.io/en/latest/).

## Output

The output folder contains :

(1) prodigal_out.faa , prodigal_out.fna , gene.coords.gbk :  output of prodigal.  .faa file means protein sequence predicted by prodigal and .fna file means nucleotide sequence predicted by prodigal.

(2) hit_table.tsv : results of the query. There are five columns in the file: query_name,gene_id,align_category,gene_dna,gene_protein.

(3) genome_bin.tsv : results of times of a genome bin that genes hitting it.

(4) summary.txt : Summary of the query.



## Align_category

* EXACT : above 95% nucleotide identity with at least 95% coverage

* SIMILAR : above 80% nucleotide identity with at least 80% coverage

* MATCH : above 50% nucleotide identity with at least 50% coverage

* NO MATCH : no match in GMGC

