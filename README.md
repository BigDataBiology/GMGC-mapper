# Genome2gmgc

Command line tool to query input genome to GMGC project. 



## Install

Genome2gmgc requires [prodigal](https://github.com/hyattpd/Prodigal)

Install from source

```bash
python setup.py install
```



## Parameters

-i/--input : path to the input genome file(.fasta/.gz/.bz2).

-o/--output : Output directory (will be created if non-existent).

-nt_input : path to the input DNA gene file(.fasta/.gz/.bz2).

-aa_input : path to the input Protein gene file(.fasta/.gz/.bz2).

The input must contain a genome file or both DNA and Protein gene file.

## Examples

```bash
genome2gmgc -i input.fasta -o output
```

```bash
genome2gmgc -nt_input genes.fna -aa_input genes.faa -o output
```

## Output

The output folder contains :

(1) prodigal_out.faa , prodigal_out.fna , gene.coords.gbk :  output of prodigal.  .faa file means protein sequence predicted by prodigal and .fna file means nucleotide sequence predicted by prodigal.

(2) split_file folder : contains files splitted from original input fasta for the reason that query API can not query more than 50 sequences in one call.

(3) hit_table.tsv : results of the query. There are five columns in the file: query_name,gene_id,align_category,gene_dna,gene_protein

(4) summary.txt : Summary of the query.



## Align_category

EXACT : above 95% nucleotide identity with at least 95% coverage

SIMILAR : above 80% nucleotide identity with at least 80% coverage

MATCH : above 50% nucleotide identity with at least 50% coverage

NO MATCH : no match in GMGC

