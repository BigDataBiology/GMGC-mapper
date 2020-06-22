# Output of GMGC-mapper

Explanation of the files in the output directory

## Prodigal output

These three files are the output of prodigal (if GMGC-mapper was called in
genome mode)

- `prodigal_out.faa` protein sequence
- `prodigal_out.fna` DNA sequence
- `gene.coords.gbk` gene information in Genebank format


## Hit Table (`hit_table.tsv`)

The results of the queries to the GMGC.

There are five columns in the file.

- `query_name`: the name/id of the input gene
- `unigene_id`: the Unigene with the best score in the GMGC
- `align_category: there are four different classes of alignment (see below)
- `gene_dna`: the DNA sequence of the best hit in GMGC
- `gene_protein`: the protein sequence of the best hit in GMGC

### Alignment category

- `EXACT`: at least 95% nucleotide identity with at least 95% coverage. As
   unigenes in the GMGC represent 95% nucleotide clusterings (species-level
   threshold), this would mean that the query gene would have clustered with
   the GMGC unigene.
- `SIMILAR`: at least 80% amino acid identity with at least 80% coverage.
- `MATCH`: at least 50% amino acid identity with at least 50% coverage.
- `NO MATCH`: no match in GMGC.


## Genome bins (`genome_bin.tsv`)

Genome bins (MAGs) found in the results (and a count of how many genes are
contained in them).

There are two columns in the file.

- `genome_bin`: the name of genome bins in GMGC
- `nr_hits`: the times of input genes hitting it

Note while not all GMGC unigenes are contained in a genome bin, some are
contained in many. Thus, the total counts will not (except by coincidence)
correspond to the number of genes queried.

## Summary (`summary.txt` and `runlog.yaml`)

The file `summary.txt` provides a human-readable summary of the results, while
`runlog.yaml` is a summary of run metadata (as a YaML file, it is both machine
and human-readable).

The file `summary.txt` should be reproducible and running GMGC-mapper twice on
the same input should produce the same results. By design, though,
`runglog.yaml` includes information such as the time when the analysis was run
which is not reproducible.

