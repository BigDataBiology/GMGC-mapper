# Output

Explaination of the files in the output



## prodigal_out.faa , prodigal_out.fna , gene.coords.gbk
These three files are the output of the prodigal.

prodigal_out.faa is the protein sequence.

prodigal_out.fna is the dna sequence.

gene.coords.gbk is the gene information



## hit_table.tsv :

The results of the queries to the GMGC.

There are five columns in the file.

- query_name: the name/id of the input genome contig
- gene_id: the gene_id with the best hit_score in GMGC
- align_category: there are four different classes of alignment
- gene_dna : the dna sequence of the hitted gene in GMGC
- gene_protein : the protein sequence of the hitted gene in GMGC

#### Align_category

- EXACT : above 95% nucleotide identity with at least 95% coverage

- SIMILAR : above 80% nucleotide identity with at least 80% coverage

- MATCH : above 50% nucleotide identity with at least 50% coverage

- NO MATCH : no match in GMGC

  

## genome_bin.tsv

Times of a genome bin that input genes hitting it

There are two columns in the file.

* genome_bin : the name of genome bins in GMGC
* times_gene_hit : the times of input genes hitting it 



## summary.txt

Summary of the query

