# day 2 lunch answers
# excercise 1
# answer 1
cut -f 7  hg38-gene-metadata-feature.tsv | sort | uniq -c
19618 protein coding genes
I would like to learn more about miRNAs

# answer 2
cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort -n
273 ENSG00000168036 has the most go_ids

# Exercise 2 
# answer 1
grep "IG_._gene" genes.gtf isolates IG genes row
grep "IG_._gene" genes.gtf | cut -f 1 | uniq -c   pipe to cut the 1st column (chrm) to see number of pseudogenes on each chromosome
52 chr2
  91 chr14
  16 chr15
   6 chr16
   1 chr21
  48 chr22
to compare to IG pseudogenes
grep -e "IG_._pseudogene" -e "IG_pseudogene" genes.gtf | cut -f 1 | uniq -c | sort      
   1 chr1
   1 chr10
   1 chr18
   1 chr8
   5 chr9
   6 chr15
   8 chr16
  45 chr2
  48 chr22
  84 chr14

# answer 2
The word "psuedogene" appears in more than one columns, so using just "pseudogene" will pull more information than we want. 
use "gene_.*pseudogene.* genes.gtf" to narrow down to only "pseudogene" included in the gene_type column

# answer 3
sed "s/ /\t/g" genes.gtf > gene-tabs.gtf
cut -f 1,4,5,14 gene-tabs.gtf > gene-tabs.bed


