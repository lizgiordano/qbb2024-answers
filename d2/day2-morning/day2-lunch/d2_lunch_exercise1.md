# day 2 lunch answers
# answer 1
cut -f 7  hg38-gene-metadata-feature.tsv | sort | uniq -c
19618 protein coding genes
I would like to learn more about miRNAs
# answer 2
cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort -n
273 ENSG00000168036 has the most go_ids
