#!usr/bin/bash
# week2 
# conda activate qb24

#1.4
# Loop through each BED file in exon_chr1_bed 
files=("chr1_exons.bed" "chr1_genes.bed" "chr1_cCREs.bed")
for I in "${files[@]}"
    do
        bedtools sort -i ${I} > sorted.bed # save sorted file
        bedtools merge -i sorted.bed > merged_${I} #merge sorted file 
    done

# Merge all the individual merged BED files into one final file
#cat merged_* > final_exon_chr1_bed.bed
#bedtools merge -i final_exon_chr1_bed.bed > final_merged_exon_chr1_bed.bed

#1.5 Use bedtools to find intervals not covered by other features
#bedtools subtract - takes overlap between two sets and removes the overlap regions
# format = bedtools subtract -a 1st set -b 2nd set > new file
# remove overlap between whole file and exons -> get introns
bedtools subtract -a merged_chr1_genes.bed -b merged_chr1_exons.bed > introns_chr1.bed 
#substracts b from a (exons from genes) and saves new file as introns

#1.6
#do the same thing but subtract exons, introns, cCREs to get other
bedtools subtract -a genome_chr1.bed -b merged_chr1_exons.bed > tempfile #genome minus exons
bedtools subtract -a tempfile -b introns_chr1.bed > tempfile #now subtract introns
bedtools subtract -a tempfile -b merged_chr1_cCREs.bed > other_chr1.bed #now subtract cCREs



# 2.1 create bash script that loops through each MAF file and each feature bed file 
# use bedtools coverage to find how many SNPs fall within each set of features.

echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt
SNP=("chr1_snps_0.1" "chr1_snps_0.2" "chr1_snps_0.3" "chr1_snps_0.4" "chr1_snps_0.5")
features=("chr1_exons.bed" "chr1_cCREs.bed" "introns_chr1.bed" "other_chr1.bed")
genome=("genome_chr1.bed")

for SNP in "${SNP[@]}"
    do
        bedtools coverage -a $genome -b $SNP > SNPcoverage.txt # SNP coverage
        coverage_sum=$(awk '{s+=$4}END{print s}' SNPcoverage.txt)
        sum_total_bases=$(awk '{s+=$6}END{print s}' SNPcoverage.txt)
        # find number of SNPs for MAF per chromosome length
        background=$(echo "$coverage_sum / $sum_total_bases" | bc -l)
        #add second for loop 
        for feature in "${feature[@]}"
            do  
                bedtools coverage -a $feature -b $SNP > SNPfeatures.txt

            done    


    done        
    
    # wasn't able to finish, will continue later, ran into too many tech issues with computer