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
bedtools subtract -a tempfile -b introns_chr1.bed > tempfile1 #now subtract introns
bedtools subtract -a tempfile1 -b merged_chr1_cCREs.bed > other_chr1.bed #now subtract cCREs



# 2.1 create bash script that loops through each MAF file and each feature bed file 
# use bedtools coverage to find how many SNPs fall within each set of features.

echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt
# create variables for MAF values containing snp files
MAF=("chr1_snps_0.1" "chr1_snps_0.2" "chr1_snps_0.3" "chr1_snps_0.4" "chr1_snps_0.5")

# feature files made in step 1
features=("chr1_exons.bed" "chr1_cCREs.bed" "introns_chr1.bed" "other_chr1.bed")
genome=("genome_chr1.bed")

# loop through MAF files
for MAF in ${MAF[*]}
    do
        echo ${MAF}
        SNP=${MAF}.bed #select the snp file
        bedtools coverage -a $genome -b $SNP > SNPcoverage.txt # remove SNP from whole genome
        coverage_sum=$(awk '{s+=$4}END{print s}' SNPcoverage.txt) # take the SNPs
        bases_sum=$(awk '{s+=$6}END{print s}' SNPcoverage.txt) # take the bases
        # find number of SNPs for MAF per chromosome length
        background=$(echo "$coverage_sum / $bases_sum" | bc -l) # bc gives background density of SNPs, -l gives decimal
        #add second for loop 
        for feature in ${features[*]} # goes through the feature files
            do  
                echo ${feature}
                bedtools coverage -a $feature -b $SNP > feature_coverage.txt # goes through coverage of SNPs on the features
                feature_snps=$(awk '{s+=$4}END{print s}' feature_coverage.txt) # total the SNPs in the feature
                feature_bases=$(awk '{s+=$6} END {print s}' feature_coverage.txt) # total the bases in the feature
                feature_density=$(echo "$feature_snps / $feature_bases" | bc -l) # find SNP density in the feature
                enrichment=$(echo "$feature_density / $background" | bc -l) # find the enrichment of density to bakcground
                echo -e "${MAF}\t${feature}\t${enrichment}" >> snp_counts.txt # calculates the results and saves to new file
            done    


    done        

    # wasn't able to finish, will continue later, ran into too many issues trying to fix error messages