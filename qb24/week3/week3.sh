#!/usr/bin/env bash
# download live coding data.tar
# BWA index c_elegans
# SRR16356854_1.subset.fastq.gz & SRR16356854_2.subset.fastq.gz are one sample


for my_sample in *_1.subset.fastq.gz
do
    my_sample=`basename ${my_sample} _1.subset.fastq.gz`
    bwa index c_elegans.PRJNA13758.WS283.genomic.fa
    bwa mem c_elegans.PRJNA13758.WS283.genomic.fa ${my_sample}_1.subset.fastq.gz ${my_sample}_2.subset.fastq.gz > ${my_sample}.sam
    samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam
    samtools index ${my_sample}.bam
done


# exercise 1
# use the commands given in instructions to decompress the Dropbox files
wget https://www.dropbox.com/s/ogx7nbhlhrp3ul6/BYxRM.tar.gz
tar -xvzf BYxRM.tar.gz #tar to change the files, x to extract, v to display, z to indicate the files are compressed, f tells file name
READLENGTH=$(head -n 2 A01_09.fastq | tail -n 1 | awk '{print length}') 
# use first 2 lines, the name and the sequence, tail -n 1 takes out the 2nd line
# awk prints the sequence read

# Reads are 76 bp long


# exercise 1.2
# how many reads are in the file?

LINECOUNT=$(wc -l < A01_09.fastq)  # get the number of lines, wc -1 counts the number of lines
READCOUNT=$((LINECOUNT / 4))    # Count the lines, divide by 4 for each read in a FASTQ
# each read in a FASTQ file has 4 lines - sequence identifier, sequence, separator, and quality score

# There are 669548 reads in the file

# exercise 1.3 
# Question 1.3: Given your answers to 1 and 2, as well as knowledge of the length of the S. cerevisiae reference genome, what is the expected average depth of coverage?
# calculate coverage by adding number of reads and the readlength and dividing that by total genome length
GENOMELENGTH=$(wc -m < sacCer3.fa) # full genome length, wc -m counts the number of characters
COVERAGE=$(echo "${READCOUNT} * ${READLENGTH} / ${GENOMELENGTH}" | bc -l) # calculate coverage by readlength and readcount / genome length
echo "Coverage: $COVERAGE" # print 
# the expected average depth of coverage? 4.2X coverage

# exercise 1.4 
# Which sample has the largest file size? Which sample has the smallest file size?
# use du? to check the file sizes

du -m A01_{09,11,23,24,27,31,35,39,62,63}.fastq | awk '{print $1 " MB\t" $2}' > fastq_file_sizes.txt
# du -m pulls sizes of the files within {}
less fastq_file_sizes.txt # see file sizes to compare

# A01_62 is the largest (150MB) and A01_27 is the smallest (110MB)

# exercise 1.5
# Run the program FastQC
# FastQC A01_09.fastq
# FastQC A01_11.fastq
# FastQC A01_23.fastq
# FastQC A01_24.fastq
# FastQC A01_27.fastq
# FastQC A01_31.fastq
# FastQC A01_35.fastq
# FastQC A01_39.fastq
# FastQC A01_62.fastq
# FastQC A01_63.fastq
FastQC *.fastq
# median base quality along the read is 37
# The base quality is the accuracy of each base in a sequencing read. It is a log function of Q=−10log(P). So the median base quality means that the median probability of an error is 10e-37
# There is more variation in quality at the end of the read. Which makes sense because that usually has occurred when I have doen sequencing in the past. My sequencing is always worse at the end. 



# Excercise 2

# download and index sacCer3 genome
wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
# using bwa index, create an index for the sacCer3.fa reference
bwa index sacCer3.fa

less sacCer3.fa.amb
# Question 2.1: How many chromosomes are in the yeast genome?
## 17 chromosomes


# follow Rajiv's code from class:
# bwa index c_elegans.PRJNA13758.WS283.genomic.fa

# for my_sample in *_1.subset.fastq.gz
# do
#     my_sample=`basename ${my_sample} _1.subset.fastq.gz`
#     bwa mem c_elegans.PRJNA13758.WS283.genomic.fa ${my_sample}_1.subset.fastq.gz ${my_sample}_2.subset.fastq.gz > ${my_sample}.sam
#     samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam
#     samtools index ${my_sample}.bam
# done




# Question 2.2: How many total read alignments are recorded in the SAM file?

TOTAL_READS=$(echo "$(wc -l < sacCer_A01_09.sam) - 20" | bc)
# echo $TOTAL_READS
## there are 669548 read alignments

# Question 2.3: How many of the alignments are to loci on chromosome III?
grep -w "III" sacCer_A01_09.sam | wc -l
## 18196 alignments in chrom 3

# Step 2.4: Format and index your alignments
# use samtools, sort .sam files using samtools sort
# create an index for each sorted .bam files using samtools index
# should have 10 sorted .bam files and .bam.bai indices

for MYSAMPLE in 09 11 23 24 27 31 35 39 62 63
do
    SAMPLE=A01_${MYSAMPLE}  
    # loop generates filenames based on the list of numeric values (09, 11, 23, etc.). These filenames follow the pattern A01_<value>.fastq
    bwa mem A01_${SAMPLE}.fastq > sacCer_${SAMPLE}.sam
    # aligning sequencing reads against a reference genome. mem is algorithim in bwa for illumina seq alignment
    # aligns the sequencing reads from the FASTQ file {SAMPLE}.fastq to the reference genome. 
    samtools sort -@ 4 -O bam -o sacCer_A01_${SAMPLE}.bam sacCer_A01_${SAMPLE}.sam 
    # samtools sort the SAM file into BAM format (binary version of SAM)
    samtools index sacCer_A01_${SAMPLE}.bam
    # creates an index for the BAM file
done

# Question 2.4: Does the depth of coverage appear to match that which you estimated in Step 1.3? Why or why not?
# Yes it looks like the predicted 4X coverage

# Question 2.5: Set your window to chrI:113113-113343 (paste that string in the search bar and click enter). 
# How many SNPs do you observe in this window? 
# 3
# Are there any SNPs about which you are uncertain? Explain your answer.
# chrI:113,326 because unlike the other 2 it's not covered by all the reads

# Question 2.6: Set your window to chrIV:825548-825931. 
# What is the position of the SNP in this window? 
# chrIV:825,834
# Does this SNP fall within a gene?
# SCC2





