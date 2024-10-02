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
# the expected average depth of coverage? 4X coverage

# exercise 1.4 
# Which sample has the largest file size? Which sample has the smallest file size?
# use du? to check the file sizes