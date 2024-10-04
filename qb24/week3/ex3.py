#!/usr/bin/env python3
import sys
import numpy
import re


# Step 3.1 parse the VCR file
# pseudocode
# for line in open(<vcf_file_name>):
#     if line.startswith('#'):
#         continue
#     fields = line.rstrip('\n').split('\t')

# for line in open(biallelic.vcf):
#     if line.startswith('#'):
#         continue
#     # only keep header with column info

# for line in open(biallelic.vcf):
#     if line.startswith('#'):
#         header = line.rstrip("\n").split("\t") # strip last line, split by tab 
#         info_pos = header.index("header") # find the position of header column
#         format_pos = header.index("format") # find position of "format" in the header
#         sample_pos = list(range(len(header[:format_pos + 1]),len(header))) 
#         # header format pos +1 refers to pos after format
#         # everything beforre sample column is left out
#         # len gives length of list
#         # len header gives number of columns 
#     else:
#     fields = line.rstrip('\n').split('\t')
    
#     # Extracting allele frequency (AF) from the INFO field
#     allele_freq.append(fields[info_pos].split(";")[fields[info_pos].split(";").index(
#         [i for i in fields[info_pos].split(";") if re.findall(r'^AF=', i)][0])])
    
#     # Extracting depth distribution (DP) for each sample
#     for sample in sample_pos:
#         depth_dist.append(fields[sample].split(":")[fields[format_pos].split(":").index("DP")])

allele_freq = []
depth_dist = []
# Open and read the VCF file
with open("biallelic.vcf") as f:
    for line in f:
        # only keep header with column info
        if line.startswith('#'):
            # keep line with the column names
            if line.startswith('#CHROM'):
                header = line.rstrip("\n").split("\t")
                # Find the position of the "INFO" and "FORMAT" fields
                info_pos = header.index("INFO")
                format_pos = header.index("FORMAT")
                # Get the position of the sample columns
                sample_pos = list(range(format_pos + 1, len(header)))
            continue
        
        # for the sample lines
        fields = line.rstrip('\n').split('\t')
        
        # get allele frequency (AF) from INFO 
        info_field = fields[info_pos].split(";")
        allele_freq.append([item for item in info_field if re.findall(r'^AF=', item)][0])
        
        # get the depth value DP for each sample
        for sample in sample_pos:
            format_field = fields[format_pos].split(":")
            sample_data = fields[sample].split(":")
            dp_index = format_field.index("DP")  # Find DP in FORMAT
            depth_dist.append(sample_data[dp_index])  # Append DP for each sample
# step 3.2
    with open(f"AF.txt", "w") as file: # create AF.txt file
        file.write("Frequency\n") # file contain allele freq
        for frequency in allele_freq:
            file.write(f"{frequency[3:]}\n") # starts at 4th position

    with open(f"DP.txt", "w") as file: # creates DP txt file
        file.write("Depth\n") # containing Depth
        for depth in depth_dist:
            file.write(f"{depth}\n") # Writes each depth value to a file, followed by a new line 