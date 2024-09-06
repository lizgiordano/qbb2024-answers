#!/usr/bin/env python3

import sys

import numpy

#question 1 
# load gene_tissue.tsv file into d4 folder
# wget https://raw.githubusercontent.com/bxlab/cmdb-quantbio/main/assignments/bootcamp/dictionaries/extra_data/gene_tissue.tsv

#question 2
# get file name
filename=sys.argv[2]    # open file
fs = open(filename, mode='r')   # create dict to hold samples fo gene tissue pairs
fs.readline() #skip line
tissue_samples = {}       #step through file
for line in fs:         
    fields = line.rstrip("\n").split("\t")      #split line into fields
    key = fields[6]
    value = fields[0]   #create key from gene and tissue
    tissue_samples.setdefault(key, [])
    tissue_samples[key].append(value)     #initialize dict from key with list to hold samples
fs.close()
#print(tissue_samples)


# question 3
filename=sys.argv[3]
fs=open(filename, mode='r')
fs.readline()
fs.readline()
header=fs.realine().rstrip("\n").split("\t")
header=header[2:]

#question 4
tissue_columns={}
for tissue, samples in tissue_samples.items():
    tissue_columns.setdefault(tissue, [])
    for sample in samples:
        if sample in header:
            position=header.index(sample)
            tissue_columns[tissue].append(position)
print(header)

# question 5