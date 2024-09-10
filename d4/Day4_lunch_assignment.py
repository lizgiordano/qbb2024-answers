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
    tissue_samples.setdefault(key, [ ]) # set default so if it doesn't exist in the list, then it doesn't count it
    tissue_samples[key].append(value)     #initialize dict from key with list to hold samples
fs.close()
#print(tissue_samples)


# question 3
# show column names, skipping the first two columns (don't need name and description)
filename=sys.argv[3] # 3 is the 3rd file in the list on the command line, allows you to load in file
fs=open(filename, mode='r')
fs.readline() # skip row
fs.readline() # skip row
header=fs.realine().rstrip("\n").split("\t") # header is a list, shows column names, split by tabs
header=header[2:] # 2: shows everything after column 2

# question 4 and 5
# pull out the columns for all the tissue types
tissue_columns={}
for tissue, samples in tissue_samples.items(): # .items gives keys,values (tissues,sampleID)
    tissue_columns.setdefault(tissue, [])
    for sample in samples:
        if sample in header:
            position=header.index(sample) # index tells you where in the 'key' do you see 'value'
            tissue_columns[tissue].append(position) # append allows you to add to list
# print(header)
# which tissue types have the largest number of samples
maxValue=0
maxValue=""
for tissue, samples in tissue_samples.items():
    if len(samples) > maxValue: 
        maxValue=len(sample)
        maxValue(key) = tissue # if the sample is bigger than the max value, it becomes the new max value, at the end of the run the max value is the largest num of samples
        print(maxValueKey)
        #muscle skeletal has the most number of samples

minValue=20000000
minValue=""
for tissue, samples in tissue_samples.items():
    if len(samples) < minValue:
        minValue=len(sample)
        minValue(key) = tissue #if the sample is smaller than the min value, it becomes the new min value
        # print(minValueKey) 
        # cells-leukemia cell line has the least number of samples


# question 6
# fix Q1
filename = sys.argv[1]
# open file - did before for loop (which would autoclose it) so have to explicitly close it later
fs = open(filename, mode ='r')
# create dict to hold samples for gene-tissue pairs
relevant_samples = {}
# step through file
for line in fs:
    # split line into fields: remove end of line character and tab
    fields = line.rstrip("\n").split("\t")
    # create key from gene and tissue in () to make tuple so we make it an immutable
    key = (fields[0])
    # initalize dict from key with list to hold samples, use [] to make it a list
    relevant_samples[key] = fields[2]
fs.close()


# Q6
f= open("test_data.gct", "r")
for l in f:
    l.strip().split("\t")
    
    geneName=l[0]
    #print(geneName)

    if geneName in relevant_samples.keys():
        myTissue = relevant_samples[geneName] # find the tissue your gene is expressed in
        print(tissue_columns[myTissue[0]]) # keys of tissue_columns are tissues


# question 7
# saved in R file Day4_lunch_Q7
# I expected highly expressed genes to show higher variability becuase they are expressed in more tissue types.
# Some tissues had differences in variances. For example, the pancreas had the least variability. Perhaps certain tissues showed lower variability becuase it is vital for them to be more regulated. This could also just be because generally as numbers get higher, their range for variance increases.

