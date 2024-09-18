#!/usr/bin/env python3

# open file
# skip 2 lines
# split column header by tabs and skip firt two entries
# for each line
    # split line
    # save field 0 into gene IDs
    # save field 1 into gene names
    # save 2+ into expression values


import sys

import numpy

# Question 1

fs = open(sys.argv[1], mode='r')    # sys.argv is a list, [1] is row, mode r is read
fs.readline()
fs.readline()       #skip 2 lines
line = fs.readline()
fields = line.strip("\n").split("\t")           # split column header by tabs  
tissues = fields [2:]       # skip firt two entries
gene_names = []         # create way to hold gene names 
gene_IDs  = []          # create way to hold gene IDs
expression = []         # create way to hold expression values

for line in fs:         # for each line
    fields = line.strip("\n").split("\t")           # split line
    gene_IDs.append(fields[0])          # save field 0 into gene IDs
    gene_names.append(fields[1])          # save field 1 into gene names
    expression.append(fields[2:])       # save 2+ into expression values
fs.close()

# Question 2

tissues=numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs)
gene_names = numpy.array(gene_names)
expression = numpy.array(expression, dtype=float)

# need to tell numpy what type of data expression data is because it is reading it as strings so we need to tell numpy to read as floats

# Question 3 told to skip


# Question 4
expression10 = expression[:10,:]
mean_expression10 = numpy.mean(expression10, axis=1)
# print(mean_expression10)


# Question 5
mean_expression = numpy.mean(expression)
median_expression = numpy.median(expression)
#print(mean_expression.tolist())
#print(median_expression.tolist())
# The mean is 16.558 and median is 0.027. The mean gives a better view of the data because the median is so skewed


# Question 6
expression1=expression+1        # add 1 to everything to avoid zero
log2_expression1 = numpy.log2(expression1)
mean_log2_expression1 = numpy.mean(log2_expression1)
median_log2_expression1 = numpy.median(log2_expression1)
#print(mean_log2_expression1)
#print(median_log2_expression1)
#Now the mean and median values are closer. They now look more similar to the untransformed mean. So transforming the data got rid of the median bias from outliers.

# Question 7 
# numpy.sort returns a sorted array, it doesn't sort an array in place
cp_log2_expression1 = numpy.copy(log2_expression1)
sort_cp_log2_expression1 = numpy.sort(cp_log2_expression1, axis=1)
# print(sort_cp_log2_expression1[0:5,:])
diff_array = sort_cp_log2_expression1[:,-1] -sort_cp_log2_expression1[:,-2]
# print(diff_array)

#Question 8
where_diff_array = numpy.where(diff_array >=10)
print(numpy.shape(where_diff_array))
# 33 genes whose difference between the highest and second highest tissue expression is greater than 10.






