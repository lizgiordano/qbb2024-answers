#!/usr/bin/env python3

#importing packages into the current python environment 
import sys
import numpy
import scipy
# from google - "scipy is a python library used for scientific and technical computing. It builds on the foundational NumPy library"

# use python to simulate sequencing 3X coverage 
# then input python code results into R to plot the coverage and analyze distribution
# Q1.2 3X coverage
# genome_size = 1000000
# read_size = 100
# coverage = 3

# from the pseudocode
# calculate number of reads needed to get the coverage you want for the genome
num_reads = int((genome_size * coverage)/read_size)

#print(num_reads)
# num_reads = 30000


genome_coverage = numpy.zeros(genome_size, int)
# creates an array of genome_size with zeros, make each entry an integer 
# track how many times each base pair has been sequenced

for _ in range(num_reads): # create a for loop to fine number of reads
     # selects a random start position and random integer so it doesn't go keep going past the end of the genome
     start_pos = numpy.random.randint(0, genome_size - read_size + 1) 
      # find the end position
     end_pos = start_pos + read_size
     # Each time a read covers a part of the genome, the coverage values increase by 1.
     genome_coverage[start_pos:end_pos] += 1
     #  if start_pos = 10 and read_size = 20
     #  coverage of positions 10 through 30 by adding 1 to the respective positions in the genome_coverage array.

# numpy.savetxt("genome_coverage_3x.txt", genome_coverage)


# # ## get the range of coverages observed
# needed to plot histogram
# max_coverage = max(genome_coverage) # find the max coverage for all positions in genome_coverage array
# xs Creates a list from 0 to max_coverage
# covert to a list
# xs = list(range(0, max_coverage +1)) 

# ## Get the poisson pmf at each of these
# Poisson probability mass function (PMF) for each value in xs
# probability of different coverages across the genome
# poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# ## Get normal pdf at each of these 
# normal probability density function (PDF) for the coverage values based on the mean and standard deviation 
# normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))




# # step 1.5 10X coverage
# repeat 3X coverage, change to 10X

genome_size = 1000000
read_size = 100
coverage = 10

# #pseudocode
num_reads = int((genome_size * coverage)/read_size)

# print(num_reads)
# # num_reads = 30000

# creates an array of genome_size with zeros, make each entry an integer 
# track how many times each base pair has been sequenced

genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads):
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    # find the end position
    end_pos = start_pos + read_size 
    # Each time a read covers a part of the genome, the coverage values increase by 1.
    genome_coverage[start_pos:end_pos] += 1

# numpy.savetxt("genome_coverage_10x.txt", genome_coverage)

# # ## get the range of coverages observed
# needed to plot histogram
# max_coverage = max(genome_coverage) # find the max coverage for all positions in genome_coverage array
# xs Creates a list from 0 to max_coverage
# covert to a list
# xs = list(range(0, max_coverage +1)) 

# ## Get the poisson pmf at each of these
# Poisson probability mass function (PMF) for each value in xs
# probability of different coverages across the genome
# poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# ## Get normal pdf at each of these 
# normal probability density function (PDF) for the coverage values based on the mean and standard deviation 
# normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))





# step 1.6 30X

genome_size = 1000000
read_size = 100
coverage = 30

#pseudocode
num_reads = int((genome_size * coverage)/read_size)

#print(num_reads)
# num_reads = 30000

# array to track coverage at each position in the genome

genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads):
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1

numpy.savetxt("genome_coverage_30x.txt", genome_coverage)

# ## get the range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage +1))

## Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

## Get normal pdf at each of these 
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))






# step 2.1
# reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

# edges = {}
# k = 3
# for read in reads:
#   for i in range(len(read) - k):  #go through each read, look at 1st 3 letter, see if next read matches, keep track of how many times a pairing shows up
#      kmer1 = read[i: i+k] # kmer is fragment of length k
#      kmer2 = read[i+1: i+1+k] 
#      edges.setdefault((kmer1,kmer2), 0)
#      edges[(kmer1,kmer2)] += 1


# # Step 2.2 done in class

# # step 2.3 print("digraph{")

# for left, right in edges:
#    print(f"{left}->{right}")
# print("}")
# # edge is overlap and the base on either side

    
# # Step 2.2 done in class

# # step 2.4 





