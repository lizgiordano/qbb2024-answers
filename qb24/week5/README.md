# exercise 1
# Step 1.1 
# Can you think of a reason why this sample does not match the expected GC content distribution and base content at the beginning of sequences?
## Fastqc is meant for DNA sequences from across the whole genome in all tissues, but our samples are from the gut. So perhaps the gut just tends to be more GC rich.

# step 1.2
# What is the origin of the most overrepresented sequence in this sample? Does it make sense?
## Serine protease genes in Canton-S Drosophila. This does make sense because serine proteases are a big classes of enzymes invlolved in a lot of processes. So it would make sens for them to be overexpressed, espescially in the gut where they'd be more needed for digestion and immune response.

# exercise 2
# If you were to reject any samples with the percentage of unique reads less than 45%, how many samples would you keep?
## 15 samples are above 45%
# can you see the blocks of triplicates clearly? Try adjusting the min slider up. Does this suggest anything about consistency between replicates?
## Sliding the min slider up ceates more clustering. There triplicate blocks cluster closely together suggesting that there is high consistency between replicates. 

# exercise 3
# step 3.3 What does the PCA plot suggest to you about the data? Why?
## Fe replicate 1 and LFC-Fe replicate 3 may be switched.

# step 3.6 Do the categories of enrichments make sense? Why?
## Yes because there's a lot of metabolic and anabolic processes. Which makes sense for the gut. 