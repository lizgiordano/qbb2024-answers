BiocManager::install("DropletUtils")
BiocManager::install("zellkonverter")

library(tidyverse)
library(DropletUtils)
library(zellkonverter)
library(ggthemes)
library(scuttle)
library(scater)
library(scran)

gut = readH5AD("v2_fca_biohub_gut_10x_raw.h5ad")
assayNames(gut) = "counts"
gut = logNormCounts(gut)

print(gut)

# Question 1
# how many genes in the dataset?
  # 13407 
# How many cells?
  # 11788 cells 
# What dimension reduction datasets are present?
  # X-pca, X_tsne, and X_umap 

# Question 2 Inspect the available cell metadata

# How many columns are there in colData(gut)?
colData(gut)
colnames(colData(gut))
# 39 columns
# Which three column names reported by colnames() seem most interesting? Briefly explain why.
  # broad annotation gives the cell type. 
  # sizeFactor can normalize the data so they are comparable between cell types. 
  # n_gene gives the number of genes for each cell so you can identify cells with high or low gene counts


#Plot cells according to X_umap using plotReducedDim() and colouring by broad_annotation













