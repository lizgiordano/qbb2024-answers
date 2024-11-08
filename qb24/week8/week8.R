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
plotReducedDim(gut,"X_umap",color="broad_annotation")



# Explore data
# Sum the expression of each gene across all cells
# Create a vector named genecounts by using rowSums() on the counts matrix returned by assay(gut)

genecounts = rowSums(assay(gut)) # finds gene counts to see expression levels
summary(genecounts)
sort(genecounts, TRUE) # TRUE sets in descending order, use to find highest expressed genes

#Question 3: Explore the genecounts distribution (1 pt)
#What is the mean and median genecount according to summary()?
  ## mean=3185 median=254
# What might you conclude from these numbers?
  ## most genes are on the low end of expression but a few have dramatically higher expression that skew the mean
# What are the three genes with the highest expression after using sort()? 
  ## Hsromega, CR45845, and roX1
# What do they share in common?
  ## all 3 are RNAs



# Question 4a: Explore the total expression in each cell across all genes 
# Create a vector named cellcounts using colSums()
cellcounts = colSums(assay(gut))
# Create a histogram of cellcounts using hist()
hist(cellcounts)
# What is the mean number of counts per cell?
mean_counts <- mean(cellcounts)
mean_counts
  ## mean is 3622.082 counts per cell
# How would you interpret the cells with much higher total counts (>10,000)?
  ## you could look at the cell types with >10,000 expression. 
  ## I would expect them to be in really active cell types like immune cells or metabolic cells
  ## If they are in unexpected cells then there might have been a problem with the sample prep, amplification, or sequencing that gives incorrect extra counts in those cells


## Question 4b: Explore the number of genes detected in each cell
# Create a vector named celldetected using colSums() but this time on assay(gut)>0
celldetected = colSums(assay(gut)>0)
# Create a histogram of celldetected using hist()
hist(celldetected)
# What is the mean number of genes detected per cell?
mean_counts <- mean(celldetected)
mean_counts
  ## mean is 1059.392 genes detected per cell
# What fraction of the total number of genes does this represent?
  # total genes = 13407, 
  # 1059/13407
  # Represents about 7.9% of total number of genes 


# explore mitochondrial reads
mito = grep("^mt:",rownames(gut),value=TRUE)
df = perCellQCMetrics(gut,subsets=list(Mito=mito))
df = as.data.frame(df) 
summary(df)
# to confirm mean and detected match earlier answers
colData(gut) <- cbind( colData(gut), df )

# question 5
# plot subsets mito percent vs broad annotation 
plotColData(object = gut, y = "subsets_Mito_percent", x = "broad_annotation") + 
  theme( axis.text.x=element_text( angle=90 ) ) +
  labs(x= "broad annotation", y="subset_Mito_percent",
      title = "Percentage of Mitochondrial Gene Expression by Cell Type")
# Which cell types may have a higher percentage of mitochondrial reads? 
  ## Metabolically active cells, immune cells, stem cells
# Why might this be the case?
  ## because these cells have high energy demands (and the mitochondria is the powerhouse of the cell :) )


# question 6a epithelial cells
coi = colData(gut)$broad_annotation == "epithelial cell"
epi = gut[,coi]
# now plot epi according to X_umap, color by annotation
plotReducedDim(epi,"X_umap",color="annotation") + 
  labs(title = "Epithelial Cells UMAP") 

# ID marker genes in the anterior midgut
marker.info = scoreMarkers( epi, colData(epi)$annotation )
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

# what are the 6 top marker genes?
  # Mal-A6, Men-b, vnd, betaTry, Mal-A1, Nhe2
# what macromolecule does this region of the gut appear to specialize in metabolizing?
  # sugars

# plot expression of the top marker gene across cell types
plotExpression(epi,c("Mal-A6"), x="annotation") + 
  labs(x = "cell type", y="Expression",title = "Expression of Mal-A6") + 
  theme(axis.text.x=element_text( angle=90 ) ) 


# question 7
# repeat for somatic precursors
precursors = colData(gut)$broad_annotation == "somatic precursor cell"
spc = gut[,precursors]
spc.marker.info = scoreMarkers( spc, colData(spc)$annotation )
chosen = spc.marker.info[["intestinal stem cell"]]
ordered <- chosen[order(chosenSPCs$mean.AUC, decreasing=TRUE),]
goi = rownames(ordered)[1:6]
head(goi)

plotExpression(spc,goi, x = "annotation")  + 
  labs(x= "cell type", y="Expression", title = "Expression of Somatic Precursor Cell Markers") + 
  theme( axis.text.x=element_text( angle=90 ) ) 

# Which two cell types have more similar expression based on these markers?
  # Intestinal Stem Cells and enteroblasts
# Which marker looks most specific for intestinal stem cells?
  # DI
  




