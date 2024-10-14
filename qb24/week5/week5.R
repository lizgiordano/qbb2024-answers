# code from live coding session
# install new packages
BiocManager::install("DESeq2")
BiocManager::install("vsn")
install.packages("ggfortify")
install.packages("hexbin")

# Load libraries we'll need
library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(hexbin)
library(ggfortify)



### exercise 3 ###
### step 3.1 ###
# Load tab-separated data file
data = readr::read_tsv('salmon.merged.gene_counts.tsv')

# Change gene names into row names
data = column_to_rownames(data, var="gene_name")

# Get rid of gene id column
data = data %>% select(-gene_id)

# Change data to integers
data = data %>% mutate_if(is.numeric, as.integer)

# Remove low coverage samples
data = data[rowSums(data) > 100,]

# Pull out narrow-midgut-tissues
narrow = data %>% select("A1_Rep1":"P2-4_Rep3")

# Create metadata tibble with tissues and replicate numbers based on sample names
narrow_metadata = tibble(Tissue=as.factor(c("A1", "A1", "A1", "A2-3", "A2-3", "A2-3", "Cu", "Cu", "Cu", "LFC-Fe", "LFC-Fe", "Fe", "LFC-Fe", "Fe", "Fe", "P1", "P1", "P1", "P2-4","P2-4","P2-4")),
                         Replicate=as.factor(c("Replicate1", "Replicate2", "Replicate3", "Replicate1", "Replicate2", "Replicate3", "Replicate1", "Replicate2", "Replicate3",
                                               "Replicate1", "Replicate2", "Replicate1", "Replicate3", "Replicate2", "Replicate3", "Replicate1", "Replicate2", "Replicate3", "Replicate1", "Replicate2", "Replicate3")))


### Step 3.2 ####
# Create a DESeq data object
narrowdata = DESeqDataSetFromMatrix(countData=as.matrix(narrow), colData=narrow_metadata, design=~Tissue)

narrowvstdata=vst(narrowdata)
                  
# Plot variance by average
meanSdPlot(assay(narrowvstdata))

narrowvstPCAdata = plotPCA(narrowvstdata,intgroup=c("Replicate","Tissue"), returnData=TRUE)
ggplot(narrowvstPCAdata, aes(PC1, PC2, color=Tissue, shape=Replicate)) +
  geom_point(size=5) +
  labs(title = "PCA - Midgut Tissues")
  ggsave("~/qbb2024-answers/week5/PCAplot.png")



### step 3.3 ###

# Create PCA data
narrowPcaData = plotPCA(narrowLogdata,intgroup=c("rep","tissue"), returnData=TRUE)

# Plot PCA data
ggplot(narrowPcaData, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)

# Batch-correct data to remove excess variance with variance stabilizing transformation
narrowVstdata = vst(narrowdata)

# Plot variance by average to verify the removal of batch-effects
meanSdPlot(assay(narrowVstdata))

# Perform PCA and plot to check batch-correction
narrowPcaData = plotPCA(narrowVstdata,intgroup=c("rep","tissue"), returnData=TRUE)
ggplot(narrowPcaData, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)

# Convert into a matrix
narrowVstdata = as.matrix(assay(narrowVstdata))

# Find replicate means
combined = narrowVstdata[,seq(1, 21, 3)]
combined = combined + narrowVstdata[,seq(2, 21, 3)]
combined = combined + narrowVstdata[,seq(3, 21, 3)]
combined = combined / 3

# Use replicate means to filter low variance genes out
filt = rowSds(combined) > 1
narrowVstdata = narrowVstdata[filt,]

# Plot expression values with hierarchical clustering
heatmap(narrowVstdata, Colv=NA)

# Perform new hierarchical clustering with different clustering method
distance = dist(narrowVstdata)
Z = hclust(distance, method='ave')

# Plot expression values with new hierarchical clustering
heatmap(narrowVstdata, Colv=NA, Rowv=as.dendrogram(Z))

# Set seed so this clustering is reproducible
set.seed(42)

# Cluster genes using k-means clustering
k=kmeans(narrowVstdata, centers=12)$cluster

# Find ordering of samples to put them into their clusters
ordering = order(k)

# Reorder genes
k = k[ordering]

# Plot heatmap of expressions and clusters
heatmap(narrowVstdata[ordering,],Rowv=NA,Colv=NA,RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])

# Save heatmap
png("heatmap.jpg")
heatmap(narrowVstdata[ordering,],Rowv=NA,Colv=NA,RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])
dev.off()

# Pull out gene names from a specific cluster
genes = rownames(narrowVstdata[k == 9,])

# Same gene names to a text file
write.table(genes, "cluster_genes.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

