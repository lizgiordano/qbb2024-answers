library(tidyverse)
library(broom)
library(DESeq2)


##### LIVE CODING FROM CLASS #####


# set your working directory to where your data and output will be stored
# setwd("~/cmdb/qbb2024-answers/qb24/week7/")


# STEP 1.1


# load the gene expression counts
counts_df <- read_delim("gtex_whole_blood_counts_downsample.txt")
# move the gene_id column to rownames, so that the contents of the
# tibble is entirely numeric

counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")

# look at first five rows

# counts_df[1:5,]

# load the metadata
metadata_df <- read_delim("gtex_metadata_downsample.txt")
# metadata_df[1:5,]

# move the sample IDs from the first column to row names
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")

# confirm that column and row names are correct 
# check column names grabs column names
# check that the columns of the counts are identical and in the same order as the rows of the metadata
colnames(counts_df)
rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))

# STEP 1.2

# create a DESeq object using the counts and meta data using variables sex, age, and dthhrdy
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ SEX + AGE + DTHHRDY)



# STEP 1.3

# apply VST normalization
# give genes with different expression levels equal weights
# store as variance stabilized data
vsd <- vst(dds)

# plot to visualize the principal component analysis (PCA)
# PCA allows you to group genes and create clusters 

plotPCA(vsd, intgroup = "SEX") + labs(title= "PCA Analysis by Sex")

plotPCA(vsd, intgroup = "AGE") + labs(title = "PCA Analysis by Age")

plotPCA(vsd, intgroup = "DTHHRDY") + labs(title = "PCA Analysis by Cause of Death")


##### 1.3.3 Question:
# What proportion of variance in the gene expression data is explained by each of the first two principal components? 
# 7% and 48%
# Which principal components appear to be associated with which subject-level variables? 
# Cause of death clusters with PC1, suggesting it explains 48% of variance. Age clusters with PC2 (~7%).


## Exercise 2: Perform differential expression analysis

# STEP 2.1 test for differential expression between the sexes

# get data in a format to allow to fit a linear model

# extract the normalized expression data and bind to metadata

#pull variance stabilized values, match the columns and rows of the metadata, convert to tibble
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

#test for differential expression of gene WASH7P by performing linear regression of WASHP7
WASH7P <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()
### question 2.1.1
## 2.1.1: Does WASH7P show significant evidence of sex-differential expression (and if so, in which direction)? Explain your answer.
# There's not much change in expression based on sex because the p value 2.792437e-01 is above 0.05


# 2.1.2: Now repeat your analysis for the gene SLC25A47. 
SLC25A47 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

t_test <- t.test(SLC25A47 ~ SEX, data = vsd_df)
t_test
### question 2.1.2
## Does this gene show evidence of sex-differential expression (and if so, in which direction)? Explain your answer.
# There is significant changes in expression in SLC25A47 based on sex because the p value 2.569926e-02 is below 0.05. 
# The gene is upregulated in males because the mean expression is greater in males (female=3.34, male=3.87)

# STEP 2.2.1

# use DESeq2 to perform differential expression analysis across all genes
dds <- DESeq(dds)

# estimating dispersion and fitting model
# if you only run results(dds) it assumes you only care about the last thing
# need to specify what variable you want

# STEP 2.3
# extract differential expression results for the variable SEX
# padj is adjusted p values, used to find the genes that have the most differences in expression between male and female
# need to convert to a tibble with gene names column
sex_res <- results(dds, name = "SEX_male_vs_female")  %>%
  as_tibble(rownames = "GENE_NAME")

sex_res <- sex_res %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# Question 2.3.2 
# How many genes exhibit significant differential expression between males and females at a 10% FDR?
# 262 genes

# STEP 2.3.3
# load mappings of genes to chromosomes from dropbox download and merge to differential expression results
gene_location= read.delim("gene_locations.txt")
location = left_join(gene_loc,sex_res,by="GENE_NAME") %>% arrange(padj)
# Which chromosomes encode the genes that are most strongly upregulated in males versus females, respectively? 
# Are there more male-upregulated genes or female-upregulated genes near the top of the list? Interpret these results in your own words.
# Genes on the Y chromosome are more strongly upregulated in males
# Genes on the sex chromosomes have the most difference. 
# Male-unregulated genes are near the top probably because they are found on the Y chromosome.

# question 2.3.4
# Examine the results for the two genes (WASH7P and SLC25A47) that you had previously tested with the basic linear regression model in step 2.1. 
# Are the results broadly consistent?
# linear regression showed that WASH7P had no sex differences but SLC25A47 did. 
# The differential expression results were consistent with that because SLC25A47 was in 10% FDR but WASH7P wasn't.

### STEP 2.4
# differential expression by death classification
# repeat step 2.2 with DTHHRDY

DTHHRDY_res <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes")  %>%
  as_tibble(rownames = "GENE_NAME")

DTHHRDY_res <- DTHHRDY_res %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# How many genes are differentially expressed according to death classification at a 10% FDR?
# 16069

# Question 2.4.2
# Interpret this result in your own words. Given your previous analyses, does it make sense that there would be more genes differentially expressed based on type of death compared to the number of genes differentially expressed according to sex?
# The PCA showed death classification attributed to 48% of the variance, so it matches with the differential expression results would show more differentially expressed genes than according to sex.


# Exercise 3 Visualization

# create a volcano plot of results from 2.3 
# include significant data = 10% FDR
# plot threshold lines at log2FoldChange> 1 and -log10(padj) > 1.

ggplot(data = sex_res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & -log10(padj)>1))) +
  geom_text(data = sex_res %>% filter(abs(log2FoldChange) > 2 & -log10(padj) > 10),
            aes(x = log2FoldChange, y = -log10(padj) + 5, label = GENE_NAME), size = 2,) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(title="Differential Expression by Sex") +
  labs(y = expression(-log[10]("p-value")), x = expression(log[2]("fold change"))) +
  geom_line(mapping = aes(1), linetype = "dashed", color = "red")+
  geom_line(mapping = aes(-1), linetype = "dashed", color = "red") +
  geom_line(mapping = aes(y=1), linetype = "dashed", color = "red")


