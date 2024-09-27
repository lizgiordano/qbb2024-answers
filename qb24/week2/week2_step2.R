library(tidyverse)
library(ggplot2)

# add in data
snp_data <- read.delim("~/qbb2024-answers/qb24/week2/snp_counts.txt") # pull the snp counts file and read through it
# log2 transform enrichment on the Y-axis
snp_data$log_enrichment <- log2(snp_data$enrichment)

# plot with ggplot
ggplot <- ggplot() +
  geom_line(data=snp_data, aes(x = MAF, y = log_enrichment, color = feature, group = Feature)) +
  # graphs data
  scale_color_colorblind() + # colorblind setting
  labs(title = "SNP Enrichment by MAF and Feature", # for legend
       x = "Minor Allele Frequency (MAF)",
       y = "Log2(Enrichment)",
       color = "Feature") + 
  ggsave(filename=~/qbb2024-answers/qb24/week2/enrichments_plot.pdf)

enrichments_plot.pdf

