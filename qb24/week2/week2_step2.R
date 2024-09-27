library(tidyverse)
library(ggplot2)

# add in data
snp_data <- readr::read_tsv("~/qbb2024-answers/qb24/week2/snp_counts.txt") # pull the snp counts file and read through it
# log2 transform enrichment on the Y-axis
# snp_data$enrichment <- as.numeric(snp_data$enrichment)

snp_data <- snp_data %>% mutate(log_enrichment = log2(Enrichment)) 

# plot with ggplot
ggplot() +
  geom_line(data = snp_data, aes(x = MAF, y = log_enrichment, color = Feature, group = Feature)) +
  geom_line(linewidth = 1) +  # Add lines for each Feature
  labs(title = "SNP Enrichment by MAF and Feature",  # for legend
       x = "Minor Allele Frequency (MAF)",
       y = "Log2(Enrichment)",
       color = "Feature")
ggsave(filename = path.expand("~/qbb2024-answers/qb24/week2/enrichments_plot.pdf"))
