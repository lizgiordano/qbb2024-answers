
# exercise 3

library(tidyverse)
library(dplyr)
library(ggplot2)

mean_nucleus_signal <- read.table("~/User/cmdb/qbb2024-answers/qb24/week10/mean_nucleus_signal.txt")

# Plot Nascent RNA 
nRNA <- ggplot(mean_nucleus_signal, aes(gene, nascentRNA)) +
  geom_violin(fill = "green") +
  ggtitle("Nascent RNA Signal") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))  # Rotate axis labels for better readability
ggsave("nascent_violin_plot.pdf", plot = nRNA)

# Plot PCNA Signal
PCNA <- ggplot(mean_nucleus_signal, aes(gene, PCNA)) +
  geom_violin(fill = "blue") +
  ggtitle("PCNA Signal") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("PCNA_violin_plot.pdf", plot = PCNA)

# Plot Ratio
ratio <- ggplot(mean_nucleus_signal, aes(gene, Ratio)) +
  geom_violin(fill = "orange") +
  ggtitle("Ratio of Nascent RNA to PCNA") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("ratio_violin_plot.pdf", plot = ratio)

