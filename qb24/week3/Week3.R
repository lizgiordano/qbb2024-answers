library(tidyverse)
library(ggplot2)

# input data
frequency <- read.delim("~/qbb2024-answers/qb24/week3/AF.txt")
depth <- read.delim("~/qbb2024-answers/qb24/week3/DP.txt")

# plot on histogram showing the allele frequency spectrum (distribution) of the variants in the VCF.
ggplot(frequency, aes(x = Frequency)) +
  geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), bins = 11) + # bin to 11
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Allele Frequency Spectrum",
       y = "Proportion of Variants",
       x = "Allele Frequency") 
  ggsave("~/qbb2024-answers/qb24/week3/Frequency.png") # Save the plot 
  

  # Question 3.1: Interpret this figure in two or three sentences in your own words. Does it look as expected? Why or why not? Bonus: what is the name of this distribution?
  # it's a normal distribution. It looks as expected because I would expect there to be a higher frequency at the average coverage.
  
# step 3.2 
# Plot a histogram showing the distribution of read depth
  ggplot(depth, aes(x = Depth)) +
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), bins = 21) + 
    scale_x_continuous(breaks = c(0:20), limits = c(0, 20)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Read Depth",
         y = "Proportion of Variants",   
         x = "Read Depth") +
    ggsave("~/qbb2024-answers/qb24/week3/Depth.png")
  
# Question 3.2: Interpret this figure in two or three sentences in your own words. 
# Does it look as expected? Why or why not? 
# Bonus: what is the name of this distribution?
  ### it is a poisson distribution. I would expect the depths to increase around the average coverage. The depths are lower in regions with less coverage because those regions may have sequencing bias.
  
  
  
  