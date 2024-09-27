library(tidyverse)
library(ggplot2)

# step 1.3 plot the histogram of coverage across the genome. Overlay the histogram with a Poisson distribution with lambda = 3. 
# Also overlay the distribution with a Normal distribution with a mean of 3 and a std. dev. of 1.73 (which is the square root of 3).

# 3x coverage
# import data
genome_coverage <- read.delim("~/qb24/week1/genome_coverage_3x.txt")

# rename the coverage column in the data frame to "coverage"
colnames(genome_coverage) <- c("Coverage")

# calculate frequencies of different coverages
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>% # Group the data by each coverage value
  summarize(frequency = n()) # counts how many times that coverage value appears

# coverage values for poisson and normal distribution calculations
coverage_values <- coverage_freq$Coverage

# lambda is the mean of the distribution, average coverage
# calculate poisson pmf
poisson_pmf <- dpois(coverage_values, lambda =3)

# calculate normal pdf
normal_pdf <- dnorm(coverage_values, mean = 3, sd = sqrt(3))
# dnorm() calculates the probability density function (PDF) of a normal distribution for coverage_values, 
# mean of 3 and a standard deviation of sqrt 3


# histogram plot of genome coverage with poisson distribution and normal distribution overlayed
c3x_plot <- ggplot() +
  geom_histogram(data = genome_coverage, 
                 mapping = aes(x = Coverage, fill = "Genome 3x Coverage"),
                 binwidth = 1, color = "black", alpha = 0.5) +
  # add genome coverage as histogram
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = poisson_pmf*nrow(genome_coverage), 
                          color = "Poisson Distribution"), size = 1) +
  
  # add normal distribution as a line plot
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = normal_pdf*nrow(genome_coverage), 
                          color = "Normal Distribution"), size = 1) +
  labs(title="Genome 3x Coverage Distribution with Poisson \n and Normal Distribution", 
       x = "Coverage", 
       y = "Frequency") +
  scale_fill_manual(name = "Legend", values = c("Genome 3x Coverage" = "pink")) +
  # adds histogram to the legend
  scale_color_manual(name = "Legend", values = c("Poisson Distribution" = "blue",
                                                 "Normal Distribution" = "green"))

ggsave(filename = "~/qb24/week1/ex1_3x_cov.png", 
       plot = c3x_plot, width = 10, height = 6, dpi = 300)

c3x_plot





# step 1.5

# import data
genome_coverage <- read.delim("~/qb24/week1/genome_coverage_10x.txt")

# rename the coverage column in the data frame
colnames(genome_coverage) <- c("Coverage")

# calculate frequencies of different coverages
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# coverage values for poisson and normal distribution calculations
coverage_values <- coverage_freq$Coverage

# lambda is the mean of the distribution, average coverage
# calculate poisson pmf
poisson_pmf <- dpois(coverage_values, lambda =10)

# calculate normal pdf
normal_pdf <- dnorm(coverage_values, mean = 10, sd = 3.16)

# histogram plot of genome coverage with poisson distribution and normal distribution overlayed
c10x_plot <- ggplot() +
  geom_histogram(data = genome_coverage, 
                 mapping = aes(x = Coverage, fill = "Genome 10x Coverage"),
                 binwidth = 1, color = "black", alpha = 0.5) +
  # add genome coverage as histogram
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = poisson_pmf*nrow(genome_coverage), 
                          color = "Poisson Distribution"), size = 1) +
  
  # add normal distribution as a line plot
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = normal_pdf*nrow(genome_coverage), 
                          color = "Normal Distribution"), size = 1) +
  labs(title="Genome 10x Coverage Distribution with Poisson \n and Normal Distribution", 
       x = "Coverage", 
       y = "Frequency") +
  scale_fill_manual(name = "Legend", values = c("Genome 10x Coverage" = "pink")) +
  # adds histogram to the legend
  scale_color_manual(name = "Legend", values = c("Poisson Distribution" = "blue",
                                                 "Normal Distribution" = "green"))
ggsave(filename = "~/qb24/week1/ex1_10x_cov.png",
  plot = c10x_plot, width = 10, height = 6, dpi = 300)

c10x_plot


# Step 1.6 30X coverage

# import data
genome_coverage <- read.delim("~/qb24/week1/genome_coverage_30x.txt")

# rename the coverage column in the data frame
colnames(genome_coverage) <- c("Coverage")

# calculate frequencies of different coverages
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# coverage values for poisson and normal distribution calculations
coverage_values <- coverage_freq$Coverage

# lambda is the mean of the distribution, average coverage
# calculate poisson pmf
poisson_pmf <- dpois(coverage_values, lambda =30)

# calculate normal pdf
normal_pdf <- dnorm(coverage_values, mean = 30, sd = 5.47)

# histogram plot of genome coverage with poisson distribution and normal distribution overlayed
c30x_plot <- ggplot() +
  geom_histogram(data = genome_coverage, 
                 mapping = aes(x = Coverage, fill = "Genome 10x Coverage"),
                 binwidth = 1, color = "black", alpha = 0.5) +
  # add genome coverage as histogram
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = poisson_pmf*nrow(genome_coverage), 
                          color = "Poisson Distribution"), size = 1) +
  
  # add normal distribution as a line plot
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = normal_pdf*nrow(genome_coverage), 
                          color = "Normal Distribution"), size = 1) +
  labs(title="Genome 30x Coverage Distribution with Poisson \n and Normal Distribution", 
       x = "Coverage", 
       y = "Frequency") +
  scale_fill_manual(name = "Legend", values = c("Genome 30x Coverage" = "pink")) +
  # adds histogram to the legend
  scale_color_manual(name = "Legend", values = c("Poisson Distribution" = "blue",
                                                 "Normal Distribution" = "green"))
ggsave(filename = "~/qb24/week1/ex1_30x_cov.png",
       plot = c30x_plot, width = 10, height = 6, dpi = 300)

c30x_plot






