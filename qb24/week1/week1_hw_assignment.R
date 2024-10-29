library(tidyverse)
library(ggplot2)

# step 1.3 plot the histogram of coverage across the genome. Overlay the histogram with a Poisson distribution with lambda = 3. 
# Also overlay the distribution with a Normal distribution with a mean of 3 and a std. dev. of 1.73 (which is the square root of 3).



### GOALS:
# use probability distributions to see how much coverage you need to assess the genome
### Calculating frequencies of genome coverage and then compare them to poisson and normal distribution
# renaming the column to coverage specifies what we're looking at to calculate the frequency of each line aka each coverage value
# then store the coverage values to be able to look at distributions. 



# 3x coverage

# compare the data to poisson distribution with a coverage of 3
# organize and group coverage data
# use dpois to find probability in a poisson distribution

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
# calculate poisson probability mass function with a mean of 3
# calculates the probability of a coverage value
poisson_pmf <- dpois(coverage_values, lambda =3)


# calculate normal pdf
# # calculate normal pdf- probability density function
# makes the data to fit a normal distribution
normal_pdf <- dnorm(coverage_values, mean = 3, sd = sqrt(3))
# dnorm() calculates the probability density function (PDF) of a normal distribution for coverage_values, 
# mean of 3 and a standard deviation of sqrt 3

# plot the coverage from the poisson and normal distribution results to visualize the data and see how well it fits the poisson and normal distributions.
# histogram plot of genome coverage with poisson distribution and normal distribution overlayed
# histogram shows how ofetn each coverage value appears
c3x_plot <- ggplot() +
  geom_histogram(data = genome_coverage, 
                 mapping = aes(x = Coverage, fill = "Genome 3x Coverage"),
                 binwidth = 1, color = "black", alpha = 0.5) +
  # add genome coverage as histogram
  geom_line(data = coverage_freq, 
            mapping = aes(x = Coverage, y = poisson_pmf*nrow(genome_coverage), 
                          color = "Poisson Distribution"), size = 1) +

  
# Poisson distribution allows us to see the frequency at each coverage value (basically the depth of coverage of each base). We set lambda to 3, meaning the average coverage value is 3
# Poisson models the probability at each coverage value, but isn't a great representation if there is a lot of variance or a large mean.
  
# add normal distribution as a line plot
# Scale fill manual adds plot legend and assigns color
  
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



# normal distribution works better with a larger mean. And it is easier for statistics in large datasets. 


ggsave(filename = "~/qb24/week1/ex1_3x_cov.png", 
       plot = c3x_plot, width = 10, height = 6, dpi = 300)

c3x_plot

# now you can visualize th coverage and see if it matches the distribution
# by comparing coverage to poisson and normal distribution we can see if the sequencing reads are random or uniformly distributed
# can show us if the sequencing data is reliable


# step 1.5
# do all the same steps but switch to 10X coverage
# have to change lambda to 10 and standard deviation
# compare the data to poisson distribution with a coverage of 10
# organize and group coverage data
# use dpois to find probability in a poisson distribution
# then plot data to see how they fit with a poisson and normal distribution

# import data
genome_coverage <- read.delim("~/qb24/week1/genome_coverage_10x.txt")

# rename the coverage column in the data frame
# make sure you are pulling coverage
colnames(genome_coverage) <- c("Coverage")

# calculate frequencies of different coverages and store as coverage_freq
coverage_freq <- genome_coverage %>%
  group_by(Coverage) %>%
  summarize(frequency = n())

# coverage values for poisson and normal distribution calculations
coverage_values <- coverage_freq$Coverage

# lambda is the mean of the distribution, average coverage
# calculate poisson pmf
poisson_pmf <- dpois(coverage_values, lambda =10)

# calculate normal pdf- probability density function
# makes the data to fit a normal distribution
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
## repeat again but with 30X coverage 
# do all the same steps and logic

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



# basically we organized the coverage data in python and saved as a text file. 
# then upload the text file into R where we can better do statistical analysis
# In R we compared our data to theoretical models: normal and poisson distributions. 
# by comparing to theoretical distribution models, we can see if things stick out about our sequencing data.
# for example is there more variance or strange distributions that suggest something went wrong with the sequencing. 