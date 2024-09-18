#2
library(tidyverse)

#3
df <- read_tsv ("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )

#4
df %>%
  group_by(SUBJECT) %>%
  summarize(num_of_samples = n()) %>%
  arrange(num_of_samples) 
#K-562 and GTEX_NPJ8 has the most, GTEX-1JMI6 and 1PAR6 have the least

#5
df %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples = n()) %>%
  arrange(num_of_samples) 
#whole blood & muscle have most, kidney & cervix have least

#6
df_NPJ8 <- df %>% filter(SUBJECT=="GTEX-NPJ8") 
df_NPJ8 %>%
  group_by(SMTSD) %>%
  summarize(counts=n()) %>%
  arrange(-counts) 
view(df_NPJ8)
# whole blood has the most samples. Difference between samples (columns 15-20) due to different sequencing methods

#7
SMATSSCR <- df %>%
  filter(!is.na(SMATSSCR)) %>%
  group_by(SUBJECT) %>%
  summarize(mean=mean(SMATSSCR))
sum(SMATSSCR$mean == 0) 
#15 subjects have a mean SMTSSCR of 1
#standard deviation
#scatterplots of different tissue types