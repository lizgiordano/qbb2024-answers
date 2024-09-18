library(tidyverse)
#Q1
library(tidyverse)
df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#Q2
df
df <- df %>%
  glimpse()

#Q3
RNAseqOnly <- df %>%
  filter(SMGEBTCHT=="TruSeq.v1")
glimpse(RNAseqOnly)

#Q4
ggplot(data=RNAseqOnly, mapping = aes(x=SMTSD))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Tissue Type")

#Q5
ggplot(data=RNAseqOnly, mapping = aes(x=SMRIN)+
  geom_histogram(bins=30) +
  xlab("RIN")
#unimodal

#Q6
ggplot(data=RNAseqOnly, mapping = aes(x=SMRIN, y=SMTSD))+
  geom_boxplot()
#the cell culture lines and cancer lines have much higher RNA integrity perhaps because those cells are actively proliferating 

#Q7
ggplot(data=RNAseqOnly, mapping = aes(x=SMGNSDTC, y=SMTSD))+
  geom_violin()+
  xlab("Genes per sample")+
  ylab("Tissue Type")
#the testes have more gene expression because transcriptional scanning to protect the germline

#Q8
ggplot(data=RNAseqOnly, mapping = aes(x=SMRIN, y=SMTSISCH))+
  geom_point()+
  xlab("RNA Integrity Number")+
  ylab("Ischemic Time")+
  facet_wrap(SMTSD~.)+
  geom_smooth(method="lm")
#different tissues have differences in distribution, suggesting relationship depends on tissue

#Q9
ggplot(data=RNAseqOnly, mapping = aes(x=SMRIN, y=SMTSISCH, color = SMATSSCR))+
  geom_point()+
  xlab("RNA Integrity Number")+
  ylab("Ischemic Time")+
  facet_wrap(SMTSD~.)+
  geom_smooth(method="lm")
#there are tissue differences in autolysis. For example, autolysis is higher in the colon


