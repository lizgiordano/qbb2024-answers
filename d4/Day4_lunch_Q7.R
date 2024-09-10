library(tidyverse)
library(ggthemes)

df <- read.delim("~/qbb2024-answers/d4-afternooon/relevant_expression_values.tsv")
mutant <- df %>% dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID))
mutant <- mutant %>% dplyr::mutate(log_Expression_Values = log2(Expression_Values + 1))

ggplot(data=mutant, mapping=aes(x=log_Expression_Values, y=Tissue_Gene)) + 
  geom_violin() +
  ylab("Tissue, Gene Pair") +
  xlab("Expression Value")

# I expected highly expressed genes to show higher variability becuase they are expressed in more tissue types.
# Some tissues had differences in variances. For example, the pancreas had the least variability. Perhaps certain tissues showed lower variability becuase it is vital for them to be more regulated. This could also just be because generally as numbers get higher, their range for variance increases.