library(tidyverse)
library(broom)
library(ggplot2)

# Exercise 1
# step 1.1
dnm <- read_csv(file= "~/qbb2024-answers/d4/aau1043_dnm.csv")
head(dnm)

ages <- read_csv(file = "~/qbb2024-answers/d4/aau1043_parental_age.csv")
head(ages)

# step 1.2
dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))
# step 1.3
ages <- read_csv(file = "~/qbb2024-answers/d4/aau1043_parental_age.csv")

# step 1.4
# create combined data frame
dnm_by_parental_age <- left_join(dnm_summary, ages, by = "Proband_id")
print (dnm_by_parental_age)

# Exercise 2
# step 2.1
# 1 maternal de novo mutations vs maternal age
ggplot(data=dnm_by_parental_age, mapping = aes(x=Mother_age, y=n_maternal_dnm))+
  geom_point()
# 2 paternal de novo mutations vs paternal age
ggplot(data=dnm_by_parental_age, mapping = aes(x=Father_age, y=n_paternal_dnm))+
  geom_point()
# step 2.2
# 1 Fit a linear regression model to maternal age and maternal DNM model using lm()
lm(data=dnm_by_parental_age,
   formula = n_maternal_dnm ~ 1 + Mother_age) %>%
  summary()
ggplot(data=dnm_by_parental_age,
       mapping=aes(x=Mother_age, 
                   y=n_maternal_dnm))+
  geom_point()+
  stat_smooth(method = "lm")
# discuss the size and significance of the relationship. 
# The average child has 2.5 DNMs. Approximately 1 mutation is gained for every 3 years the mother ages.
# increase in mother age increases the chance of maternal de novo mutations. This does match with plots in step 2.1
# This is significant and of high size because the p-value is so small (exponent of -16), suggesting there is a relationship between maternal age and mutations.

# 2 
# the relationship is significant because the p-value is so small. So this means that the observed relationship is not due to chance.

# step 2.3
# test for an association between paternal age and paternally inherited de novo mutations.
paternal_model <- lm(data=dnm_by_parental_age, formula = n_paternal_dnm ~ 1 + Father_age)
summary(paternal_model)
ggplot(data=dnm_by_parental_age,
       mapping=aes(x=Father_age, 
                   y=n_paternal_dnm))+
  geom_point()+
  stat_smooth(method = "lm")
# discuss the size and significance of the relationship
# 1 As paternal age increases, the chance of paternal de novo mutations also increases. The average child has 10 dnm, and gains 1 mutation with every year of paternal age.
# 2 the relationship is significant with high size because the p-value is so small (exponent of -16). So this means that the observed relationship, that mutations increase with paternal age increase, is not due to chance.

# Step 2.4
father_a <- tibble(
  paternal_age = 50.5
)
print(predict(paternal_model,father_a))


# Step 2.5
ggplot(data=dnm_by_parental_age) + 
  geom_histogram(mapping = aes(x = n_paternal_dnm, alpha = 0.5)) +
  ylab("Frequency")+
  xlab("Number of Paternal de Novo Mutations")

ggplot(data=dnm_by_parental_age) + 
  geom_histogram(mapping = aes(x = n_maternal_dnm, alpha = 0.5)) +
  ylab("Frequency")+
  xlab("Number of Maternal de Novo Mutations")


# Step 2.6
t.test(dnm_by_parental_age$n_paternal_dnm, dnm_by_parental_age$n_maternal_dnm, paired = TRUE)
# I chose a paired t-test because I am trying to find the significant difference between two measurements. In this case, the difference between maternal and paternal mutations.
# My test was statistically significant. This means that there is a difference between maternal and faternal mutation rates.
