# get prior for effect of bacterial infection
library(tidyverse)
library(lme4)
library(broom.mixed)

# GET DATA-----------------

# set working directory
setwd(dirname(file.choose()))

# DATA FROM STUDY 1 
# get unpublished raw data on allogrooming and self-grooming from non-fasted bats (in a small cage)
# Stockmaier S, Bolnick DI, Page RA, Carter GG. 2018. An immune challenge reduces social grooming in vampire bats. Animal Behaviour 140: 141-149.
# raw data download from dropbox here:
# https://www.dropbox.com/scl/fi/bf4wujtec4luyiicjmi7g/Stockmaier_2018_behavior_data.csv?rlkey=v2xxpj3fqhw2y04swnhuouz8a&dl=0

# get self-grooming data
s <- read.csv("Stockmaier_2018_behavior_data.csv")

# DATA FROM STUDY 2
# get published allogrooming data from fasted bats injected with LPS (in a flight cage)
# Stockmaier, Sebastian (2019). Dataset and R.code: "Sickness effects on social interactions depend on the type of behaviour and relationship". figshare. Dataset. 
# download from here:
# https://doi.org/10.6084/m9.figshare.7726256.v7

# get allogrooming data
a <- read.csv('rates_LPS_trials_treatment.csv')

# ANALYZE SOCIAL GROOMING DATA ---------------------------------

# tidy allogrooming data
a2 <- 
  a %>% 
  filter(behav=="g") %>% 
  group_by(actor, treatment) %>% 
  summarize(grooming= sum(rate)) %>% 
  ungroup() %>% 
  mutate(total= 60*60) 

# fit binomial GLMM
fit <- glmer(cbind(grooming, total-grooming) ~ treatment + (1|actor),
             data= a2,
             family= binomial)

# get mean and the standard deviation of the mean (log odds)
summary(fit)
median <- summary(fit)$coefficients[2,1]
median_SD <- summary(fit)$coefficients[2,2]

# mean and SD of effect of infection on allogrooming by fasted bats in flight cage (Stockmaier et al. 2020 J Anim Ecol)
median # 0.634 log odds
median_SD # 0.012 log odds

# these are prior estimates of the model coefficient

# get odds ratio
OR <- tidy(fit,conf.int=TRUE,exponentiate=T,effects="fixed", conf.method="profile")
OR
# Odds Ratio = 1.89, 95% CI = 1.84, 1.93
# LPS injection changes the odds of self-grooming by a factor of 1.89 for fasted bats in a flight cage (Stockmaier et al. 2020 J Anim Ecol)


# now get estimates from non-fasted bats in small cages

# fit binomial GLMM
fit1 <- glmer(cbind(socialgrooming.count, obs.count-socialgrooming.count) ~ treatment + (1|focal.bat),
              data= s,
              family= binomial)

# get mean and the standard deviation of the mean
summary(fit1)
median1 <- summary(fit1)$coefficients[2,1]
median_SD1 <- summary(fit1)$coefficients[2,2]

# mean and its SD of effect of infection on allogrooming given in non-fasted bats in small cage (Stockmaier et al. 2018 Anim Behav) 
median1 # 1.316
median_SD1 # 0.269

# get odds ratio
OR1 <- tidy(fit1,conf.int=TRUE,exponentiate=T,effects="fixed", conf.method="profile")
OR1
# Odds Ratio = 3.73, 95% CI = 2.25, 6.50
# LPS injection changes the odds of self-grooming by a factor of 3.73 for non-fasted bats in a small cage (Stockmaier et al. 2018 Anim Behav) 

# take average of two effects (log odds)
mean(c(median, median1)) # 0.975 
mean(c(median_SD, median_SD1)) # 0.140

# NOW ANALYZE SELF_GROOMING---------------
     
# fit binomial GLMM
fit2 <- glmer(cbind(grooming.count, obs.count-grooming.count) ~ treatment + (1|focal.bat),
             data= s,
             family= binomial)

# get mean and the standard deviation of the mean
summary(fit2)
mean2 <- summary(fit2)$coefficients[2,1]
mean_SD2 <- summary(fit2)$coefficients[2,2]

# get odds ratio
OR2 <- tidy(fit2,conf.int=TRUE,exponentiate=T,effects="fixed", conf.method="profile")
OR2
# Odds Ratio = 16.1, 95% CI = 11, 24.7

# mean and its SD of effect of infection self-grooming in non-fasted bats in small cage (Stockmaier et al. 2018 Anim Behav)  
mean2  # 2.78
mean_SD2 # 0.204


# To estimate a prior probability distribution for the effect of an infection on allogrooming, we took the estimated effects from data from two controlled experiments measuring the effect of a simulated bacterial infection on allogrooming in vampire bats. To estimate odds ratios for the effects of simulated infections on allogrooming, we fit binomial mixed effects models to the raw data (cite supplementary data). In the first study, vampire bats were injected with either saline or lipopolysaccharide (LPS) and observed in small groups within small cages (Stockmaier et al. 2018 Anim Behav). Under these conditions, LPS injection reduced the odds of allogrooming by a factor of 3.73 (Odds Ratio, 95% CI = 2.25, 6.50) . In the second study, vampire bats were fasted to induce food sharing then injected with either saline or LPS and observed in larger groups in a flight cage (Stockmaier et al. 2018 Anim Behav). Under these conditions, LPS injection reduced the odds of allogrooming by a factor of 1.89 (Odds Ratio, 95% CI = 1.84, 1.93). This difference in effect size could be due to the fasting of the bats or due to the differences in group and cage size. In the present study, bats were unfasted like the first study but held in a large group in a flight cage, like the second study. As our prior, we therefore averaged the two means (log odds = 0.975) and averaged their corresponding standard errors (0.140). To estimate a prior probability distribution for the effect of an infection on self-grooming and on movement, we used data from the first study (Stockmaier et al. 2018 Anim Behav), where LPS injection reduced the odds of self-grooming by a factor of 16 (Odds Ratio, 95% CI = 11, 25; log odds = 2.78, SE = 0.204). Self-grooming and movement were not measured in the second study (Stockmaier et al. 2018 Anim Behav).
