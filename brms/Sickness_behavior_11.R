# Compare individual behavior in bats with and without Staph infection

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(boot)

# set working directory manually
setwd(dirname(file.choose())) 

############### clean data ####################
# get raw data 
raw <- read.csv("Data_sickness_behavior_02.csv")

# look at data
head(raw)
str(raw)

# start data at row 6
d <- raw[6:nrow(raw),1:13]  
head(d)
colnames(d)
# remove raw data
rm(raw)

# relabel column names
colnames(d) <- c("obs", "observer", "file", "camera", "year", "month", "day", "video.time", "event.time", "bat", "infected", "behav", "comment")
head(d)

# check observations
d$obs <- as.numeric(d$obs)
unique(d$observer)
unique(d$file)
d$file[which(d$file=="na")] <- NA
d$file[which(d$file=="")] <- NA
unique(d$camera)
d$camera[which(d$camera=="na")] <- NA
d$camera <- paste0("camera", d$camera)
unique(d$year)
d$month[which(d$month== "7")] <- "07"
unique(d$month)
unique(d$day)
d$date <- paste(d$year, d$month, d$day, sep= "-")
d$date2 <- as.Date(d$date)
unique(d$date2)
unique(d$video.time)
unique(d$event.time)
unique(d$bat)
d$bat[which(d$bat=="")] <- NA
unique(d$infected)
d$infected <- as.logical(d$infected) 

# fix behaviors
unique(d$behav)
d$behav[which(d$behav== " sg")] <- "sg" 
d$behav[which(d$behav== "")] <- NA 
d$behav[which(d$behav== "m ")] <- "m"
d$behav[which(d$behav== " m")] <- "m"

# label offscreen as missing
d$behav[which(d$behav== "n")] <- "off-camera" 
d$behav[which(d$behav== "na")] <- "off-camera"

# label treatment period
d$period = ifelse(d$date > "2019-08-05", "post-treatment", "pre-treatment")

# get infected bats
infected.bats <- 
  d %>% 
  filter(infected) %>% 
  pull(bat) %>% 
  unique()

# get injured bats
injured.bats <- c("cnone", 'cldw', 'cfd', 'cww')

# get all combinations of bats and behaviors and set to zero
zeros <- 
  expand_grid(date2= unique(d$date2), bat= unique(d$bat), behav= unique(d$behav), period= unique(d$period)) %>% 
  mutate(infected = bat %in% infected.bats) %>% 
  filter(!is.na(behav)) %>% 
  mutate(behav= case_when(
    behav == 'g' ~ "self-grooming",
    behav == 's' ~ "resting", # combine sleeping and resting
    behav == 'r' ~ "resting",
    behav == 'sg' ~ "allogrooming",
    behav == 'sh' ~ "allogrooming", # include mouthlicking as allogrooming
    behav == 'm' ~ "moving",
    behav == "off-camera" ~ NA_character_,
    TRUE ~ behav))  %>% 
  group_by(bat, behav, period, date2) %>% 
  filter(!is.na(behav)) %>% 
  mutate(n= 0) 

# summarize data into proportions (this assumes equal sampling across bats)
d2 <- 
  d %>% 
  # relabel variables
  mutate(behav= case_when(
    behav == 'g' ~ "self-grooming",
    behav == 's' ~ "resting", # combine sleeping and resting
    behav == 'r' ~ "resting",
    behav == 'sg' ~ "allogrooming",
    behav == 'sh' ~ "allogrooming", # combine mouthlicking and allogrooming
    behav == 'm' ~ "moving",
    behav == "off-camera" ~ NA_character_,
    TRUE ~ behav)) %>%
  # get on camera time
  filter(!is.na(behav)) %>% 
  # get counts per day
  group_by(infected, bat, period, behav, date2) %>% 
  summarize(n=n()) %>% 
  full_join(zeros) %>% 
  # get per bat
  group_by(infected, bat, period, behav) %>% 
  summarize(n= sum(n, na.rm=T)) %>%
  # get total sampling time on camera
  group_by(infected, bat, period) %>% 
  mutate(total= sum(n, na.rm=T)) %>% 
  mutate(prop = n/total) %>% 
  ungroup() %>% 
  mutate(period = factor(period, levels= c("pre-treatment", "post-treatment"))) %>% 
  mutate(injured = bat %in% injured.bats) %>% 
  mutate(infected2 = case_when(
    injured ~ "injured",
    infected ~ "treated",
    !infected ~ "asymptomatic")) %>% 
  # ignore resting for simplicity
  filter(behav!= 'resting') %>% 
  filter(total>0)

#################### get bootstrapped 95% CIs #######################
# create functions for bootstrapping for 95% CI
# get mean and 95% CI of values x via bootstrapping
boot_ci <- function(x, perms=5000, bca=F) {
  get_mean <- function(x, d) {
    return(mean(x[d]))
  }
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  if(bca){
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="bca")
    low <- boot$bca[1,4]
    high <- boot$bca[1,5] 
  }else{
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="perc")
    low <- boot$perc[1,4]
    high <- boot$perc[1,5] 
  }
  c(low=low,mean=mean,high=high, N=round(length(x)))
}


# get mean and 95% CI via bootstrapping of values y within grouping variable x
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000, bca=F){
  df <- data.frame(effect=unique(x))
  df$low <- NA
  df$mean <- NA
  df$high <- NA
  df$n.obs <- NA
  for (i in 1:nrow(df)) {
    ys <- y[which(x==df$effect[i])]
    if (length(ys)>1 & var(ys)>0 ){
      b <- boot_ci(y[which(x==df$effect[i])], perms=perms, bca=bca) 
      df$low[i] <- b[1]
      df$mean[i] <- b[2]
      df$high[i] <- b[3]
      df$n.obs[i] <- b[4]
    }else{
      df$low[i] <- min(ys)
      df$mean[i] <- mean(ys)
      df$high[i] <- max(ys)
      df$n.obs[i] <- length(ys)
    }
  }
  df
}

# get means and 95%CI for infected vs not
means <- 
  d2 %>% 
  mutate(temp= paste(infected, period, behav, sep= "_")) %>% 
  select(temp, prop) %>% 
  boot_ci2(x= .$temp, y= .$prop) %>% 
  separate(effect, into= c("infected", "period", "behav"), sep= "_") %>% 
  mutate(period = factor(period, levels= c("pre-treatment", "post-treatment")))

# make table of means and 95% CIs
table <- 
  means %>% 
  mutate(infected= ifelse(infected, "infected", 'asymptomatic')) %>% 
  select(infected, period, behav, low, mean, high) %>% 
  pivot_longer(low:high) %>% 
  mutate(value= signif(value,2)) %>% 
  mutate(label= paste(period, behav)) %>% 
  mutate(name= paste(infected, name)) %>% 
  select(-period, -behav, -infected) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  separate(label, into=c("period", "behavior"), sep= " ") %>% 
  arrange(desc(period)) %>% 
  mutate(significant_results = case_when(
    `asymptomatic low` > `infected high` ~ "infected is lower",
    `asymptomatic high` < `infected low` ~ "infected is higher",
    TRUE ~ "no clear difference"))

write.csv(table, "staph_behav.csv")

# get means and 95% CI for injured vs not
means2 <- 
  d2 %>% 
  mutate(temp= paste(infected2, period, behav, sep= "_")) %>% 
  select(temp, prop) %>% 
  boot_ci2(x= .$temp, y= .$prop) %>% 
  separate(effect, into= c("infected2", "period", "behav"), sep= "_") %>% 
  mutate(period = factor(period, levels= c("pre-treatment", "post-treatment")))

# make table of means and 95% CI
table2 <- 
  means2 %>% 
  select(infected2, period, behav, low, mean, high) %>% 
  pivot_longer(low:high) %>% 
  mutate(value= signif(value,2)) %>% 
  mutate(label= paste(period, behav)) %>% 
  mutate(name= paste(infected2, name)) %>% 
  select(-period, -behav, -infected2) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  separate(label, into=c("period", "behavior"), sep= " ") %>% 
  arrange(desc(period)) %>% 
  mutate(significant= case_when(
    `asymptomatic low` > `treated high` ~ "treated is lower",
    `asymptomatic high` < `treated low` ~ "treated is higher",
    `asymptomatic low` > `injured high` ~ "injured is lower",
    `asymptomatic high` < `injured low` ~ "injured is higher",
    TRUE ~ "no clear difference"))

write.csv(table2, "injured_behav.csv")

# FIT MODELS #######################################################

# load packages
library(brms)
library(tidybayes)
library(bayesplot)
library(patchwork)

# choose whether to analyze infected or injured bats 
INFECTED <- TRUE
INJURED <- FALSE

# INFECTED #########################################################

if (INFECTED == TRUE){
  
  # ALLOGROOMING-----
  # interaction -----
  a <- d2 %>% filter(behav=="allogrooming")  
  
  # model fit
  fita <- brm(n | trials(total) ~ infected*period + (1|bat),
              data = a,
              family = binomial(),
              chains = 4,
              #prior = prior(normal(-0.975, 0.140), class = b, coef = infectedTRUE), # set prior for infection
              iter = 5000, warmup = 1000)
  statsa <- as_tibble(summary(fita, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "interaction")
  statsa <- 
    statsa %>% 
    mutate(effect = row.names(summary(fita)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # there is a significant interaction between infection and period
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita)
  pp_check(fita)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa <- as.matrix(fita)
  
  dimnames(posterior_samplesa)$variable[1] <- "intercept"
  dimnames(posterior_samplesa)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesa)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_samplesa)$variable[4] <- "infectedTRUE:periodPOST-TREATMENT"
  dimnames(posterior_samplesa)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t <- mcmc_intervals(posterior_samplesa, 
                 regex_pars = c("intercept", "infectedTRUE", "periodPOST", "infectedTRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(c(-4, 2))+
    labs(title = "Effects on allogrooming"))
    
  # Calculate median and HDI for each combination of 'infected' and 'period'
  pointsa <- 
    a %>%
    add_predicted_draws(fita) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(infected, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fita <- a %>%
      add_predicted_draws(fita) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = infected, color = infected))+  
      geom_jitter(data = a, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointsa, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointsa, aes(x = period, 
                                          y = median, 
                                          ymin = l_hdi, 
                                          ymax = u_hdi,
                                          width = 0.1), 
                      size = 0.75, 
                      position = position_dodge(width = 0.75))+
      geom_line(data = pointsa, aes(y = median, group = infected), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("") +
      theme_bw() +
      theme(legend.position = "none"))
  
  # # plot effects for each bat
  # fita %>%
  #   spread_rvars(b_Intercept, r_bat[bat,]) %>%
  #   mutate(condition_mean = b_Intercept + r_bat) %>%
  #   ggplot(aes(y = bat, xdist = condition_mean)) +
  #   stat_halfeye()
  
  # pre-treatment -----
  a1 <- d2 %>% filter(behav=="allogrooming") %>%  filter(period=="pre-treatment") 
  
  # model fit
  fita1 <- brm(n | trials(total) ~ infected + (1|bat),
              data = a1,
              family = binomial(),
              chains = 4,
              #prior = prior(normal(-0.975, 0.140), class = b, coef = infectedTRUE), # set prior for infection
              iter = 5000, warmup = 1000)
  statsa1 <- as_tibble(summary(fita1, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "pre-treatment")
  statsa1 <- 
    statsa1 %>% 
    mutate(effect = row.names(summary(fita1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # significant negative effect of infection
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita1)
  pp_check(fita1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa1 <- as.matrix(fita1)
  
  dimnames(posterior_samplesa1)$variable[1] <- "intercept"
  dimnames(posterior_samplesa1)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesa1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t2 <- mcmc_intervals(posterior_samplesa1, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(c(-4, 2))+
    labs(title = "Effects on allogrooming (pre-treatment)"))
  
  # post-treatment -----
  a2 <- d2 %>% filter(behav=="allogrooming") %>%  filter(period=="post-treatment")
  
  # model fit
  fita2 <- brm(n | trials(total) ~ infected + (1|bat),
               data = a2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statsa2 <- as_tibble(summary(fita2, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "post-treatment")
  statsa2 <- 
    statsa2 %>% 
    mutate(effect = row.names(summary(fita2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # no significant effect of infection 
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita2)
  pp_check(fita2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa2 <- as.matrix(fita2)
  
  dimnames(posterior_samplesa2)$variable[1] <- "intercept"
  dimnames(posterior_samplesa2)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesa2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t3 <- mcmc_intervals(posterior_samplesa2, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    xlim(c(-4, 2))+
    labs(title = "Effects on allogrooming (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, infected pre-treatment, and infected post-treatment
  (effectsa <- 
    rbind(as_tibble(summary(fita, prob = 0.90)$fixed)[4,],
        as_tibble(summary(fita1, prob = 0.90)$fixed)[2,],
        as_tibble(summary(fita2, prob = 0.90)$fixed)[2,]) %>% 
      mutate(effect = c("interaction term (infected x treatment period)", "infected effect (pre-treatment)", "infected effect (post-treatment)")) %>% 
    rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
    ggplot() +
    geom_hline(yintercept = -1.316, linetype="dashed", 
                 color = "purple", size=0.75)+
    geom_hline(yintercept = -0.634, linetype="dashed", 
               color = "orange", size=0.75)+
    geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
    hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")+
    xlab("")+
    ylab("") +
    #ylim(c(-2.5, 2))+
    ylim(c(-3, 3)) +
    ggtitle("allogrooming"))
  
  # plot all effects for allogrooming -----
  plot(t / t2 / t3)
  
  #ggsave("allogrooming_effects_priors.PDF", width = 7, height = 5)
  
  # SELF - GROOMING -----
  
  # interaction -----
  s <- d2 %>% filter(behav=="self-grooming")   
  
  # model fit
  fits <- brm(n | trials(total) ~ infected*period + (1|bat),
              data = s,
              family = binomial(),
              chains = 4,
              #prior = prior(normal(-2.78, 0.204), class = b, coef = infectedTRUE), # set prior for infection
              iter = 5000, warmup = 1000)
  statss <- as_tibble(summary(fits, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "interaction")
  statss <- 
    statss %>% 
    mutate(effect = row.names(summary(fits)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # there is a significant interaction between infection and period (lower-95%CI close to zero, though)
  # no clear effect of infection
  # negative effect of period
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fits)
  pp_check(fits)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless <- as.matrix(fits)
  
  dimnames(posterior_sampless)$variable[1] <- "intercept"
  dimnames(posterior_sampless)$variable[2] <- "infectedTRUE"
  dimnames(posterior_sampless)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_sampless)$variable[4] <- "infectedTRUE:periodPOST-TREATMENT"
  dimnames(posterior_sampless)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t4 <- mcmc_intervals(posterior_sampless, 
                 regex_pars = c("intercept", "infectedTRUE", "periodPOST", "infectedTRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(-2,1)+
    labs(title = "Effects on self-grooming"))
  
  # Calculate median and HDI for each combination of 'infected' and 'period'
  pointss <- 
    s %>%
    add_predicted_draws(fits) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(infected, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fits <- s %>%
      add_predicted_draws(fits) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = infected, color = infected))+  
      geom_jitter(data = s, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointss, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointss, aes(x = period, 
                                        y = median, 
                                        ymin = l_hdi, 
                                        ymax = u_hdi,
                                        width = 0.1), 
                    size = 0.75, 
                    position = position_dodge(width = 0.75))+
      geom_line(data = pointss, aes(y = median, group = infected), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("") +
      theme(legend.title = element_blank())+
      theme_bw())
  
  # pre-treatment -----
  s1 <- d2 %>% filter(behav=="self-grooming") %>% filter(period == "pre-treatment")
  
  # model fit
  fits1 <- brm(n | trials(total) ~ infected + (1|bat),
              data = s1,
              family = binomial(),
              chains = 4,
              #prior = prior(normal(-2.78, 0.204), class = b, coef = infectedTRUE), # set prior for infection
              iter = 5000, warmup = 1000)
  statss1 <- as_tibble(summary(fits1, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "pre-treatment")
  statss1 <- 
    statss1 %>% 
    mutate(effect = row.names(summary(fits1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # no significant effect of infection 
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fits1)
  pp_check(fits1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless1 <- as.matrix(fits1)
  
  dimnames(posterior_sampless1)$variable[1] <- "intercept"
  dimnames(posterior_sampless1)$variable[2] <- "infectedTRUE"
  dimnames(posterior_sampless1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t5 <- mcmc_intervals(posterior_sampless1, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(-2, 1)+
    labs(title = "Effects on self-grooming (pre-treatment)"))
  
  # post-treatment -----
  s2 <- d2 %>% filter(behav=="self-grooming") %>% filter(period == "post-treatment")
  
  # model fit
  fits2 <- brm(n | trials(total) ~ infected + (1|bat),
               data = s2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statss2 <- as_tibble(summary(fits2, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "post-treatment")
  statss2 <- 
    statss2 %>% 
    mutate(effect = row.names(summary(fits2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # no significant effect of infection, but there is a trend towards a positive effect
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fits2)
  pp_check(fits2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless2 <- as.matrix(fits2)
  
  dimnames(posterior_sampless2)$variable[1] <- "intercept"
  dimnames(posterior_sampless2)$variable[2] <- "infectedTRUE"
  dimnames(posterior_sampless2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t6 <- mcmc_intervals(posterior_sampless2, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(-2,1)+
    labs(title = "Effects on self-grooming (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, infected pre-treatment, and infected post-treatment
  (effectss <- 
      rbind(as_tibble(summary(fits, prob = 0.90)$fixed)[4,],
        as_tibble(summary(fits1, prob = 0.90)$fixed)[2,],
        as_tibble(summary(fits2, prob = 0.90)$fixed)[2,]) %>% 
    mutate(effect = c("interaction term (infected x treatment period)", "infected effect (pre-treatment)", "infected effect (post-treatment)")) %>% 
    rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
    ggplot() +
    geom_hline(yintercept = -2.78, linetype="dashed", 
                 color = "purple", size=0.75)+
    geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
    hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")+
    xlab("")+
    ylab("") +
    #ylim(c(-2.5, 2))+
    ylim(c(-3, 3)) +
    ggtitle("self-grooming"))
  
  # plot all effects for selfgrooming -----
  plot(t4 / t5 / t6)
  
  #ggsave("selfgrooming_effects_priors.PDF", width = 7, height = 5)
  
  # MOVING -----
  
  # interaction -----
  m <- d2 %>% filter(behav=="moving")   
  
  # model fit
  fitm <- brm(n | trials(total) ~ infected*period + (1|bat),
              data = m,
              family = binomial(),
              chains = 4,
              #prior = prior(normal(-2.78, 0.204), class = b, coef = infectedTRUE), # set prior for infection
              iter = 5000, warmup = 1000)
  statsm <- as_tibble(summary(fitm, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "interaction")
  statsm <- 
    statsm %>% 
    mutate(effect = row.names(summary(fitm)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # significant interaction between infection and period
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm)
  pp_check(fitm)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm <- as.matrix(fitm)
  
  dimnames(posterior_samplesm)$variable[1] <- "intercept"
  dimnames(posterior_samplesm)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesm)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_samplesm)$variable[4] <- "infectedTRUE:periodPOST-TREATMENT"
  dimnames(posterior_samplesm)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t7 <- mcmc_intervals(posterior_samplesm, 
                 regex_pars = c("intercept", "infectedTRUE", "periodPOST", "infectedTRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(-2,1.5)+
    labs(title = "Effects on moving"))
  
  # Calculate median and HDI for each combination of 'infected' and 'period'
  pointsm <- 
    m %>%
    add_predicted_draws(fitm) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(infected, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fitm <- m %>%
      add_predicted_draws(fitm) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = infected, color = infected))+  
      geom_jitter(data = m, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointsm, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointsm, aes(x = period, 
                                        y = median, 
                                        ymin = l_hdi, 
                                        ymax = u_hdi,
                                        width = 0.1), 
                    size = 0.75, 
                    position = position_dodge(width = 0.75))+
      geom_line(data = pointsm, aes(y = median, group = infected), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("period") +
      theme_bw() +
      theme(legend.position = "none"))
  
  # pre-treatment -----
  m1 <- d2 %>% filter(behav=="moving") %>% filter(period == "pre-treatment")  
  
  # model fit
  fitm1 <- brm(n | trials(total) ~ 0 + infected + (1|bat),
              data = m1,
              family = binomial(),
              #prior = prior(normal(-2.78, 0.204), class = b, coef = infectedTRUE), # set prior for infection
              chains = 4,
              iter = 5000, warmup = 1000)
  statsm1 <- as_tibble(summary(fitm1, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "pre-treatment")
  statsm1 <- 
    statsm1 %>% 
    mutate(effect = row.names(summary(fitm1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # significant negative effect of infection
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm1)
  pp_check(fitm1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm1 <- as.matrix(fitm1)
  
  dimnames(posterior_samplesm1)$variable[1] <- "intercept"
  dimnames(posterior_samplesm1)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesm1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t8 <- mcmc_intervals(posterior_samplesm1, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    #xlim(-2,1.5)+
    labs(title = "Effects on moving (pre-treatment)"))
  
  # post-treatment -----
  m2 <- d2 %>% filter(behav=="moving") %>% filter(period == "post-treatment")  
  
  # model fit
  fitm2 <- brm(n | trials(total) ~ infected + (1|bat),
               data = m2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statsm2 <- as_tibble(summary(fitm2, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "post-treatment")
  statsm2 <- 
    statsm2 %>% 
    mutate(effect = row.names(summary(fitm2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  # no significant effect of infection
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm2)
  pp_check(fitm2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm2 <- as.matrix(fitm2)
  
  dimnames(posterior_samplesm2)$variable[1] <- "intercept"
  dimnames(posterior_samplesm2)$variable[2] <- "infectedTRUE"
  dimnames(posterior_samplesm2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t9 <- mcmc_intervals(posterior_samplesm2, 
                 regex_pars = c("intercept", "infectedTRUE", "SD bat intercept"), 
                 prob = 0.90) +
    theme_bw()+
    xlim(-2,1.5)+
    labs(title = "Effects on moving (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, infected pre-treatment, and infected post-treatment
  (effectsm <- 
      rbind(as_tibble(summary(fitm, prob = 0.90)$fixed)[4,],
        as_tibble(summary(fitm1, prob = 0.90)$fixed)[2,],
        as_tibble(summary(fitm2, prob = 0.90)$fixed)[2,]) %>% 
    mutate(effect = c("interaction term (infected x treatment period)", "infected effect (pre-treatment)", "infected effect (post-treatment)")) %>% 
    rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
    ggplot() +
    geom_hline(yintercept = -2.78, linetype="dashed", 
                 color = "purple", size=0.75)+
    geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
    hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")+
    xlab("")+
    ylab("posterior probability") +
    #ylim(c(-2.5, 2)) +
    ylim(c(-3, 3)) +
    ggtitle("moving"))
  
  # plot all effects for moving ----
  plot(t7 / t8 / t9)
  
  #ggsave("moving_effects_priors.PDF", width = 7, height = 5)
  
  # plot figure -----
  p1 <- plot(effectsa + model_fita)
  p2 <- plot(effectss + model_fits)
  p3 <- plot(effectsm + model_fitm)
  
  plot(p1/ p2/ p3) + plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "")))
  
  ggsave("Fig1_infected.pdf", width = 9, height = 7)
  
  # save results -----
  brms_results <- rbind(statsa, statsa1, statsa2, statss, statss1, statss2, statsm, statsm1, statsm2)
  
  write.csv(brms_results, "brms_results_infected.csv")
}else{NULL}
  
# INJURED ###############################################################

if (INJURED == TRUE){
  
  # ALLOGROOMING -----
  # interaction
  fita <- brm(n | trials(total) ~ injured*period + (1|bat),
              data = a,
              family = binomial(),
              #prior = prior(normal(-0.975, 0.140), class = b, coef = injuredTRUE), # set prior for infection
              chains = 4,
              iter = 5000, warmup = 1000)
  statsa <- as_tibble(summary(fita, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "interaction")
  statsa <- 
    statsa %>% 
    mutate(effect = row.names(summary(fita)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita)
  pp_check(fita)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa <- as.matrix(fita)
  
  dimnames(posterior_samplesa)$variable[1] <- "intercept"
  dimnames(posterior_samplesa)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesa)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_samplesa)$variable[4] <- "injuredTRUE:periodPOST-TREATMENT"
  dimnames(posterior_samplesa)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t <- mcmc_intervals(posterior_samplesa, 
                       regex_pars = c("intercept", "injuredTRUE", "periodPOST", "TRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                       prob = 0.90) +
      theme_bw()+
      #xlim(c(-4, 2))+
      labs(title = "Effects on allogrooming"))
  
  # Calculate median and HDI for each combination of 'injured' and 'period'
  pointsa <- 
    a %>%
    add_predicted_draws(fita) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(injured, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fita <- a %>%
      add_predicted_draws(fita) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = injured, color = injured))+  
      geom_jitter(data = a, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointsa, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointsa, aes(x = period, 
                                        y = median, 
                                        ymin = l_hdi, 
                                        ymax = u_hdi,
                                        width = 0.1), 
                    size = 0.75, 
                    position = position_dodge(width = 0.75))+
      geom_line(data = pointsa, aes(y = median, group = injured), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("") +
      theme_bw() +
      theme(legend.position = "none"))
  
  # # plot effects for each bat
  # fita %>%
  #   spread_rvars(b_Intercept, r_bat[bat,]) %>%
  #   mutate(condition_mean = b_Intercept + r_bat) %>%
  #   ggplot(aes(y = bat, xdist = condition_mean)) +
  #   stat_halfeye()
  
  # pre-treatment 
  a1 <- d2 %>% filter(behav=="allogrooming") %>%  filter(period=="pre-treatment") 
  
  # model fit
  fita1 <- brm(n | trials(total) ~ injured + (1|bat),
               data = a1,
               family = binomial(),
               #prior = prior(normal(-0.975, 0.140), class = b, coef = injuredTRUE), # set prior for infection
               chains = 4,
               iter = 5000, warmup = 1000)
  statsa1 <- as_tibble(summary(fita1, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "pre-treatment")
  statsa1 <- 
    statsa1 %>% 
    mutate(effect = row.names(summary(fita1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita1)
  pp_check(fita1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa1 <- as.matrix(fita1)
  
  dimnames(posterior_samplesa1)$variable[1] <- "intercept"
  dimnames(posterior_samplesa1)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesa1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t2 <- mcmc_intervals(posterior_samplesa1, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(c(-4, 2))+
      labs(title = "Effects on allogrooming (pre-treatment)"))
  
  # post-treatment 
  a2 <- d2 %>% filter(behav=="allogrooming") %>%  filter(period=="post-treatment")
  
  # model fit
  fita2 <- brm(n | trials(total) ~ injured + (1|bat),
               data = a2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statsa2 <- as_tibble(summary(fita2, prob = 0.90)$fixed) %>% mutate(behavior = "grooming") %>% mutate(model = "post-treatment")
  statsa2 <- 
    statsa2 %>% 
    mutate(effect = row.names(summary(fita2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fita2)
  pp_check(fita2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesa2 <- as.matrix(fita2)
  
  dimnames(posterior_samplesa2)$variable[1] <- "intercept"
  dimnames(posterior_samplesa2)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesa2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t3 <- mcmc_intervals(posterior_samplesa2, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(c(-4, 2))+
      labs(title = "Effects on allogrooming (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, injured pre-treatment, and injured post-treatment
  (effectsa <- 
      rbind(as_tibble(summary(fita, prob = 0.90)$fixed)[4,],
            as_tibble(summary(fita1, prob = 0.90)$fixed)[2,],
            as_tibble(summary(fita2, prob = 0.90)$fixed)[2,]) %>% 
      mutate(effect = c("interaction term (injured x treatment period)", "injured effect (pre-treatment)", "injured effect (post-treatment)")) %>% 
      rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
      ggplot() +
      geom_hline(yintercept = -1.316, linetype="dashed", 
                 color = "purple", size=0.75)+
      geom_hline(yintercept = -0.634, linetype="dashed", 
                 color = "orange", size=0.75)+
      geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
      hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none")+
      xlab("")+
      ylab("") +
      #ylim(c(-2.5, 2))+
      ylim(c(-3, 3)) +
      ggtitle("allogrooming"))
  
  # plot all effects for allogrooming 
  plot(t / t2 / t3)
  
  #ggsave("allogrooming_effects.PDF", width = 7, height = 5)
  
  # SELF - GROOMING -----
  
  # interaction
  s <- d2 %>% filter(behav=="self-grooming")   
  
  # model fit
  fits <- brm(n | trials(total) ~ injured*period + (1|bat),
              data = s,
              family = binomial(),
              #prior = prior(normal(-2.78, 0.204), class = b, coef = injuredTRUE), # set prior for infection
              chains = 4,
              iter = 5000, warmup = 1000)
  statss <- as_tibble(summary(fits, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "interaction")
  statss <- 
    statss %>% 
    mutate(effect = row.names(summary(fits)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fits)
  pp_check(fits)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless <- as.matrix(fits)
  
  dimnames(posterior_sampless)$variable[1] <- "intercept"
  dimnames(posterior_sampless)$variable[2] <- "injuredTRUE"
  dimnames(posterior_sampless)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_sampless)$variable[4] <- "injuredTRUE:periodPOST-TREATMENT"
  dimnames(posterior_sampless)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t4 <- mcmc_intervals(posterior_sampless, 
                        regex_pars = c("intercept", "injuredTRUE", "periodPOST", "injuredTRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2,1)+
      labs(title = "Effects on self-grooming"))
  
  # Calculate median and HDI for each combination of 'injured' and 'period'
  pointss <- 
    s %>%
    add_predicted_draws(fits) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(injured, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fits <- s %>%
      add_predicted_draws(fits) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = injured, color = injured))+  
      geom_jitter(data = s, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointss, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointss, aes(x = period, 
                                        y = median, 
                                        ymin = l_hdi, 
                                        ymax = u_hdi,
                                        width = 0.1), 
                    size = 0.75, 
                    position = position_dodge(width = 0.75))+
      geom_line(data = pointss, aes(y = median, group = injured), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("") +
      theme(legend.title = element_blank())+
      theme_bw())
  
  # pre-treatment
  s1 <- d2 %>% filter(behav=="self-grooming") %>% filter(period == "pre-treatment")
  
  # model fit
  fits1 <- brm(n | trials(total) ~ injured + (1|bat),
               data = s1,
               family = binomial(),
               #prior = prior(normal(-2.78, 0.204), class = b, coef = injuredTRUE), # set prior for infection
               chains = 4,
               iter = 5000, warmup = 1000)
  statss1 <- as_tibble(summary(fits1, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "pre-treatment")
  statss1 <- 
    statss1 %>% 
    mutate(effect = row.names(summary(fits1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fits1)
  pp_check(fits1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless1 <- as.matrix(fits1)
  
  dimnames(posterior_sampless1)$variable[1] <- "intercept"
  dimnames(posterior_sampless1)$variable[2] <- "injuredTRUE"
  dimnames(posterior_sampless1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t5 <- mcmc_intervals(posterior_sampless1, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2, 1)+
      labs(title = "Effects on self-grooming (pre-treatment)"))
  
  # post-treatment 
  s2 <- d2 %>% filter(behav=="self-grooming") %>% filter(period == "post-treatment")
  
  # model fit
  fits2 <- brm(n | trials(total) ~ injured + (1|bat),
               data = s2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statss2 <- as_tibble(summary(fits2, prob = 0.90)$fixed) %>% mutate(behavior = "self-grooming") %>% mutate(model = "post-treatment")
  statss2 <- 
    statss2 %>% 
    mutate(effect = row.names(summary(fits2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)

  # assess model fit and obtain posterior probability distribution plots
  plot(fits2)
  pp_check(fits2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_sampless2 <- as.matrix(fits2)
  
  dimnames(posterior_sampless2)$variable[1] <- "intercept"
  dimnames(posterior_sampless2)$variable[2] <- "injuredTRUE"
  dimnames(posterior_sampless2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t6 <- mcmc_intervals(posterior_sampless2, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2,1)+
      labs(title = "Effects on self-grooming (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, injured pre-treatment, and injured post-treatment
  (effectss <- 
      rbind(as_tibble(summary(fits, prob = 0.90)$fixed)[4,],
            as_tibble(summary(fits1, prob = 0.90)$fixed)[2,],
            as_tibble(summary(fits2, prob = 0.90)$fixed)[2,]) %>% 
      mutate(effect = c("interaction term (injured x treatment period)", "injured effect (pre-treatment)", "injured effect (post-treatment)")) %>% 
      rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
      ggplot() +
      geom_hline(yintercept = -2.78, linetype="dashed", 
                 color = "purple", size=0.75)+
      geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
      hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none")+
      xlab("")+
      ylab("") +
      #ylim(c(-2.5, 2))+
      ylim(c(-3, 3)) +
      ggtitle("self-grooming"))
  
  # plot all effects for selfgrooming
  plot(t4 / t5 / t6)
  
  #ggsave("selfgrooming_effects.PDF", width = 7, height = 5)
  
  # MOVING -----
  
  # interaction 
  m <- d2 %>% filter(behav=="moving")   
  
  # model fit
  fitm <- brm(n | trials(total) ~ 0 + injured*period + (1|bat),
              data = m,
              family = binomial(),
              #prior = prior(normal(-2.78, 0.204), class = b, coef = injuredTRUE), # set prior for infection
              chains = 4,
              iter = 5000, warmup = 1000)
  statsm <- as_tibble(summary(fitm, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "interaction")
  statsm <- 
    statsm %>% 
    mutate(effect = row.names(summary(fitm)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm)
  pp_check(fitm)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm <- as.matrix(fitm)
  
  dimnames(posterior_samplesm)$variable[1] <- "intercept"
  dimnames(posterior_samplesm)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesm)$variable[3] <- "periodPOST-TREATMENT"
  dimnames(posterior_samplesm)$variable[4] <- "injuredTRUE:periodPOST-TREATMENT"
  dimnames(posterior_samplesm)$variable[5] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t7 <- mcmc_intervals(posterior_samplesm, 
                        regex_pars = c("intercept", "injuredTRUE", "periodPOST", "injuredTRUE:periodPOST-TREATMENT", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2,1.5)+
      labs(title = "Effects on moving"))
  
  # Calculate median and HDI for each combination of 'injured' and 'period'
  pointsm <- 
    m %>%
    add_predicted_draws(fitm) %>% 
    mutate(.prediction = .prediction/total) %>% 
    group_by(injured, period) %>% 
    reframe(median = median(.prediction),
            l_hdi = quantile(.prediction, c(0.05)),
            u_hdi = quantile(.prediction, c(0.95)))
  
  # plot posterior probabilities pre- vs post-treatment
  (model_fitm <- m %>%
      add_predicted_draws(fitm) %>%  # adding the posterior distribution
      ggplot(aes(x = period, y = n/total, shape = injured, color = injured))+  
      geom_jitter(data = m, size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
      geom_point(data = pointsm, aes(x = period, y = median), size = 3, position = position_dodge(width = 0.75))+
      geom_errorbar(data = pointsm, aes(x = period, 
                                        y = median, 
                                        ymin = l_hdi, 
                                        ymax = u_hdi,
                                        width = 0.1), 
                    size = 0.75, 
                    position = position_dodge(width = 0.75))+
      geom_line(data = pointsm, aes(y = median, group = injured), position = position_dodge(width = 0.75))+
      scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
      scale_fill_brewer(palette = "Greys") +
      scale_color_manual(values=c("black", "red"))+
      ylab("proportion of time") +  
      xlab("period") +
      theme_bw() +
      theme(legend.position = "none"))
  
  # pre-treatment 
  m1 <- d2 %>% filter(behav=="moving") %>% filter(period == "pre-treatment")  
  
  # model fit
  fitm1 <- brm(n | trials(total) ~ 0 + injured + (1|bat),
               data = m1,
               family = binomial(),
               #prior = prior(normal(-2.78, 0.204), class = b, coef = injuredTRUE), # set prior for infection
               chains = 4,
               iter = 5000, warmup = 1000)
  statsm1 <- as_tibble(summary(fitm1, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "pre-treatment")
  statsm1 <- 
    statsm1 %>% 
    mutate(effect = row.names(summary(fitm1)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm1)
  pp_check(fitm1)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm1 <- as.matrix(fitm1)
  
  dimnames(posterior_samplesm1)$variable[1] <- "intercept"
  dimnames(posterior_samplesm1)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesm1)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t8 <- mcmc_intervals(posterior_samplesm1, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2,1.5)+
      labs(title = "Effects on moving (pre-treatment)"))
  
  # post-treatment 
  m2 <- d2 %>% filter(behav=="moving") %>% filter(period == "post-treatment")  
  
  # model fit
  fitm2 <- brm(n | trials(total) ~ injured + (1|bat),
               data = m2,
               family = binomial(),
               chains = 4,
               iter = 5000, warmup = 1000)
  statsm2 <- as_tibble(summary(fitm2, prob = 0.90)$fixed) %>% mutate(behavior = "moving") %>% mutate(model = "post-treatment")
  statsm2 <- 
    statsm2 %>% 
    mutate(effect = row.names(summary(fitm2)$fixed)) %>% 
    select(behavior, model, effect, Estimate:Tail_ESS)
  
  # assess model fit and obtain posterior probability distribution plots
  plot(fitm2)
  pp_check(fitm2)  # posterior predictive checks
  
  # plot posterior probability distributions for effects in a different format
  posterior_samplesm2 <- as.matrix(fitm2)
  
  dimnames(posterior_samplesm2)$variable[1] <- "intercept"
  dimnames(posterior_samplesm2)$variable[2] <- "injuredTRUE"
  dimnames(posterior_samplesm2)$variable[3] <- "SD bat intercept"
  
  # Plot posterior probability distributions for each effect
  color_scheme_set("gray")
  (t9 <- mcmc_intervals(posterior_samplesm2, 
                        regex_pars = c("intercept", "injuredTRUE", "SD bat intercept"), 
                        prob = 0.90) +
      theme_bw()+
      #xlim(-2,1.5)+
      labs(title = "Effects on moving (post-treatment)"))
  
  # Plot posterior probabilities for the interaction, injured pre-treatment, and injured post-treatment
  (effectsm <- 
      rbind(as_tibble(summary(fitm, prob = 0.90)$fixed)[4,],
            as_tibble(summary(fitm1, prob = 0.90)$fixed)[2,],
            as_tibble(summary(fitm2, prob = 0.90)$fixed)[2,]) %>% 
      mutate(effect = c("interaction term (injured x treatment period)", "injured effect (pre-treatment)", "injured effect (post-treatment)")) %>% 
      rename(lower = "l-90% CI", upper = "u-90% CI") %>% 
      ggplot() +
      geom_hline(yintercept = -2.78, linetype="dashed", 
                 color = "purple", size=0.75)+
      geom_pointrange(aes(x = effect, y = Estimate, ymin = lower, ymax = upper), linewidth = 1) +
      hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none")+
      xlab("")+
      ylab("posterior probability") +
      #ylim(c(-2.5, 2)) +
      ylim(c(-3, 3)) +
      ggtitle("moving"))
  
  # plot all effects for moving 
  plot(t7 / t8 / t9)
  
  #ggsave("moving_effects.PDF", width = 7, height = 5)
  
  # plot figure -----
  p1 <- plot(effectsa + model_fita)
  p2 <- plot(effectss + model_fits)
  p3 <- plot(effectsm + model_fitm)
  
  plot(p1/ p2/ p3) + plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "")))
  
  ggsave("Fig1_injured.pdf", width = 9, height = 7)
  
  # save results -----
  brms_results <- rbind(statsa, statsa1, statsa2, statss, statss1, statss2, statsm, statsm1, statsm2)
  
  write.csv(brms_results, "brms_results_injured.csv")
}else{NULL}




