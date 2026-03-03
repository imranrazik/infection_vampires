# Compare individual behavior in bats with and without Staph infection

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(boot)
library(emmeans)
library(broom)

# set working directory manually
setwd("~/Dropbox (Personal)/Dropbox/_working/_ACTIVE/staph/2026/Revised version + response to reviewers/scripts/") 

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

# FIT MODELS #######################################################

# load packages
library(brms)
library(tidybayes)
library(bayesplot)
library(patchwork)

# ALLOGROOMING-----
# interaction -----
a <- 
  d2 %>% 
  filter(behav=="allogrooming") %>% 
  mutate(post.treatment = period== "post-treatment")

# model fit
fita <- brm(n | trials(total) ~ infected*post.treatment + (1|bat),
            data = a,
            family = binomial(),
            chains = 4,
            iter = 5000, warmup = 1000)
statsa <- 
  as_tibble(summary(fita, prob = 0.90)$fixed) %>% 
  mutate(behavior = "grooming") %>% 
  mutate(model = "interaction") %>% 
  mutate(effect = row.names(summary(fita)$fixed)) %>% 
  select(behavior, model, effect, Estimate:Tail_ESS)

# interaction between infection and period
# interaction = 1.39    error= 0.139        lower 90%CI = 1.16      upper 90% CI = 1.62

# assess model fit and obtain posterior probability distribution plots
plot(fita)
pp_check(fita)  # posterior predictive checks

# Calculate mean and 95% credible intervals for each combination of 'infected' and 'period'

new_data <- 
  expand.grid(
    infected = c(TRUE, FALSE),
    post.treatment   = c(TRUE, FALSE),
    total= round(mean(a$total)))

pointsa <- 
  new_data %>%
  add_epred_draws(fita, 
                  re_formula = NA,           
                  allow_new_levels = TRUE) %>% 
  mutate(.epred = .epred/total) %>% 
  group_by(infected, post.treatment) %>% 
  mean_qi(.epred, .width = 0.90) %>% 
  mutate(period = ifelse(post.treatment, "post-treatment", "pre-treatment"))


# plot pre- vs post-treatment----
(median_fita <- 
    pointsa %>%
    ggplot(aes(x = period, y = .epred, shape = infected, color = infected))+  
    geom_jitter(data = a, aes(y= prop), size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
    geom_point(size = 3, position = position_dodge(width = 0.75))+
    geom_errorbar(aes(ymin = .lower,ymax = .upper), 
                  width = 0.1, 
                  size = 0.75, 
                  position = position_dodge(width = 0.75))+
    geom_line(data = pointsa, aes(y = .epred, group = infected), position = position_dodge(width = 0.75))+
    scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
    scale_fill_brewer(palette = "Greys") +
    scale_color_manual(values=c("black", "red"))+
    ylab("proportion of time") +  
    xlab("period") +
    theme_bw() +
    theme(legend.key.size = unit(0.5, "cm"),   
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8), 
          legend.background = element_rect(fill = "transparent"),
          legend.position= "inside",
          legend.position.inside = c(0.75, 0.75)))


# get interaction and effects within each period
t1 <- 
  emmeans(fita, ~ infected * post.treatment, level = 0.90) %>%
  contrast(interaction = "pairwise") %>% 
  tidy() %>% 
  mutate(post.treatment = NA) %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)
t2 <- 
  emmeans(fita, ~ infected | post.treatment, level = 0.90) %>%
  contrast("revpairwise") %>% 
  tidy() %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)

# full plot----
# plot estimates for the interaction, infected pre-treatment, and infected post-treatment
(effectsa <- 
   rbind(t1, t2) %>% 
   mutate(effect = case_when(
     is.na(post.treatment) ~ "interaction term (infected x treatment period)", 
     post.treatment == FALSE ~ "infected effect (pre-treatment)",
     post.treatment == TRUE ~ "infected effect (post-treatment)")) %>% 
   ggplot() +
   geom_hline(yintercept = -1.316, linetype="dashed", 
              color = "purple", size=0.75)+
   geom_hline(yintercept = -0.634, linetype="dashed", 
              color = "orange", size=0.75)+
   geom_pointrange(aes(x = effect, y = estimate, ymin = lower.HPD, ymax = upper.HPD), 
                   linewidth = 1) +
   hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
   coord_flip() +
   theme_bw() +
   theme(legend.position = "none")+
   xlab("")+
   ylab("log odds") +
   ylim(c(-3, 3)) +
   ggtitle("allogrooming"))

# SELF - GROOMING -----

# interaction -----
s <- 
  d2 %>% 
  filter(behav=="self-grooming") %>% 
  mutate(post.treatment = period== "post-treatment")

# model fit
fits <- brm(n | trials(total) ~ infected*post.treatment + (1|bat),
            data = s,
            family = binomial(),
            chains = 4,
            iter = 5000, warmup = 1000)
statss <- 
  as_tibble(summary(fita, prob = 0.90)$fixed) %>% 
  mutate(behavior = "grooming") %>% 
  mutate(model = "interaction") %>% 
  mutate(effect = row.names(summary(fita)$fixed)) %>% 
  select(behavior, model, effect, Estimate:Tail_ESS)

# Calculate mean and 95% credible intervals for each combination of 'infected' and 'period'
new_data <- 
  expand.grid(
    infected = c(TRUE, FALSE),
    post.treatment   = c(TRUE, FALSE),
    total= round(mean(s$total)))

pointss <- 
  new_data %>%
  add_epred_draws(fits, 
                  re_formula = NA,           
                  allow_new_levels = TRUE) %>% 
  mutate(.epred = .epred/total) %>% 
  group_by(infected, post.treatment) %>% 
  mean_qi(.epred, .width = 0.90) %>% 
  mutate(period = ifelse(post.treatment, "post-treatment", "pre-treatment"))


# plot pre- vs post-treatment-----
(median_fits <- 
   pointss %>%
   ggplot(aes(x = period, y = .epred, shape = infected, color = infected))+  
   geom_jitter(data = s, aes(y= prop), size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
   geom_point(size = 3, position = position_dodge(width = 0.75))+
   geom_errorbar(aes(ymin = .lower,ymax = .upper), 
                 width = 0.1, 
                 size = 0.75, 
                 position = position_dodge(width = 0.75))+
   geom_line(data = pointss, aes(y = .epred, group = infected), position = position_dodge(width = 0.75))+
   scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
   scale_fill_brewer(palette = "Greys") +
   scale_color_manual(values=c("black", "red"))+
   ylab("proportion of time") +  
   xlab("period") +
   theme_bw() +
   theme(legend.position = "none"))

# get interaction and effects within each period
t3 <- 
  emmeans(fits, ~ infected * post.treatment, level = 0.90) %>%
  contrast(interaction = "pairwise") %>% 
  tidy() %>% 
  mutate(post.treatment = NA) %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)
t4 <- 
  emmeans(fits, ~ infected | post.treatment, level = 0.90) %>%
  contrast("revpairwise") %>% 
  tidy() %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)

# full plot----
# plot estimates for the interaction, infected pre-treatment, and infected post-treatment
(effectss <- 
   rbind(t3, t4) %>% 
   mutate(effect = case_when(
     is.na(post.treatment) ~ "interaction term (infected x treatment period)", 
     post.treatment == FALSE ~ "infected effect (pre-treatment)",
     post.treatment == TRUE ~ "infected effect (post-treatment)")) %>% 
   ggplot() +
   geom_hline(yintercept = -2.78, linetype="dashed", 
              color = "purple", size=0.75)+
   geom_pointrange(aes(x = effect, y = estimate, ymin = lower.HPD, ymax = upper.HPD), 
                   linewidth = 1) +
   hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
   coord_flip() +
   theme_bw() +
   theme(legend.position = "none")+
   xlab("")+
   ylab("log odds") +
   ylim(c(-3, 3)) +
   ggtitle("self-grooming"))

# MOVING -----

# interaction -----
m <- 
  d2 %>% 
  filter(behav=="moving") %>% 
  mutate(post.treatment = period== "post-treatment")

# model fit
fitm <- brm(n | trials(total) ~ infected*post.treatment + (1|bat),
            data = m,
            family = binomial(),
            chains = 4,
            iter = 5000, warmup = 1000)
statsm <- 
  as_tibble(summary(fitm, prob = 0.90)$fixed) %>% 
  mutate(behavior = "moving") %>% 
  mutate(model = "interaction") %>% 
  mutate(effect = row.names(summary(fita)$fixed)) %>% 
  select(behavior, model, effect, Estimate:Tail_ESS)

# Calculate mean and 95% credible intervals for each combination of 'infected' and 'period'
new_data <- 
  expand.grid(
    infected = c(TRUE, FALSE),
    post.treatment   = c(TRUE, FALSE),
    total= round(mean(m$total)))

pointsm <- 
  new_data %>%
  add_epred_draws(fitm, 
                  re_formula = NA,           
                  allow_new_levels = TRUE) %>% 
  mutate(.epred = .epred/total) %>% 
  group_by(infected, post.treatment) %>% 
  mean_qi(.epred, .width = 0.90) %>% 
  mutate(period = ifelse(post.treatment, "post-treatment", "pre-treatment"))


# plot pre- vs post-treatment-----
(median_fitm <- 
   pointsm %>%
   ggplot(aes(x = period, y = .epred, shape = infected, color = infected))+  
   geom_jitter(data = m, aes(y= prop), size = 2, width = 0.05, height = 0, alpha = 0.5) +   # raw data
   geom_point(size = 3, position = position_dodge(width = 0.75))+
   geom_errorbar(aes(ymin = .lower,ymax = .upper), 
                 width = 0.1, 
                 size = 0.75, 
                 position = position_dodge(width = 0.75))+
   geom_line(data = pointsm, aes(y = .epred, group = infected), position = position_dodge(width = 0.75))+
   scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 17)) +
   scale_fill_brewer(palette = "Greys") +
   scale_color_manual(values=c("black", "red"))+
   ylab("proportion of time") +  
   xlab("period") +
   theme_bw() +
   theme(legend.position = "none"))


# get interaction and effects within each period
t5 <- 
  emmeans(fitm, ~ infected * post.treatment, level = 0.90) %>%
  contrast(interaction = "pairwise") %>% 
  tidy() %>% 
  mutate(post.treatment = NA) %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)
t6 <- 
  emmeans(fitm, ~ infected | post.treatment, level = 0.90) %>%
  contrast("revpairwise") %>% 
  tidy() %>% 
  select(post.treatment, term, estimate, lower.HPD, upper.HPD)

# full plot----
# plot estimates for the interaction, infected pre-treatment, and infected post-treatment
(effectsm <- 
   rbind(t5, t6) %>% 
   mutate(effect = case_when(
     is.na(post.treatment) ~ "interaction term (infected x treatment period)", 
     post.treatment == FALSE ~ "infected effect (pre-treatment)",
     post.treatment == TRUE ~ "infected effect (post-treatment)")) %>% 
   ggplot() +
   geom_hline(yintercept = -2.78, linetype="dashed", 
              color = "purple", size=0.75)+
   geom_pointrange(aes(x = effect, y = estimate, ymin = lower.HPD, ymax = upper.HPD), 
                   linewidth = 1) +
   hline_at(v = 0, alpha = 0.5, linewidth = 0.5, linetype = 1) +
   coord_flip() +
   theme_bw() +
   theme(legend.position = "none")+
   xlab("")+
   ylab("log odds") +
   ylim(c(-3, 3)) +
   ggtitle("moving"))

# plot figure -----
p1 <- plot(effectsa + median_fita)
p2 <- plot(effectss + median_fits)
p3 <- plot(effectsm + median_fitm)


plot1 <- plot(p1/ p2/ p3) + plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "")))

# save to PDF
ggsave(plot1, filename= "Fig2_sickness_behavior.pdf", width = 9, height = 7)

# save results -----
brms_results1 <- 
  rbind(t1, t2, t3, t4, t5, t6) %>% 
  mutate(model= rep(c("allogrooming", "self-grooming", "moving"), each=3)) %>% 
  mutate(effect = case_when(
    is.na(post.treatment) ~ "interaction term (infected x treatment period)", 
    post.treatment == FALSE ~ "infected effect (pre-treatment)",
    post.treatment == TRUE ~ "infected effect (post-treatment)")) %>% 
  select(model, effect, estimate, lower.HPD, upper.HPD)
write.csv(brms_results1, "brms_results1.csv")

brms_results2 <- rbind(statsa, statss, statsm)
write.csv(brms_results2, "brms_results2.csv")

