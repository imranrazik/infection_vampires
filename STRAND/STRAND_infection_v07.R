# first, be sure to install "RStan" and "CmdStanR" from the following link if you have not done so already: https://mc-stan.org/

# Clear working space
rm(list = ls()) 

# set working directory
setwd(dirname(file.choose()))

# Load libraries ----
library(tidyverse)
library(devtools)
library(STRAND)
library(boot)

# functions-----
# bootstrap the mean of a vector
# get 95% CI around the mean of a vector by bootstrapping
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


# get mean and 95% CI of values y within grouping variable x via bootstrapping 
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

# plot permutation test results as histogram
hist_perm <- function(exp=exp, obs=obs, perms=perms, label=''){
  exp.range <- round(quantile(exp, probs= c(0.025, 0.975)),3)
  ggplot()+
    geom_histogram(aes(x=exp), color="black",fill="light blue")+
    geom_vline(aes(xintercept=obs), color="red", size=1)+
    xlab("expected values from null model")+
    ggtitle(label, subtitle = paste('obs = ',round(obs,3), ', exp = ', exp.range[1], ' to ', exp.range[2], ", Prob exp >= obs: p", ifelse(mean(exp>=obs)==0,paste("<",1/perms), paste("=",signif(mean(exp>=obs),digits=2))),", permutations=",perms, sep=""))
}

# Load data -----
# load("rates2019_2021-05-18.RData")
load("STRAND_results_v07.RData") 

if (!exists("res")){ # if the STRAND results objects doesn't exist, run the code (takes several days)
  # wrangle data ----
  
  # get the data from the forced proximity phase (phase2)
  small_cages <- 
    rates2019 %>% 
    filter(phase == 2) %>% 
    group_by(actor) %>% 
    summarize(cage= first(cage))
  
  # get bats with population and forced proximity cages
  bats <- 
    rates2019 %>% 
    # remove bats that are not subjects or were not sampled during each phase of the experiment
    filter(!actor %in% c("ald", "ad", "bd", "bscs", "cwd")) %>% 
    group_by(actor) %>% 
    summarize(population= first(actor.pop)) %>% 
    filter(population != "captive-born") %>% 
    mutate(cage = small_cages$cage[match(.$actor, small_cages$actor)]) %>% 
    arrange(cage, actor) %>% 
    ungroup() %>% 
    mutate(pop.cage = paste(population, cage, sep = ".")) 
  
  # get simple list of bat IDs
  bats2 <- 
    bats %>% 
    arrange(actor) %>% 
    pull(actor) %>%
    unique() 
  
  # get infected bats
  treated <- c("adldd", "aldd", "adld", "blx", "cdw", "cww", "cfd", "cldw", "cnone")
  
  # get injured bats
  injured <- c("cww", "cfd", "cldw", "cnone")
  
  # get grooming durations -----
  d <- 
    rates2019 %>% 
    # get only bats of interest (adult females sampled throughout each phase of the study)
    filter(behav == "g") %>% 
    filter(actor %in% bats$actor & receiver %in% bats$actor) %>% 
    # remove data from phase2, which is in a different format
    filter(phase!="phase2") %>% 
    # get new measures of time at the level of sampling period, in datetime format
    mutate(hour2 = gsub('^(.{2})(.*)$', '\\1:\\2', hour)) %>% 
    mutate(hour2 = gsub('^(.{5})(.*)$', '\\1:00\\2', hour2)) %>% 
    mutate(period2 = paste(date, hour2, sep = " ")) %>% 
    mutate(period2 = as.POSIXct(period2, format="%Y-%m-%d %H:%M:%S", tz=Sys.timezone())) %>% 
    # relabel phases with text instead
    mutate(phase = case_when(
      phase == 2 ~ "phase2",
      phase == 1 ~ "phase1",
      phase == 3 ~ "phase3")) %>% 
    mutate(bins = case_when(
      period2 > '2019-06-23 00:00:00' & period2 < '2019-06-29 11:59:59' ~ 1, 
      period2 > '2019-06-30 00:00:00' & period2 < '2019-07-06 11:59:59' ~ 2, 
      period2 > '2019-07-07 00:00:00' & period2 < '2019-07-13 11:59:59' ~ 3, 
      period2 > '2019-07-14 00:00:00' & period2 < '2019-07-20 11:59:59' ~ 4, 
      period2 > '2019-07-21 00:00:00' & period2 < '2019-07-27 11:59:59' ~ 5, 
      period2 > '2019-07-28 00:00:00' & period2 < '2019-08-03 11:59:59' ~ 6, 
      period2 > '2019-08-04 00:00:00' & period2 < '2019-08-10 11:59:59' ~ 7, 
      period2 > '2019-08-11 00:00:00' & period2 < '2019-08-17 11:59:59' ~ 8, 
      period2 > '2019-08-18 00:00:00' & period2 < '2019-08-24 11:59:59' ~ 9, 
      period2 > '2019-08-25 00:00:00' & period2 < '2019-08-31 11:59:59' ~ 10, 
      period2 > '2019-09-01 00:00:00' & period2 < '2019-09-07 11:59:59' ~ 11,
      period2 > '2019-09-08 00:00:00' & period2 < '2019-09-14 11:59:59' ~ 12,
      period2 > '2019-09-15 00:00:00' & period2 < '2019-09-21 11:59:59' ~ 13,
      period2 > '2019-09-22 00:00:00' & period2 < '2019-09-28 11:59:59' ~ 14,
      period2 > '2019-09-29 00:00:00' & period2 < '2019-10-05 11:59:59' ~ 15,
      period2 > '2019-10-06 00:00:00' & period2 < '2019-10-14 11:59:59' ~ 16 # this one includes two extra days
    )) %>%
    # filter out the bin that matches up with phase 2 of the experiment. 
    # There are only a few interactions on 20190804, from right before the bats were transferred to the forced proximity trials
    filter(!bins == 7) %>% 
    # label dyads from forced proximity phase (phase2)
    mutate(a.cage= bats$cage[match(.$actor, bats$actor)],
           r.cage= bats$cage[match(.$receiver, bats$actor)]) %>% 
    mutate(dyad.type= case_when(
      a.cage == r.cage & new.dyad ~ "test dyads",
      a.cage == r.cage & !new.dyad ~ "something_went_wrong",
      a.cage != r.cage & new.dyad ~ "control dyads",
      a.cage != r.cage & !new.dyad ~ "familiar dyads")) %>% 
    # label status of bats as asymptomatic (0), sick and treated (1), sick with severe symptoms (2)
    mutate(actor.status= case_when(
      actor %in% injured ~ "injured",
      actor %in% treated ~ "treated",
      T ~ "asymptomatic"
    )) %>% 
    mutate(receiver.status= case_when(
      receiver %in% injured ~ "injured",
      receiver %in% treated ~ "treated",
      T ~ "asymptomatic"
    )) %>% 
    mutate(actor.infection = case_when(
      actor.status == "asymptomatic" ~ 0,
      actor.status == "treated" ~ 1,
      actor.status == "injured" ~ 2
    )) %>% 
    mutate(receiver.infection = case_when(
      receiver.status == "asymptomatic" ~ 0,
      receiver.status == "treated" ~ 1,
      receiver.status == "injured" ~ 2
    )) %>% 
    # label undirected dyads
    mutate(dyad2= ifelse(actor<receiver,
                         paste(actor, receiver, sep="_"),
                         paste(receiver, actor, sep="_"))) %>%
    # convert grooming rates to simply presence/absence of grooming per sampling period
    mutate(rate2 = if_else(rate > 0, 1, rate)) %>% 
    # relevel dyad types
    mutate(dyad.type = as_factor(dyad.type)) %>% 
    mutate(dyad.type = fct_relevel(dyad.type, c("familiar dyads", "control dyads", "test dyads"))) %>% 
    # create combination of dyad and period
    mutate(dyad.period2 = paste(dyad, period2, sep = "_")) 
  
  # get edgelist with missing values for grooming durations
  pos.rates <- expand.grid(bats2, bats2, unique(d$date))
  pos.rates <- pos.rates %>% 
    rename(actor = Var1, receiver = Var2, date = Var3) %>% 
    arrange(actor) %>% 
    mutate(dyad = paste(actor, receiver, sep = "_")) %>% 
    mutate(dyad.date = paste(dyad, date, sep = "_"))
  pos.rates$rate <- NA
  
  # insert actual grooming durations
  t <- 
    d %>% 
    group_by(date, bins, actor, receiver, dyad.type, new.dyad, actor.infection, receiver.infection) %>% 
    summarise(rate = sum(rate, na.rm =T), exposure = length(unique(period))*60*60) %>% 
    #mutate(log.rate = log(rate + 1)) %>% 
    mutate(dyad = paste(actor, receiver, sep = "_")) %>% 
    mutate(dyad.date = paste(dyad, date, sep = "_"))
  
  # match actual durations into pos.rates
  pos.rates$rate <- t$rate[match(pos.rates$dyad.date, t$dyad.date)]
  
  pos.rates$rate[which(is.na(pos.rates$rate))] <- 0
  
  # match exposures
  pos.rates$exposure <- t$exposure[match(pos.rates$dyad.date, t$dyad.date)]
  
  pos.rates$exposure[which(is.na(pos.rates$exposure))] <- 0
  
  # match bins
  pos.rates$bins <- t$bins[match(pos.rates$date, t$date)]
  
  # add in other predictors
  # dyad type
  pos.rates$dyad.type <- t$dyad.type[match(pos.rates$dyad, t$dyad)]
  pos.rates$dyad.type <- if_else(is.na(pos.rates$dyad.type), "none", pos.rates$dyad.type)
  
  # familiarity
  pos.rates$new.dyad <- t$new.dyad[match(pos.rates$dyad, t$dyad)]
  pos.rates$new.dyad <- if_else(is.na(pos.rates$new.dyad), FALSE, pos.rates$new.dyad)
  
  # forced into proximity
  pos.rates$forced <- if_else(pos.rates$dyad.type == "test dyads", 1, 0)
  
  # actor infection 
  pos.rates$actor.infection <- t$actor.infection[match(pos.rates$actor, t$actor)]
  
  # receiver infection
  pos.rates$receiver.infection <- t$receiver.infection[match(pos.rates$receiver, t$receiver)]
  
  # sum the scale of infection across actor and receiver
  rates <- 
    pos.rates %>% 
    mutate(infection = (actor.infection + receiver.infection)) %>% 
    as_tibble()
  
  # Make sure levels are the same between actor and receiver
  actor <- 
    rates %>% 
    arrange(actor) %>% 
    pull(actor) %>% 
    unique()
  
  rates$actor <- factor(rates$actor, levels = actor)
  rates$receiver <- factor(rates$receiver, levels = actor)
  
  # weekly time bins
  
  # get vector of bins
  bins <- unique(rates$bins)
  
  # create empty lists of 16 elements to populate with results across 16 weeks
  results_dyads <- vector("list", length = 16) # dyadic effects
  results_focal <- vector("list", length = 16) # focal effects
  results_target <- vector("list", length = 16) # target effects
  results_other <- vector("list", length = 16) # other estimates (generalized reciprocity and dyadic affinity)
  
  diag_summary <- vector("list", length = 16)
  diag_dyad_effects <- vector("list", length = 16)
  plot_list <- vector("list", length = 16)
  
  for (i in c(1:6, 8:16)) {
    
    # filter data by week/bin
    t <- 
      rates %>% 
      filter(bins == i) %>% 
      group_by(actor, receiver, dyad, dyad.type, new.dyad, forced, infection) %>% 
      summarise(rate = sum(rate), exposure = sum(exposure))
    
    # grooming matrix
    m.groom <- matrix(t$rate, nrow = 21, ncol = 21) 
    rownames(m.groom) <- unique(t$actor)
    colnames(m.groom) <- unique(t$actor)
    m.groom <- as.matrix(m.groom) # this may be unecessary
    # actors are along rows and receivers along columns
    
    # exposure matrix
    m.periods <- matrix(t$exposure, nrow = 21, ncol = 21)
    rownames(m.periods) <- unique(t$actor)
    colnames(m.periods) <- unique(t$actor)
    m.periods <- as.matrix(m.periods)
    
    # new dyad matrix
    t$new.dyad <- as.numeric(t$new.dyad)
    m.newdyads <- matrix(t$new.dyad, nrow = 21, ncol = 21)
    rownames(m.newdyads) <- unique(t$actor)
    colnames(m.newdyads) <- unique(t$actor)
    m.newdyads <- as.matrix(m.newdyads)
    
    # forced proximity matrix
    m.forced <- matrix(t$forced, nrow = 21, ncol = 21)
    rownames(m.forced) <- unique(t$actor)
    colnames(m.forced) <- unique(t$actor)
    m.forced <- as.matrix(m.forced)
    
    # infection matrix
    m.infection <- matrix(t$infection, nrow = 21, ncol = 21)
    rownames(m.infection) <- unique(t$actor)
    colnames(m.infection) <- unique(t$actor)
    m.infection <- as.matrix(m.infection)
    
    # create data object
    Bat_Data_2019 <- list(groom = m.groom, periods = m.periods, new.dyad = m.newdyads, forced = m.forced, infection = m.infection)
    #Bat_Data_2019 <- list(groom = m.groom, periods = m.periods, dyad.type = m.dyad.type, pop = pop)
    
    # Run the STRAND model -----
    
    # grooming network
    nets = list(groom = Bat_Data_2019$groom)
    
    # Dyadic variables
    dyad = list(new.dyad = Bat_Data_2019$new.dyad, forced = Bat_Data_2019$forced, infection = Bat_Data_2019$infection)
    
    # exposure
    exposure = list(periods = Bat_Data_2019$periods)
    
    # get model data together
    model_dat = make_strand_data(self_report = nets,
                                 #block_covariates = group_ids,
                                 individual_covariates = NULL, 
                                 dyadic_covariates = dyad,
                                 outcome_mode = "binomial",
                                 exposure = exposure
    )
    
    # Model
    fit = fit_social_relations_model(data=model_dat,
                                     #block_regression = ~ sickness,
                                     focal_regression = ~ 1,
                                     target_regression = ~ 1,
                                     dyad_regression = ~ new.dyad*infection + forced,
                                     mode="mcmc",
                                     return_predicted_network = TRUE,
                                     stan_mcmc_parameters = list(chains = 1,
                                                                 iter_warmup = 2000, iter_sampling = 15000,
                                                                 max_treedepth = NULL, adapt_delta = .90)
    )
    
    res = summarize_strand_results(fit)
    
    # extract dataframes of results
    t2 <- as_tibble(res$summary_list$`Dyadic effects`)
    t3 <- as_tibble(res$Qsummary_list$`Focal effects: Out-degree`)
    t4 <- as_tibble(res$summary_list$`Target effects: In-degree`)
    t5 <- as_tibble(res$summary_list$`Other estimates`)
    
    t6 <- as_tibble(fit$fit$summary())
    t7 <- as_tibble(fit$fit$summary("dyad_effects"))
    
    # store dataframes in their respective lists
    results_dyads[[i]] <- t2
    results_focal[[i]] <- t3
    results_target[[i]] <- t4
    results_other[[i]] <- t5
    
    diag_summary[[i]] <- t6
    diag_dyad_effects[[i]] <- t7
    bayesplot::color_scheme_set("mix-blue-red")
    plot_list[[i]] <- bayesplot::mcmc_trace(fit$fit$draws(), pars = c("dyad_effects[1]", "dyad_effects[2]", "dyad_effects[3]", "dyad_effects[4]", "sr_L[2,1]","dr_L[2,1]"))
  }
}else{NULL}

# Extract diagnostics for correlations between sender/receiver random effects and dyadic random effects -----
diag_summary_list <- list()

for (i in 1:length(diag_summary)){
  
  if (is.null(diag_summary[[i]])) {
    # Extract column names
    column_names <- names(diag_summary[[1]])
    
    # Create a new dataframe with one row filled with NA values
    t <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
    
    # Assign column names to the new dataframe
    colnames(t) <- column_names
    
    # specify the week number / time bin
    diag_summary_list[[i]] <-
      t %>% 
      mutate(bin = paste(i))
    
  }else{
    # sender/receiver random effects correlation
    sr <- 
      diag_summary[[i]][6,] %>% 
      mutate(bin = paste(i))
    
    # dyadic random effects correlation
    dr <- 
      diag_summary[[i]][53,] %>% 
      mutate(bin = paste(i))
    
    t <- rbind(sr, dr)
    
    diag_summary_list[[i]] <- t
  }
}

diag_summary_df <- bind_rows(diag_summary_list)

# save
write.csv(diag_summary_df, "Random_effects_correlation_diagnostics.csv", row.names = FALSE)

# Extract diagnostics for fixed dyadic effects -----
diag_dyad_effects_list <- list()

for (i in 1:length(diag_dyad_effects)){
  
  if (is.null(diag_dyad_effects[[i]])) {
    # Extract column names
    column_names <- names(diag_dyad_effects[[1]])
    
    # Create a new dataframe with one row filled with NA values
    t <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
    
    # Assign column names to the new dataframe
    colnames(t) <- column_names
    
    # specify week number 
    diag_dyad_effects_list[[i]] <-
      t %>% 
      mutate(bin = paste(i))
    
  }else{
    t <- 
      diag_dyad_effects[[i]][1:4,] %>% # we have four estimates from the model
      mutate(bin = paste(i))
    
    diag_dyad_effects_list[[i]] <- t
  }
}

diag_dyad_effects_df <- bind_rows(diag_dyad_effects_list)

# save
write.csv(diag_dyad_effects_df, "Dyadic_effects_diagnostics.csv", row.names = FALSE)

# Extract diagnostics plots -----
t <- plot_list[[2]]

# PLOT EFFECTS THROUGH TIME -----
# new dyad estimate ------
# create vector to new.dyad coefficient medians of each week
new.dyad <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][2, 2])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][2, 2])))
  }
  # Store the value in the result vector
  new.dyad[i] <- value
}

# low HDPI
new.dyadlow <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][2, 3])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][2, 3])))
  }
  # Store the value in the result vector
  new.dyadlow[i] <- value
}

# high HDPI
new.dyadhigh <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][2, 4])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][2, 4])))
  }
  # Store the value in the result vector
  new.dyadhigh[i] <- value
}

t <- tibble(week = 1:16, median = new.dyad, low = new.dyadlow, high = new.dyadhigh)
t %>% 
  ggplot(aes(x = week, y = median))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  geom_hline(aes(yintercept = 0))+
  ggtitle("effect of new dyad")

# infection estimate ------
infection <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][3, 2])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][3, 2])))
  }
  # Store the value in the result vector
  infection[i] <- value
}

infectionlow <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][3, 3])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][3, 3])))
  }
  # Store the value in the result vector
  infectionlow[i] <- value
}

infectionhigh <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][3, 4])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][3, 4])))
  }
  # Store the value in the result vector
  infectionhigh[i] <- value
}

t2 <- tibble(week = 1:16, median = infection, low = infectionlow, high = infectionhigh)
t2 %>% 
  ggplot(aes(x = week, y = median))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  geom_hline(aes(yintercept = 0))+
  ggtitle("effect of infection")

# forced estimate -----
forced <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][4, 2])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][4, 2])))
  }
  # Store the value in the result vector
  forced[i] <- value
}

forcedlow <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][4, 3])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][4, 3])))
  }
  # Store the value in the result vector
  forcedlow[i] <- value
}

forcedhigh <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][4, 4])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][4, 4])))
  }
  # Store the value in the result vector
  forcedhigh[i] <- value
}

t3 <- tibble(week = 1:16, median = forced, low = forcedlow, high = forcedhigh)
t3 %>% 
  ggplot(aes(x = week, y = median))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  geom_hline(aes(yintercept = 0))+
  ggtitle("effect of forced proximity")

# new.dyad:infection estimate -----
newinfection <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][5, 2])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][5, 2])))
  }
  # Store the value in the result vector
  newinfection[i] <- value
}

newinfectionlow <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][5, 3])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][5, 3])))
  }
  # Store the value in the result vector
  newinfectionlow[i] <- value
}

newinfectionhigh <- numeric(length(results_dyads))
# Loop through each data frame in the list
for (i in 1:length(results_dyads)) {
  if (length(as.numeric(unname(unlist(results_dyads[[i]][5, 4])))) == 0) {
    value <- NA
  }else{
    # Extract the value at the specified row and column
    value <- as.numeric(unname(unlist(results_dyads[[i]][5, 4])))
  }
  # Store the value in the result vector
  newinfectionhigh[i] <- value
}

t4 <- tibble(week = 1:16, median = newinfection, low = newinfectionlow, high = newinfectionhigh)
t4 %>% 
  ggplot(aes(x = week, y = median))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  geom_hline(aes(yintercept = 0))+
  ylab("posterior probability")+
  xlab("time (weeks)")

# save the values for estimates across weeks 
write.csv(t4, "interaction_estimates.csv", row.names = FALSE)

ggsave(
  'Fig2.pdf',
  plot = last_plot(),
  scale = 1,
  width = 6,
  height = 3,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

ggsave(
  'Fig2.png',
  plot = last_plot(),
  scale = 1,
  width = 6,
  height = 3,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# calculate slope for change in effect of interaction through time 
summary(glm(data = t4, median ~ week))

# interpret interaction between new.dyad and infection -----
theme_set(theme_bw(base_size = 14))
(plot <- 
    d %>% 
    group_by(actor, receiver, dyad, dyad2, date, bins, new.dyad, actor.infection, receiver.infection) %>% 
    summarise(rate = mean(rate, na.rm = T), .groups = "drop") %>% 
    group_by(actor, receiver, dyad, dyad2, bins, new.dyad, actor.infection, receiver.infection) %>% 
    summarise(rate = mean(rate, na.rm = T), .groups = "drop") %>% 
    mutate(infection = actor.infection + receiver.infection) %>% 
    group_by(dyad2, bins, new.dyad, infection) %>% 
    summarise(rate = mean(rate, na.rm = T), .groups = "drop") %>% 
    group_by(new.dyad, bins, infection) %>% 
    summarise(rate = mean(rate, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(lograte = log(rate + 1)) %>% 
    mutate(new.dyad = if_else(new.dyad == T, "new dyad", "familiar dyad")) %>% 
    mutate(infection2= as.character(infection>=1)) %>% 
    ggplot(aes(x = bins, y = lograte, color = infection2, shape= infection2, linetype= infection2))+
    geom_smooth(method = "lm", size = 2, se = T)+
    geom_point(size=2)+
    facet_wrap(~new.dyad)+
    ylab("log grooming rate")+
    xlab('time (week)')+
    scale_color_manual(values=c(
      'orange',
      'darkred'))+
    theme(legend.position= "none"))

# save plot
ggsave(
  'Fig3.pdf',
  plot = last_plot(), 
  scale = 1,
  width = 6,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

ggsave(
  'Fig3.png',
  plot = last_plot(),
  scale = 1,
  width = 6,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)



 


