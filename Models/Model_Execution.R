################################################################################################################################################################################

#Script for model execution for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Library packages to load

library(rethinking)

#Load data 

setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))

SSL_d <- read.csv("Data/SSL_Test_Data_Processed.csv") 
IMP_d <- read.csv("Data/SSL_IMP_Data_Processed.csv") 

#Minor data wrangling for first-choice model

FC_d <- subset(SSL_d, trial == 1)            #Keep only rows where trial is 1
FC_d$choice[FC_d$choice == 2] <- 0           #Convert choice 2 to 0

#Build lists for all models...

FC_stan_list <- list( #First choice
  N = nrow(FC_d),
  treat = as.integer(FC_d$treat),  
  choice = as.integer(FC_d$choice)
)

All_choice_stan_list <- list( #All choices
  N = nrow(SSL_d),
  N_id = length(unique(SSL_d$id)),
  N_treat = length(unique(SSL_d$treat)),
  id = SSL_d$id,
  treat = SSL_d$treat,
  choice = as.integer(ifelse(SSL_d$choice == 2, 0, 1)), #Convert choice 2 to 0
  trial = SSL_d$trial
)

EWA_stan_list <- list( #Baseline & monotonic
  N = nrow(SSL_d),
  N_id = length(unique(SSL_d$id)),
  id = SSL_d$id,
  manip = SSL_d$Experiment,
  sat = SSL_d$Satisfaction,
  treat = SSL_d$treat,
  trial = SSL_d$trial,
  choice = SSL_d$choice,
  cum_soc = SSL_d$cum_soc,
  cum_non_soc = SSL_d$cum_non_soc,
  social = SSL_d$social, 
  touch_P = IMP_d$IMP_P_TT_Count,
  touch_O = IMP_d$IMP_O_TT_Count,
  log_dur = IMP_d$IMP_Sec
)

#Quick looks

str(FC_stan_list)
str(All_choice_stan_list)
str(EWA_stan_list)

#
##
###Run models
##
#

FC_m <- cstan( #First choice
  file = "Models/LR_FC_Model.stan", 
  data = FC_stan_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

AC_m <- cstan( #All choices
  file = "Models/LR_All_Choices_Model.stan", #Your path here
  data = All_choice_stan_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

EWA_b <- cstan( #EWA baseline model
  file = "Models/EWA_Baseline_Model.stan", #Your path here
  data = EWA_stan_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

EWA_m <- cstan( #EWA monotonic model
  file = "Models/EWA_Monotonic_Model.stan", #Your path here
  data = EWA_stan_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

#
##
### PRE-STUDY MODEL VALIDATION
##
#

#Note Base_And_Mono_Pre_Study_Sim.R script must have been run prior to executing the following!

d_base <- Pre_study_sim_fct(Mono = 0)
d_mono <- Pre_study_sim_fct(Mono = 1)

#Put both data sets into lists for Stan

EWA_PS_base_list <- list(
  N = nrow(d_base),
  N_id = length(unique(d_base$id)),
  id = as.integer(d_base$id),
  manip = d_base$experiment,
  sat = as.integer(d_base$sat_level),
  trial = d_base$trial,
  choice = d_base$choice,
  cum_soc = as.integer(d_base$cum_soc),
  cum_non_soc = as.integer(d_base$cum_non_soc),
  social = as.integer(d_base$social), 
  touch_P = sample(10:100, 47, replace = TRUE),
  touch_O = sample(10:100, 47, replace = TRUE),
  log_dur = log(runif(47, 14400, (14400 * 4))) 
)

EWA_PS_mono_list <- list(
  N = nrow(d_mono),
  N_id = length(unique(d_mono$id)),
  id = as.integer(d_mono$id),
  manip = d_mono$experiment,
  sat = as.integer(d_mono$sat_level),
  trial = d_mono$trial,
  choice = d_mono$choice,
  cum_soc = as.integer(d_mono$cum_soc),
  cum_non_soc = as.integer(d_mono$cum_non_soc),
  social = as.integer(d_mono$social), 
  touch_P = sample(10:100, 47, replace = TRUE),
  touch_O = sample(10:100, 47, replace = TRUE),
  log_dur = log(runif(47, 14400, (14400 * 4))) 
)

#Quick looks

str(EWA_PS_base_list)
str(EWA_PS_mono_list)

#Run EWA models on simulated data

PS_EWA_b <- cstan(
  file = "Models/EWA_Baseline_Model.stan", #Your path here
  data = EWA_PS_base_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

PS_EWA_m <- cstan(
  file = "Models/EWA_Monotonic_Model.stan", #Your path here
  data = EWA_PS_mono_list,
  cores = 4,
  chains = 4,
  refresh = 10,
  iter = 2000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)
