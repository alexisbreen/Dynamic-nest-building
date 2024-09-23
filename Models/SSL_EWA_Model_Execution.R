################################################################################################################################################################################

#Stan EWA model execution script for the manuscript

#Strategic social learning in avian nest construction and potentially beyond

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de)

################################################################################################################################################################################

#Library packages to load

library(tidyverse)
library(rethinking)

#Load data 

SSL_d <- read.csv(file.choose(), header = T) #file: SSL_Test_Data_Processed
IMP_d <- read.csv(file.choose(), header = T) #file: SSL_IMP_Data_Processed

#Put data into a list for Stan

d_stan_list <- list(
  
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

#Run stan model

m <- cstan(file = "SSL_EWA_Model.stan", data = d_stan_list, cores = 4, chains = 4, refresh = 10, iter = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))

#Inspect effective model samples 

precis(m, depth = 3, pars = c("phi", "lambda", "rho", "sigma"))

#Extract posterior

s <- extract.samples(m) 

#Extract latent parameters from posterior

#Phi

s_phi <- as.data.frame(
  list(
    S_1 = s$phi[ , 1, 1], #Exp 1, sat
    D_1 = s$phi[ , 1, 2], #Exp 1, diss
    S_2 = s$phi[ , 2, 1], #Exp 2, sat
    D_2 = s$phi[ , 2, 2]  #Exp 2, diss
  )
)

#Lambda

s_lambda <- as.data.frame(
  list(
    S_1 = s$lambda[ , 1, 1], #Exp 1, sat
    D_1 = s$lambda[ , 1, 2], #Exp 1, diss
    S_2 = s$lambda[ , 2, 1], #Exp 2, sat
    D_2 = s$lambda[ , 2, 2]  #Exp 2, diss
  )
)

#Rho

s_rho <- as.data.frame(
  list(
    S_1 = s$rho[ , 1, 1], #Exp 1, sat
    D_1 = s$rho[ , 1, 2], #Exp 1, diss
    S_2 = s$rho[ , 2, 1], #Exp 2, sat
    D_2 = s$rho[ , 2, 2]  #Exp 2, diss
  )
)

#Sigma

s_sig <- as.data.frame(
  list(
    S_1 = s$sigma[ , 1, 1], #Exp 1, sat
    D_1 = s$sigma[ , 1, 2], #Exp 1, diss
    S_2 = s$sigma[ , 2, 1], #Exp 2, sat
    D_2 = s$sigma[ , 2, 2]  #Exp 2, diss
  )
)

