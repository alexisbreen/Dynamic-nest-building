################################################################################################################################################################################

#Script for model prediction model power checks (Pareto-Smoothed Importance Sampling Leave-Future-Out Cross-Validation) for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Library packages to load
library(cmdstanr)           #Interface for fitting Stan models using CmdStan backend
library(loo)                #For approximate leave-one-out cross-validation diagnostics
library(MASS)               #Includes statistical functions and datasets e.g., multivariate normal
library(posterior)          #Tools for working with posterior draws from Bayesian models
library(rethinking)         #Companion package to the Statistical Rethinking book by McElreath; includes convenient Bayesian wrappers

#Note Model_Execution.R script must have been run to run the below!

#Configuration
L <- 5                      #Number of future trials to consider in leave-future-out cross-validation
k_thresh <- 0.7             #Threshold for Pareto-k diagnostic to identify bad importance sampling estimates
M <- 1                      #Number of chains or possible replications 
test_birds <- 47            #Bird ID index being used for prediction 

#Helper functions
log_mean_exp <- function(x) {
  max_x <- max(x)                                 #Numerical stabilisation: subtract max to avoid exponent overflow
  return(max_x + log(mean(exp(x - max_x))))       #Return log of mean exponentiated values adjusted for max
}

log_sum_exp <- function(x) {
  max_x <- max(x)                                 #Numerical stabilisation same as above
  max_x + log(sum(exp(x - max_x)))                #Return log of sum of exponentials adjusted for max
}

#Stan list for running the models
make_stan_data_common <- function(bird_id, bird_data, stan_list) {
  list(
    N = nrow(bird_data),                            #Number of trials for this bird
    N_id = 1,                                       #Single bird at a time is being passed to Stan
    id = rep(1, nrow(bird_data)),                   #Vector of IDs 
    manip = bird_data$manip,                        #Manipulation condition (1 = Material, 2 = Breeding)
    sat = bird_data$sat,                            #Satisfaction condition (1 = Satisfied, 2 = Dissatisfied)
    treat = bird_data$treat,                        #Combined treatment code 
    trial = bird_data$trial,                        #Trial number 
    choice = bird_data$choice,                      #Choice made by the bird 
    cum_soc = bird_data$cum_soc,                    #Cumulative social information used up to this trial
    cum_non_soc = bird_data$cum_non_soc,            #Cumulative asocial information used up to this trial
    social = bird_data$social                       #Binary indicator whether this trial had social info

  )
}

#Extract baseline target parameters
extract_baseline_components <- function(posterior) {
  posterior_mat <- as.matrix(as.data.frame(posterior))                             #Convert posterior to matrix for indexing
  list(
    logit_phi = posterior_mat[, grep("^logit_phi", colnames(posterior_mat))],      #Latent asocial information appetite
    log_lambda = posterior_mat[, grep("^log_lambda", colnames(posterior_mat))],    #Latent asocial weighting
    log_epsilon = posterior_mat[, grep("^log_epsilon", colnames(posterior_mat))],  #Latent social appetite
    logit_sigma = posterior_mat[, grep("^logit_sigma", colnames(posterior_mat))],  #Latent social weighting
    z_ID = posterior_mat[, grep("^z_ID", colnames(posterior_mat))],                #ID-specific offsets
    sigma_ID = posterior_mat[, grep("^sigma_ID", colnames(posterior_mat))],        #Scale of ID-specific variability
    Rho_ID = posterior_mat[, grep("^Rho_ID", colnames(posterior_mat))],            #Correlation matrix between latents
    A_init = posterior_mat[, grep("^A_init", colnames(posterior_mat))]             #Initial attractions per option
  )
}

#Custom log-likelihood function for baseline
compute_log_lik_baseline <- function(posterior, trial_row, id = 1) {
  pars <- extract_baseline_components(posterior)                                   #Extract all posterior parameters needed
  S <- nrow(posterior)                                                             #Number of posterior samples
  trial <- trial_row$trial                                                         #Current trial number
  log_lik <- numeric(S)                                                            #Preallocate log-likelihood vector
  
  for (s in 1:S) {
    offset <- pars$z_ID[s, ((id - 1) * 4 + 1):(id * 4)]                            #ID-specific offsets for phi, lambda, epsilon, sigma
    Sigma_ID <- diag(as.numeric(pars$sigma_ID[s, 1:4]))                            #Diagonal covariance matrix from std devs
    Rho <- matrix(as.numeric(pars$Rho_ID[s, 1:(4 * 4)]), 4, 4)                     #Correlation matrix among parameters
    L <- t(chol(Sigma_ID %*% Rho %*% Sigma_ID))                                    #Cholesky decomposition of covariance matrix
    v_ID <- L %*% offset                                                           #Draw of latent parameters for this ID
    
    idx <- (trial_row$manip - 1) * 2 + trial_row$sat                               #Index to access treatment-specific latents
    lambda <- min(exp(pars$log_lambda[s, idx] + v_ID[1]), 15)                      #Asocial weighting bounded above at 15
    epsilon <- min(exp(pars$log_epsilon[s, idx] + v_ID[2]), 10)                    #Social appetite bounded above at 10
    sigma <- inv_logit(pars$logit_sigma[s, idx] + v_ID[3])                         #Social weighting transformed via logistic
    phi <- inv_logit(pars$logit_phi[s, idx] + v_ID[4])                             #Asocial appetite transformed via logistic
    
    A <- inv_logit(pars$A_init[s, ((id - 1) * 2 + 1):(id * 2)])                    #Initial attraction to each option
    A <- sort(A)                                                                   #Ensure consistent ordering of options
    pA <- softmax(lambda * A)                                                      #Choice probabilities from attractions
    
    if (trial_row$social == 0) {
      pC <- pA                                                                     #If no social info, use asocial probabilities
    } else {
      safe_cum_soc <- max(trial_row$cum_soc, 1e-6)                                 #Avoid zero division by setting lower bound
      safe_cum_non_soc <- max(trial_row$cum_non_soc, 1e-6)                         #Same for asocial cumulative
      pS <- c(
        safe_cum_soc^epsilon / (safe_cum_soc^epsilon + safe_cum_non_soc^1),
        safe_cum_non_soc^1 / (safe_cum_soc^epsilon + safe_cum_non_soc^1)
      )                                                                            #Social choice probs based on cumulative counts
      pC <- (1 - sigma) * pA + sigma * pS                                          #Weighted mix of asocial and social probs
    }
    
    log_lik[s] <- log(pC[trial_row$choice])                                        #Log-likelihood of observed choice
  }
  
  return(log_lik)                                                                  #Return vector of log likelihoods over samples
}

#Extract montonic target parameters
extract_monotonic_components <- function(posterior) {
  posterior_mat <- as.matrix(posterior)                                                       #Convert posterior draws to a matrix for easy subsetting
  list(
    logit_phi_first = posterior_mat[, grep("^logit_phi_first", colnames(posterior_mat))],     #Asocial information appetite (phi) at trial 1
    logit_phi_last = posterior_mat[, grep("^logit_phi_last", colnames(posterior_mat))],       #Asocial information appetite (phi) at final trial
    log_lambda_first = posterior_mat[, grep("^log_lambda_first", colnames(posterior_mat))],   #Asocial weighting (lambda) at trial 1
    log_lambda_last = posterior_mat[, grep("^log_lambda_last", colnames(posterior_mat))],     #Asocial weighting (lambda) at final trial
    log_epsilon_first = posterior_mat[, grep("^log_epsilon_first", colnames(posterior_mat))], #Social information appetite (epsilon) at trial 1
    log_epsilon_last = posterior_mat[, grep("^log_epsilon_last", colnames(posterior_mat))],   #Social information appetite (epsilon) at final trial
    logit_sigma_first = posterior_mat[, grep("^logit_sigma_first", colnames(posterior_mat))], #Social weighting (sigma) at trial 1
    logit_sigma_last = posterior_mat[, grep("^logit_sigma_last", colnames(posterior_mat))],   #Social weighting (sigma) at final trial
    delta_phi = posterior_mat[, grep("^delta_phi", colnames(posterior_mat))],                 #Monotonic change trajectory for phi across trials
    delta_lambda = posterior_mat[, grep("^delta_lambda", colnames(posterior_mat))],           #Monotonic change trajectory for lambda
    delta_epsilon = posterior_mat[, grep("^delta_epsilon", colnames(posterior_mat))],         #Monotonic change trajectory for epsilon
    delta_sigma = posterior_mat[, grep("^delta_sigma", colnames(posterior_mat))],             #Monotonic change trajectory for sigma
    z_ID = posterior_mat[, grep("^z_ID", colnames(posterior_mat))],                           #ID-specific latent offsets (8 parameters per individual)
    sigma_ID = posterior_mat[, grep("^sigma_ID", colnames(posterior_mat))],                   #Standard deviations of random effects
    Rho_ID = posterior_mat[, grep("^Rho_ID", colnames(posterior_mat))],                       #Correlation matrix of ID-level deviations
    A_init = posterior_mat[, grep("^A_init", colnames(posterior_mat))]                        #Initial option attraction values
  )
}

#Custom log-likelihood function for monotonic
compute_log_lik_monotonic <- function(posterior, trial_row, id = 1) {
  pars <- extract_monotonic_components(posterior)                             #Extract posterior parameter components
  S <- nrow(posterior)                                                        #Number of posterior draws 
  t <- trial_row$trial                                                        #Current trial number
  log_lik <- numeric(S)                                                       #Storage for log-likelihoods per posterior sample
  
  for (s in 1:S) {
    offset <- pars$z_ID[s, ((id - 1) * 8 + 1):(id * 8)]                       #Extract ID-specific random effect vector for 8 parameters
    Sigma_ID <- diag(as.numeric(pars$sigma_ID[s, 1:8]))                       #Build diagonal matrix of standard deviations
    Rho <- matrix(as.numeric(pars$Rho_ID[s, 1:(8 * 8)]), 8, 8)                #Reshape correlation vector into full matrix
    L <- t(chol(Sigma_ID %*% Rho %*% Sigma_ID))                               #Compute Cholesky factor of full covariance matrix
    v_ID <- L %*% offset                                                      #Transform uncorrelated deviations into correlated random effects
    
    idx <- (trial_row$manip - 1) * 2 + trial_row$sat                          #Index treatment condition 
    dL <- c(0, as.numeric(pars$delta_lambda[s, ]))                            #Cumulative delta trajectory for lambda 
    dR <- c(0, as.numeric(pars$delta_epsilon[s, ]))                           #Cumulative delta trajectory for epsilon
    dS <- c(0, as.numeric(pars$delta_sigma[s, ]))                             #Cumulative delta trajectory for sigma 
    dP <- c(0, as.numeric(pars$delta_phi[s, ]))                               #Cumulative delta trajectory for phi 
    
    lambda_f <- exp(pars$log_lambda_first[s, idx] + v_ID[1])                  #Initial lambda, transformed to positive scale
    lambda_l <- exp(pars$log_lambda_last[s, idx] + v_ID[2])                   #Final lambda
    epsilon_f <- exp(pars$log_epsilon_first[s, idx] + v_ID[3])                #Initial epsilon
    epsilon_l <- exp(pars$log_epsilon_last[s, idx] + v_ID[4])                 #Final epsilon
    sigma_f <- inv_logit(pars$logit_sigma_first[s, idx] + v_ID[5])            #Initial sigma 
    sigma_l <- inv_logit(pars$logit_sigma_last[s, idx] + v_ID[6])             #Final sigma
    phi_f <- inv_logit(pars$logit_phi_first[s, idx] + v_ID[7])                #Initial phi
    phi_l <- inv_logit(pars$logit_phi_last[s, idx] + v_ID[8])                 #Final phi
    
    lambda <- min(15, lambda_f + (lambda_l - lambda_f) * sum(dL[1:t]))        #Interpolate lambda by cumulative delta trajectory
    epsilon <- min(10, epsilon_f + (epsilon_l - epsilon_f) * sum(dR[1:t]))    #Interpolate epsilon
    sigma <- sigma_f + (sigma_l - sigma_f) * sum(dS[1:t])                     #Interpolate sigma
    phi <- phi_f + (phi_l - phi_f) * sum(dP[1:t])                             #Interpolate phi
    
    A <- inv_logit(pars$A_init[s, ((id - 1) * 2 + 1):(id * 2)])               #Initial attractions to each option
    A <- sort(A)                                                              #Ensure fixed order of options for consistency
    pA <- softmax(lambda * A)                                                 #Compute asocial choice probabilities
    
    if (trial_row$social == 0) {
      pC <- pA                                                                #Use asocial probabilities if no social info
    } else {
      safe_cum_soc <- max(trial_row$cum_soc, 1e-6)                            #Avoid division by zero
      safe_cum_non_soc <- max(trial_row$cum_non_soc, 1e-6)                    #Same safeguard
      pS <- c(
        safe_cum_soc^epsilon / (safe_cum_soc^epsilon + safe_cum_non_soc^1),
        safe_cum_non_soc^1 / (safe_cum_soc^epsilon + safe_cum_non_soc^1)
      )                                                                       #Compute social choice weights
      pC <- (1 - sigma) * pA + sigma * pS                                     #Combine asocial and social influence
    }
    
    log_lik[s] <- log(pC[trial_row$choice])                                   #Log-likelihood of chosen option
  }
  
  return(log_lik)                                                             #Return vector of log likelihoods
}

#Prepare common dataset for per-bird prediction
full_data <- data.frame(
  id = EWA_stan_list$id,                          #Bird identifier
  trial = EWA_stan_list$trial,                    #Trial number
  manip = EWA_stan_list$manip,                    #Manipulation condition (1 = Material, 2 = Breeding)
  sat = EWA_stan_list$sat,                        #Satisfaction condition (1 = Satisfied, 2 = Dissatisfied)
  treat = EWA_stan_list$treat,                    #Combined treatment label 
  choice = EWA_stan_list$choice,                  #Observed choice made by bird
  cum_soc = EWA_stan_list$cum_soc,                #Cumulative social observations up to current trial
  cum_non_soc = EWA_stan_list$cum_non_soc,        #Cumulative asocial observations up to current trial
  social = EWA_stan_list$social                   #Binary flag indicating whether social info was present
)

#Extract all bird IDs
all_ids <- unique(full_data$id)                   #Vector of unique bird IDs in dataset

#Subset birds to reduce computation, if full birds not assigned - did this initially to troubleshoot code & make sure was running/spitting out what I wanted
bird_subset <- all_ids[1:min(test_birds, length(all_ids))]   #Take up to 'test_birds' individuals

#Find maximum trial count for each bird
n_trials_per_bird <- sapply(bird_subset, function(b) max(full_data$trial[full_data$id == b]))  #Trial count per bird

#Load compiled Stan models - again, note Model_Execution.R script must have been run to correctly set working directory!
model_baseline <- cmdstan_model("Models/EWA_Baseline_Model.stan")     #CmdStan model for baseline EWA
model_monotonic <- cmdstan_model("Models/EWA_Monotonic_Model.stan")   #CmdStan model for monotonic EWA

#Create storage for elpd results
elpd_results <- data.frame(
  bird_id = bird_subset,                   #Bird ID
  elpd_baseline = NA,                      #Expected log predictive density under baseline
  elpd_monotonic = NA,                     #Expected log predictive density under monotonic
  delta_elpd = NA                          #Difference between monotonic and baseline ELPD
)

#Create storage for Pareto_k diagnostics
k_df <- data.frame()                       #Empty data frame to hold PSIS Pareto-k values per bird/trial

#Initialise refit counter per bird per model
refit_counts <- data.frame(
  bird = bird_subset,                      #Bird ID
  refits_baseline = 0,                     #Number of refits for baseline model
  refits_monotonic = 0                     #Number of refits for monotonic model
)

#
##
###Forecase for both EWA Models
##
#

#Begin loop
for (b in seq_along(bird_subset)) {
  bird_id <- bird_subset[b]                                   #Get the current bird ID
  bird_data <- subset(full_data, id == bird_id)               #Subset full data to the current bird
  N_trials <- n_trials_per_bird[b]                            #Get the number of trials for this bird
  
  #Baseline model
  stan_data <- make_stan_data_common(bird_id, bird_data[bird_data$trial <= L, ], EWA_stan_list)  #Prepare data up to horizon L
  fit <- model_baseline$sample(data = stan_data, chains = 4, seed = 123)                         #Fit the baseline model using Stan
  draws <- suppressWarnings(as_draws_df(fit$draws()))                                            #Convert Stan draws to posterior data frame
  elpd_list <- rep(NA, N_trials)                                                                 #Initialize ELPD result vector
  i_refit <- L                                                                                   #Index of last model refit
  
  for (i in (L + 1):(N_trials - M)) {
    new_data <- bird_data[bird_data$trial > i_refit & bird_data$trial <= i, ]                    #Window of trials since last refit
    k_val <- if (nrow(new_data) > 0) {
      log_ratios <- rowSums(sapply(1:nrow(new_data), function(j) compute_log_lik_baseline(draws, new_data[j, ])))  #Compute pointwise log-likelihood ratios
      psis_result <- psis(log_ratios)                                                                              #Run Pareto-smoothed importance sampling
      pareto_k_values(psis_result)                                                                                 #Extract Pareto k values
    } else 0
    
    k_df <- rbind(k_df, data.frame(bird = bird_id, trial = i, model = "baseline", k = k_val))        #Store diagnostic k value
    
    if (k_val > k_thresh) {                                                                          #If k is too large, refit the model
      refit_counts$refits_baseline[b] <- refit_counts$refits_baseline[b] + 1
      i_refit <- i
      stan_data <- make_stan_data_common(bird_id, bird_data[bird_data$trial <= i, ], EWA_stan_list)  #Update data up to trial i
      fit <- model_baseline$sample(data = stan_data, chains = 4, seed = 123)                         #Refit model
      draws <- suppressWarnings(as_draws_df(fit$draws()))                                            #Update posterior draws
    }
    
    test_row <- bird_data[bird_data$trial == i + 1, ]                                                #Get the next trial for testing
    log_lik_test <- compute_log_lik_baseline(draws, test_row)                                        #Compute test trial log-likelihood
    
    elpd_list[i + 1] <- if (k_val > k_thresh) {
      log_mean_exp(log_lik_test)                                                                     #Fallback if k was high
    } else {
      lw <- weights(psis_result, normalize = TRUE)[, 1]                                              #Normalized importance weights
      log_sum_exp(lw + log_lik_test)                                                                 #Weighted log-likelihood
    }
  }
  
  elpd_b <- sum(elpd_list, na.rm = TRUE)                                                             #Sum log predictive density values
  
  #Montonic model
  stan_data <- make_stan_data_common(bird_id, bird_data[bird_data$trial <= L, ], EWA_stan_list)      #Prepare data for monotonic model
  fit <- model_monotonic$sample(data = stan_data, chains = 4, seed = 123)                            #Fit monotonic model
  draws <- suppressWarnings(as_draws_df(fit$draws()))                                                #Convert to posterior draws
  elpd_list <- rep(NA, N_trials)                                                                     #Initialize ELPD storage
  i_refit <- L                                                                                       #Reset last refit index
  
  for (i in (L + 1):(N_trials - M)) {
    new_data <- bird_data[bird_data$trial > i_refit & bird_data$trial <= i, ]                        #Window since last refit
    k_val <- if (nrow(new_data) > 0) {
      log_ratios <- rowSums(sapply(1:nrow(new_data), function(j) compute_log_lik_monotonic(draws, new_data[j, ])))  #Compute log-likelihoods
      psis_result <- psis(log_ratios)                                                                #Run PSIS
      pareto_k_values(psis_result)                                                                   #Extract k diagnostic
    } else 0
    
    k_df <- rbind(k_df, data.frame(bird = bird_id, trial = i, model = "monotonic", k = k_val))       #Record result
    
    if (k_val > k_thresh) {                                                                          #If diagnostic is poor, refit
      refit_counts$refits_monotonic[b] <- refit_counts$refits_monotonic[b] + 1
      i_refit <- i
      stan_data <- make_stan_data_common(bird_id, bird_data[bird_data$trial <= i, ], EWA_stan_list)  #Update training window
      fit <- model_monotonic$sample(data = stan_data, chains = 4, seed = 123)                        #Refit model
      draws <- suppressWarnings(as_draws_df(fit$draws()))                                            #Extract posterior
    }
    
    test_row <- bird_data[bird_data$trial == i + 1, ]                                                #Get test trial
    log_lik_test <- compute_log_lik_monotonic(draws, test_row)                                       #Compute test log-likelihood
    
    elpd_list[i + 1] <- if (k_val > k_thresh) {
      log_mean_exp(log_lik_test)                                                                     #Fallback for unstable weights
    } else {
      lw <- weights(psis_result, normalize = TRUE)[, 1]                                              #Importance weights
      log_sum_exp(lw + log_lik_test)                                                                 #Importance-weighted log predictive
    }
  }
  
  elpd_m <- sum(elpd_list, na.rm = TRUE)                                                             #Sum up ELPD for monotonic model
  
  #Record results
  elpd_results$elpd_baseline[b] <- elpd_b                                                            #Save baseline model ELPD
  elpd_results$elpd_monotonic[b] <- elpd_m                                                           #Save monotonic model ELPD
  elpd_results$delta_elpd[b] <- elpd_m - elpd_b                                                      #Compute and save difference
}

#Results summaries
print(elpd_results)                                                                                        #Print ELPD results
print(refit_counts)                                                                                        #Print refit count table
message("Total ELPD baseline: ", round(sum(elpd_results$elpd_baseline, na.rm = TRUE), 2))                  #Summary total baseline ELPD
message("Total ELPD monotonic: ", round(sum(elpd_results$elpd_monotonic, na.rm = TRUE), 2))                #Summary total monotonic ELPD
message("Total delta ELPD (monotonic - baseline): ", round(sum(elpd_results$delta_elpd, na.rm = TRUE), 2)) #Total ELPD difference
message("Total refits (baseline): ", sum(refit_counts$refits_baseline))                                    #Total refits for baseline model
message("Total refits (monotonic): ", sum(refit_counts$refits_monotonic))                                  #Total refits for monotonic model

