################################################################################################################################################################################

#Script for getting EWA baseline and monotonic model summaries for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Load required libraries

library(rethinking)
library(tibble)

#Extract samples

sb <- extract.samples(EWA_b)
sm <- extract.samples(EWA_m)

#Define treatment indices and labels
treatment_indices <- list(
  c(1:10),   #Mat-Sat
  c(11:18),  #Mat-Diss
  c(19:32),  #Inc-Sat
  c(33:47)   #Inc-Diss
)

treatment_labels <- c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss")

#Transform + cap functions
inv_logit_capped <- function(x) inv_logit(x)
exp_capped_lambda <- function(x) pmin(exp(x), 15)
exp_capped_epsilon <- function(x) pmin(exp(x), 10)

#
##
### BASELINE ESTIMATES
##
#

#Function to compute posterior summaries per treatment
summarise_baseline <- function(param_matrix, v_ID, param_index, transform_fn) {
  out <- list()
  for (t in 1:4) {
    inds <- treatment_indices[[t]]
    g1 <- ifelse(t <= 2, 1, 2)
    g2 <- ifelse(t %% 2 == 1, 1, 2)
    all_draws <- sapply(inds, function(i) {
      transform_fn(param_matrix[, g1, g2] + v_ID[, i, param_index])
    })
    # Take the mean across individuals (row-wise)
    pooled <- rowMeans(all_draws)
    out[[t]] <- c(mean = mean(pooled), HPDI(pooled, prob = 0.89))
  }
  return(out)
}

#Compute summaries
baseline_phi    <- summarise_baseline(sb$logit_phi, sb$v_ID, 4, inv_logit_capped)
baseline_lambda <- summarise_baseline(sb$log_lambda, sb$v_ID, 1, exp_capped_lambda)
baseline_epsilon    <- summarise_baseline(sb$log_epsilon, sb$v_ID, 2, exp_capped_epsilon)
baseline_sigma  <- summarise_baseline(sb$logit_sigma, sb$v_ID, 3, inv_logit_capped)

#
##
### MONOTONIC ESTIMATES
##
#

summarise_monotonic <- function(first_param, last_param, v_ID, index, transform_fn) {
  out <- list()
  for (t in 1:4) {
    inds <- treatment_indices[[t]]
    g1 <- ifelse(t <= 2, 1, 2)
    g2 <- ifelse(t %% 2 == 1, 1, 2)
    
    all_first <- sapply(inds, function(i) {
      transform_fn(first_param[, g1, g2] + v_ID[, i, index])
    })
    all_last <- sapply(inds, function(i) {
      transform_fn(last_param[, g1, g2] + v_ID[, i, index + 1])
    })
    
    pooled_first <- rowMeans(all_first)
    pooled_last  <- rowMeans(all_last)
    
    out[[t]] <- list(
      first = c(mean = mean(pooled_first), HPDI(pooled_first, prob = 0.89)),
      last  = c(mean = mean(pooled_last),  HPDI(pooled_last,  prob = 0.89))
    )
  }
  return(out)
}

#Apply to each parameter
mono_phi    <- summarise_monotonic(sm$logit_phi_first, sm$logit_phi_last,   sm$v_ID, 7, inv_logit)
mono_lambda <- summarise_monotonic(sm$log_lambda_first, sm$log_lambda_last,  sm$v_ID, 1, function(x) pmin(exp(x), 15))
mono_epsilon    <- summarise_monotonic(sm$log_epsilon_first, sm$log_epsilon_last, sm$v_ID, 3, function(x) pmin(exp(x), 10))
mono_sigma  <- summarise_monotonic(sm$logit_sigma_first, sm$logit_sigma_last, sm$v_ID, 5, inv_logit)

#
##
###COMBINE INTO FINAL TABLE
##
#

make_summary_table <- function(param_name, baseline, mono) {
  df <- tibble(
    Treatment = treatment_labels,
    Parameter = param_name,
    Baseline_Mean = sapply(baseline, `[[`, "mean"),
    Baseline_Lo = sapply(baseline, `[[`, 2),
    Baseline_Hi = sapply(baseline, `[[`, 3),
    First_Mean = sapply(mono, function(x) x$first["mean"]),
    First_Lo = sapply(mono, function(x) x$first[2]),
    First_Hi = sapply(mono, function(x) x$first[3]),
    Last_Mean = sapply(mono, function(x) x$last["mean"]),
    Last_Lo = sapply(mono, function(x) x$last[2]),
    Last_Hi = sapply(mono, function(x) x$last[3])
  )
  return(df)
}

phi_table    <- make_summary_table("phi", baseline_phi, mono_phi)
lambda_table <- make_summary_table("lambda", baseline_lambda, mono_lambda)
epsilon_table    <- make_summary_table("epsilon", baseline_epsilon, mono_epsilon)
sigma_table  <- make_summary_table("sigma", baseline_sigma, mono_sigma)

#Combine all
final_table <- rbind(phi_table, lambda_table, epsilon_table, sigma_table)

#View table
print(final_table)

