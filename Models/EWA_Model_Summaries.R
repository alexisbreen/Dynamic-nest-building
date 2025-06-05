################################################################################################################################################################################

#Script for getting EWA baseline and monotonic model summaries for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

##
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

#
##
### BASELINE ESTIMATES
##
#

#Function to compute posterior summaries per treatment
summarise_baseline <- function(param_matrix, v_ID, param_index, transform_fn) {
  out <- list()
  for(t in 1:4){
    inds <- treatment_indices[[t]]
    g1 <- ifelse(t <= 2, 1, 2)
    g2 <- ifelse(t %% 2 == 1, 1, 2)
    all_draws <- sapply(inds, function(i){
      transform_fn(param_matrix[, g1, g2] + v_ID[, i, param_index])
    })
    pooled <- rowMeans(all_draws)
    out[[t]] <- c(mean = mean(pooled), HPDI(pooled, prob = 0.89))
  }
  return(out)
}

baseline_phi     <- summarise_baseline(sb$logit_phi, sb$v_ID, 4, inv_logit)
baseline_lambda  <- summarise_baseline(sb$log_lambda, sb$v_ID, 1, function(x) pmin(exp(x), 15))
baseline_epsilon <- summarise_baseline(sb$log_epsilon, sb$v_ID, 2, function(x) pmin(exp(x), 10))
baseline_sigma   <- summarise_baseline(sb$logit_sigma, sb$v_ID, 3, inv_logit)

#
##
### MONOTONIC ESTIMATES
##
#

summarise_monotonic <- function(first_param, last_param, delta, v_ID_col, subset_j, transform = identity, cap = NULL) {
  Delta_mu <- apply(delta, 2, mean)
  Delta_mu <- c(0, Delta_mu)
  
  effect_matrix <- matrix(0, nrow = length(subset_j), ncol = 25)
  
  for(j_idx in seq_along(subset_j)){
    j <- subset_j[j_idx]
    idx <- get_treatment_indices(j)
    for (i in 1:25) {
      first <- transform(mean(first_param[, idx[1], idx[2]]) + mean(sm$v_ID[, j, v_ID_col]))
      last  <- transform(mean(last_param[, idx[1], idx[2]]) + mean(sm$v_ID[, j, v_ID_col + 1]))
      value <- first + (last - first) * sum(Delta_mu[1:i])
      if (!is.null(cap)) value <- pmin(value, cap)
      effect_matrix[j_idx, i] <- value
    }
  }
  
  first_values <- effect_matrix[, 1]
  last_values  <- effect_matrix[, 25]
  list(
    first = c(mean = mean(first_values), HPDI(first_values, prob = 0.89)),
    last  = c(mean = mean(last_values),  HPDI(last_values,  prob = 0.89))
  )
}

mono_phi     <- lapply(treatment_indices, function(j) summarise_monotonic(sm$logit_phi_first, sm$logit_phi_last, sm$delta_phi, 7, j, inv_logit))
mono_lambda  <- lapply(treatment_indices, function(j) summarise_monotonic(sm$log_lambda_first, sm$log_lambda_last, sm$delta_lambda, 1, j, exp_capped_lambda, 15))
mono_epsilon <- lapply(treatment_indices, function(j) summarise_monotonic(sm$log_epsilon_first, sm$log_epsilon_last, sm$delta_epsilon, 3, j, exp_capped_epsilon, 10))
mono_sigma   <- lapply(treatment_indices, function(j) summarise_monotonic(sm$logit_sigma_first, sm$logit_sigma_last, sm$delta_sigma, 5, j, inv_logit))

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

#Combine all & round

final_table_rounded <- final_table
num_cols <- sapply(final_table_rounded, is.numeric)
final_table_rounded[num_cols] <- lapply(final_table_rounded[num_cols], function(x) round(x, 2))

#View table

print(final_table_rounded)

