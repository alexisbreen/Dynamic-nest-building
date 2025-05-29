################################################################################################################################################################################

#Script for plotting Figure S3 for the manuscript

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

#Extract samples
#Note this overwrites EWA_b and EWA_m samples, if corresponding models already run/extracted; or anything related to Figure S1 script

sb <- extract.samples(PS_EWA_b)
sm <- extract.samples(PS_EWA_m)

#Treatment group indices

treatment_indices <- list(
  c(1:10),   #Mat-Sat
  c(11:18),  #Mat-Diss
  c(19:32),  #Inc-Sat
  c(33:47)   #Inc-Diss
)

#Treatment labels and colors

param_names <- c("Phi", "Lambda", "Epsilon", "Sigma")
treatment_labels <- c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss")
treatment_colors <- c("Mat-Sat" = rangi2,
                      "Mat-Diss" = rangi2,
                      "Inc-Sat" = "gold",
                      "Inc-Diss" = "gold")

#Y-axis limits 

y_axis_limits <- list(
  Phi = c(0, 1),
  Lambda = c(0, 15),
  Epsilon = c(0, 10),
  Sigma = c(0, 1)
)

#Transform + cap functions

inv_logit_capped <- function(x) inv_logit(x)
exp_capped_lambda <- function(x) pmin(exp(x), 15)
exp_capped_epsilon <- function(x) pmin(exp(x), 10)

#
##
### BASELINE MODEL EXTRACTION
##
#

#Function to compute individual-level baseline effects

compute_base_effect <- function(param_matrix, v_ID, param_index, treat_indices, transform_fn) {
  out <- list()
  for(t in 1:4){
    treat_inds <- treatment_indices[[t]]
    g1 <- ifelse(t <= 2, 1, 2)
    g2 <- ifelse(t %% 2 == 1, 1, 2)
    temp <- list()
    for(i in treat_inds){
      draws <- transform_fn(param_matrix[, g1, g2] + v_ID[, i, param_index])
      temp[[length(temp) + 1]] <- draws
    }
    out[[t]] <- temp
  }
  return(out)
}

#Apply base function

phi_indiv_draws    <- compute_base_effect(sb$logit_phi, sb$v_ID, 4, treatment_indices, inv_logit_capped)
lambda_indiv_draws <- compute_base_effect(sb$log_lambda, sb$v_ID, 1, treatment_indices, exp_capped_lambda)
epsilon_indiv_draws    <- compute_base_effect(sb$log_epsilon, sb$v_ID, 2, treatment_indices, exp_capped_epsilon)
sigma_indiv_draws  <- compute_base_effect(sb$logit_sigma, sb$v_ID, 3, treatment_indices, inv_logit_capped)

#Bundle results for plotting

baseline_params <- list(
  phi_indiv_draws,
  lambda_indiv_draws,
  epsilon_indiv_draws,
  sigma_indiv_draws
)

#
##
### MONOTONIC MODEL EXTRACTION
##
#

#Helper function to map subject ID to (manip, sat) index

get_treatment_indices <- function(subject_id) {
  if (subject_id >= 1 & subject_id <= 10) return(c(1, 1))
  if (subject_id >= 11 & subject_id <= 18) return(c(1, 2))
  if (subject_id >= 19 & subject_id <= 32) return(c(2, 1))
  if (subject_id >= 33 & subject_id <= 47) return(c(2, 2))
  stop("Invalid subject ID")
}

#Function to compute individual-level monotonic effects

compute_mono_effect <- function(first_param, last_param, delta, v_ID_col, subset_j, transform = identity, cap = NULL, num_draws) {
  Delta_mu <- apply(delta, 2, mean)
  Delta_mu <- c(0, Delta_mu)
  
  Delta_draws <- cbind(0, delta[1:num_draws, ])
  
  effect_matrix <- matrix(0, nrow = length(subset_j), ncol = 25)
  draws_array <- array(0, dim = c(length(subset_j), num_draws, 25))
  
  for(j_idx in seq_along(subset_j)) {
    j <- subset_j[j_idx]
    idx <- get_treatment_indices(j)
    for(i in 1:25){
      first <- transform(mean(first_param[, idx[1], idx[2]]) + mean(sm$v_ID[, j, v_ID_col]))
      last  <- transform(mean(last_param[, idx[1], idx[2]]) + mean(sm$v_ID[, j, v_ID_col + 1]))
      value <- first + (last - first) * sum(Delta_mu[1:i])
      if(!is.null(cap)) value <- pmin(value, cap)
      effect_matrix[j_idx, i] <- value
    }
    for(d in 1:num_draws){
      for(i in 1:25){
        first_d <- transform(first_param[d, idx[1], idx[2]] + sm$v_ID[d, j, v_ID_col])
        last_d  <- transform(last_param[d, idx[1], idx[2]] + sm$v_ID[d, j, v_ID_col + 1])
        value_d <- first_d + (last_d - first_d) * sum(Delta_draws[d, 1:i])
        if(!is.null(cap)) value_d <- pmin(value_d, cap)
        draws_array[j_idx, d, i] <- value_d
      }
    }
  }
  
  list(effect_matrix = effect_matrix, draws = draws_array)
}

#Apply mono function

phi_results     <- lapply(treatment_indices, function(j) compute_mono_effect(sm$logit_phi_first, sm$logit_phi_last, sm$delta_phi, 7, j, inv_logit,  NULL, 10))
lambda_results  <- lapply(treatment_indices, function(j) compute_mono_effect(sm$log_lambda_first, sm$log_lambda_last, sm$delta_lambda,1, j, exp, 15, 10))
epsilon_results <- lapply(treatment_indices, function(j) compute_mono_effect(sm$log_epsilon_first, sm$log_epsilon_last, sm$delta_epsilon, 3, j, exp, 10, 10))
sigma_results   <- lapply(treatment_indices, function(j) compute_mono_effect(sm$logit_sigma_first, sm$logit_sigma_last,sm$delta_sigma, 5, j, inv_logit, NULL, 10))

#Bundle results for plotting

monotonic_params <- list(
  phi_results,
  lambda_results,
  epsilon_results,
  sigma_results
)

#
##
### PLOT FIGURE S3
##
#

pdf(file = "Figure_S3.pdf", height = 7, width = 10)  #If want pdf of plot

#General layout

par(mfrow = c(4, 4), mar = c(1, 2, 2, 1), oma = c(4, 6, 5, 1))

#Loop over the 4 learning parameters

for(p in 1:4){
  base_mat <- baseline_params[[p]]                # Extract baseline draws
  mono_list <- monotonic_params[[p]]              # Extract monotonic results
  ylims <- y_axis_limits[[param_names[p]]]        # Set y-axis limits based on parameter
  
  treat_indices <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))  # Treatment cell ordering
  
  #Loop through the 4 treatment groups
  
  for(j in 1:4){
    idx <- treat_indices[[j]]
    indiv_draws_list <- base_mat[[j]]                    # Individual posterior draws for this treatment
    traj_mat <- mono_list[[j]]$effect_matrix             # Mean learning trajectories per individual
    
    #Set up an empty plot area
    
    plot(NULL, xlim = c(0.5, 34.5), ylim = ylims, xaxt = "n", yaxt = "n",
         xlab = "", ylab = ifelse(j == 1, param_names[p], ""))
    
    #Add treatment titles to the top row
    
    if (p == 1) {
      mtext(treatment_labels[j], side = 3, line = 0.5, cex = 1)
    }
    
    #Add y-axis for each panel
    
    axis(2, at = pretty(ylims, n = 5), las = 1, cex.axis = 0.8)
    
    #Compute mean + HPDI of individual-level posteriors
    
    indiv_means <- sapply(indiv_draws_list, mean)
    indiv_hpdi <- lapply(indiv_draws_list, HPDI, prob = 0.89)
    
    #Plot 100 jittered posterior draws as semi-transparent points
    
    jittered_x <- jitter(rep(3.5, length(indiv_means)), amount = 2.5)
    pch_type <- ifelse(j %% 2 == 1, 1, 2)
    
    #Add individual-level HPDI segments and posterior means
    
    for (i in seq_along(indiv_means)) {
      segments(jittered_x[i], indiv_hpdi[[i]][1], jittered_x[i], indiv_hpdi[[i]][2],
               col = col.alpha("black", 0.2), lwd = 1.5)
    }
    for (i in seq_along(indiv_means)) {
      points(jittered_x[i], indiv_means[i],
             col = col.alpha(treatment_colors[j], 0.8), pch = pch_type)
    }
    
    #Add vertical line to separate baseline and monotonic estimates
    
    abline(v = 8, lty = 2, col = "black")
    
    #Plot each individual's mean posterior learning trajectories (faint)
    
    for (subj in 1:nrow(traj_mat)) {
      lines(10:34, traj_mat[subj, ],
            col = col.alpha(treatment_colors[j], 0.4), lwd = 1)
    }
    
    # No overall mean line (removed)
    
    # Add horizontal reference line for epsilon plots
    #if (p == 3) {
    #abline(h = 1, lty = 2, col = "black")
    #}
    
    #Add custom x-axis labels
    
    axis(1, at = 3.5, labels = "1-25", cex.axis = 0.7, tick = TRUE, line = 0)
    axis(1, at = c(10, 15, 20, 25, 30, 34), labels = c("1", "5", "10", "15", "20", "25"), cex.axis = 0.8)
    
    #Add left-side y-axis labels only to leftmost column
    
    if (j == 1) {
      label_text <- switch(p,
                           expression(paste("updating ", italic(phi))),
                           expression(paste("weighting ", italic(lambda))),
                           expression(paste("updating ", italic(epsilon))),
                           expression(paste("weighting ", italic(sigma))))
      mtext(label_text, side = 2, line = 2.5, cex = 1)
    }
  }
}

#Add common outer labels

mtext("Choice", side = 1, line = 1.5, outer = TRUE, cex = 1)
mtext("Individual cognition", side = 3, line = 0.5, outer = TRUE, font = 2, cex = 1)

#Adjust vertical placement to align text with row centers

mtext("Asocial information", side = 2, line = 3, outer = TRUE, at = 0.74, cex = 1)
mtext("Social information",  side = 2, line = 3, outer = TRUE, at = 0.24, cex = 1)

dev.off() #Close PDF file
