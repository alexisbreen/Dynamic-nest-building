################################################################################################################################################################################

#Script for running EWA monotonic model-generated forward simulations for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) & Richard McElreath (richard_mcelreath@eva.mpg.de)

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#First derive the full posterior distribution with random effects variance-covariance matrix 

#Extract posterior from EWA model

sm <- extract.samples(EWA_m)

#1) Get standard deviation of learning parameters (variance) among individuals

sigma_id_m <- sm$sigma_ID

#2) Compute correlation matrix from Cholesky factors

Correlations_m <- array(numeric(),c(4000,8,8))
for (i in 1:4000) {
  Correlations_m[i,,] <- sm$Rho_ID[i,,] %*% t(sm$Rho_ID[i,,])
}

#3) Compute variance-covariance matrix from variances and correlations

S_m <- array(numeric(),c(4000,8,8))
for (i in 1:4000) {
  S_m[i,,] <- diag(sigma_id_m[i,]) %*% Correlations_m[i,,]  %*% diag(sigma_id_m[i,])
}

#Extract posteriors for means

#Phi

phi_first <- array(numeric(), c(4000,2,2))
phi_first[,1,1] <- sm$logit_phi_first[ ,1,1] #Mat-Sat
phi_first[,1,2] <- sm$logit_phi_first[ ,1,2] #Mat-Dis
phi_first[,2,1] <- sm$logit_phi_first[ ,2,1] #Inc-Sat
phi_first[,2,2] <- sm$logit_phi_first[ ,2,2] #Inc-Diss

phi_last <- array(numeric(), c(4000,2,2))
phi_last[,1,1] <- sm$logit_phi_last[ ,1,1] #Mat-Sat
phi_last[,1,2] <- sm$logit_phi_last[ ,1,2] #Mat-Dis
phi_last[,2,1] <- sm$logit_phi_last[ ,2,1] #Inc-Sat
phi_last[,2,2] <- sm$logit_phi_last[ ,2,2] #Inc-Diss

#Lambda

lambda_first <- array(numeric(), c(4000,2,2))
lambda_first[,1,1] <- sm$log_lambda_first[ ,1,1] #Mat-Sat
lambda_first[,1,2] <- sm$log_lambda_first[ ,1,2] #Mat-Dis
lambda_first[,2,1] <- sm$log_lambda_first[ ,2,1] #Inc-Sat
lambda_first[,2,2] <- sm$log_lambda_first[ ,2,2] #Inc-Diss

lambda_last <- array(numeric(), c(4000,2,2))
lambda_last[,1,1] <- sm$log_lambda_last[ ,1,1] #Mat-Sat
lambda_last[,1,2] <- sm$log_lambda_last[ ,1,2] #Mat-Dis
lambda_last[,2,1] <- sm$log_lambda_last[ ,2,1] #Inc-Sat
lambda_last[,2,2] <- sm$log_lambda_last[ ,2,2] #Inc-Diss

#Epsilon

epsilon_first <- array(numeric(), c(4000,2,2))
epsilon_first[,1,1] <- sm$log_epsilon_first[ ,1,1] #Mat-Sat
epsilon_first[,1,2] <- sm$log_epsilon_first[ ,1,2] #Mat-Dis
epsilon_first[,2,1] <- sm$log_epsilon_first[ ,2,1] #Inc-Sat
epsilon_first[,2,2] <- sm$log_epsilon_first[ ,2,2] #Inc-Diss

epsilon_last <- array(numeric(), c(4000,2,2))
epsilon_last[,1,1] <- sm$log_epsilon_last[ ,1,1] #Mat-Sat
epsilon_last[,1,2] <- sm$log_epsilon_last[ ,1,2] #Mat-Dis
epsilon_last[,2,1] <- sm$log_epsilon_last[ ,2,1] #Inc-Sat
epsilon_last[,2,2] <- sm$log_epsilon_last[ ,2,2] #Inc-Diss

#Sigma
sigma_first <- array(numeric(), c(4000,2,2))
sigma_first[,1,1] <- sm$logit_sigma_first[ ,1,1] #Mat-Sat
sigma_first[,1,2] <- sm$logit_sigma_first[ ,1,2] #Mat-Dis
sigma_first[,2,1] <- sm$logit_sigma_first[ ,2,1] #Inc-Sat
sigma_first[,2,2] <- sm$logit_sigma_first[ ,2,2] #Inc-Diss

sigma_last <- array(numeric(), c(4000,2,2))
sigma_last[,1,1] <- sm$logit_sigma_last[ ,1,1] #Mat-Sat
sigma_last[,1,2] <- sm$logit_sigma_last[ ,1,2] #Mat-Dis
sigma_last[,2,1] <- sm$logit_sigma_last[ ,2,1] #Inc-Sat
sigma_last[,2,2] <- sm$logit_sigma_last[ ,2,2] #Inc-Diss

#Build post-hoc forward simulation function

PH_sim_mono_fct <- function( #Begin function definition and parameter imputation
  N_sim = 1,            #Number of simulations to run
  N_choices = 25,       #Number of nest-material choices
  N_choosers = 500,     #Total number of choosers to draw from posterior - must be divisible by four! 
  N_soc_pay = 1,        #Reward payoff for choosing social material
  N_non_soc_pay = 0,    #Reward payoff for choosing non-social material
  SS_matched = 1)       #Whether to simulate within-treatment numbers matched to true sample size or use N_choosers value (0 = use N_choosers; 1 = use matched)
  
{ #Begin definition of function operations
  
  #Create output object i.e., final dataset
  
  d_Overall <- c()
  
  #Begin looping over simulations
  
  for(sim in 1:N_sim){
    
    if(SS_matched == 0){
      if(N_choosers %% 4 != 0){
        print("Simulation stopped because N_choosers must be an integer divisible by four e.g., 100 and not 101")
        stop()
      }
    }
    
    #For generating choosers within each treatment, decide whether to use matched to true sample size or not
    
    if(SS_matched == 0){ #If not, use defined number of simulated choosers divided by 4
      N_S1 <- (N_choosers/4)
      N_D1 <- (N_choosers/4)
      N_S2 <- (N_choosers/4)
      N_D2 <- (N_choosers/4)
    } else { #If yes, use true sample size
      N_S1 <- 10
      N_D1 <- 8
      N_S2 <- 14
      N_D2 <- 15
    } #End sample size determination
    
    #Define total number of individual choosers based on sample size determination above
    
    N_id <- N_S1 + N_D1 + N_S2 + N_D2
    
    #Define vector of random draws to be pulled from the posterior for each learning parameter 
    #If uneven treatment numbers, draws need this sample size information
    
    if(N_S1 != (N_choosers/4) | N_D1 != (N_choosers/4) | N_S2 != (N_choosers/4) | N_D2 != (N_choosers/4)){
      draw_S1 <- sample(1:4000, N_S1, replace = TRUE)
      draw_D1 <- sample(1:4000, N_D1, replace = TRUE)
      draw_S2 <- sample(1:4000, N_S2, replace = TRUE)
      draw_D2 <- sample(1:4000, N_D2, replace = TRUE)
    } else { #If even treatment numbers, draws can be the total number of simulated choosers divided by 4 
      draw <- sample(1:4000, (N_choosers/4), replace = TRUE)
      draw_S1 <-  draw
      draw_D1 <-  draw
      draw_S2 <-  draw
      draw_D2 <-  draw
    }
    
    #Generate choosers from multivariate normal distribution
    
    Sim_S1 <- t(sapply(1:N_S1, function(i) rmvnorm(1, c(phi_first[draw_S1[i], 1, 1], phi_last[draw_S1[i], 1, 1], 
                                                        lambda_first[draw_S1[i], 1, 1], lambda_last[draw_S1[i], 1, 1], 
                                                        epsilon_first[draw_S1[i], 1, 1], epsilon_last[draw_S1[i], 1, 1], 
                                                        sigma_first[draw_S1[i], 1, 1], sigma_last[draw_S1[i], 1, 1]), 
                                                        S_m[draw_S1[i],,]))) #Mat-Sat
    Sim_D1 <- t(sapply(1:N_D1, function(i) rmvnorm(1, c(phi_first[draw_D1[i],1, 2], phi_last[draw_D1[i], 1, 2], 
                                                        lambda_first[draw_D1[i], 1, 2], lambda_last[draw_D1[i], 1, 2], 
                                                        epsilon_first[draw_D1[i], 1, 2], epsilon_last[draw_D1[i], 1, 2], 
                                                        sigma_first[draw_D1[i], 1, 2], sigma_last[draw_D1[i], 1, 2]), 
                                                        S_m[draw_D1[i],,]))) #Mat-Dis
    Sim_S2 <- t(sapply(1:N_S2, function(i) rmvnorm(1, c(phi_first[draw_S2[i],2, 1], phi_last[draw_S2[i], 2, 1], 
                                                        lambda_first[draw_S2[i], 2, 1], lambda_last[draw_S2[i], 2, 1], 
                                                        epsilon_first[draw_S2[i], 2, 1], epsilon_last[draw_S2[i], 2, 1], 
                                                        sigma_first[draw_S2[i], 2, 1], sigma_last[draw_S2[i], 2, 1]), 
                                                        S_m[draw_S2[i],,]))) #Inc-Sat
    Sim_D2 <- t(sapply(1:N_D2, function(i) rmvnorm(1, c(phi_first[draw_D2[i],2, 2], phi_last[draw_D2[i], 2, 2], 
                                                        lambda_first[draw_D2[i], 2, 2], lambda_last[draw_D2[i], 2, 2], 
                                                        epsilon_first[draw_D2[i], 2, 2], epsilon_last[draw_D2[i], 2, 2], 
                                                        sigma_first[draw_D2[i], 2, 2], sigma_last[draw_D2[i], 2, 2]), 
                                                        S_m[draw_D2[i],,]))) #Inc-Sat
    
    #Combine draws
    
    Sim_all <- rbind(Sim_S1, Sim_D1, Sim_S2, Sim_D2)
    
    #Index draws to target parameters & transform back to original scale to be used later in the model
    
    phi_f <- inv_logit(Sim_all[, 1])
    phi_l <- inv_logit(Sim_all[, 2])
    lambda_f <- exp(Sim_all[, 3])
    lambda_l <- exp(Sim_all[, 4])
    epsilon_f <- exp(Sim_all[, 5])
    epsilon_l <- exp(Sim_all[, 6])
    sigma_f <- inv_logit(Sim_all[, 7])
    sigma_l <- inv_logit(Sim_all[, 8])
    
    #Deltas
    
    #Phi
    
    Delta_phi_S1 <- matrix(0, length(draw_S1), 24)
    Delta_phi_D1 <- matrix(0, length(draw_D1), 24)
    Delta_phi_S2 <- matrix(0, length(draw_S2), 24)
    Delta_phi_D2 <- matrix(0, length(draw_D2), 24)
    
    for(i in 1:24){
      for(j in 1:length(draw_S1))
        Delta_phi_S1[j, i] <- sm$delta_phi[draw_S1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D1))
        Delta_phi_D1[j, i] <- sm$delta_phi[draw_D1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_S2))
        Delta_phi_S2[j, i] <- sm$delta_phi[draw_S2[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D2))
        Delta_phi_D2[j, i] <- sm$delta_phi[draw_D2[j],i]
    }
    
    Delta_phi_S1 <- cbind(0, Delta_phi_S1)
    Delta_phi_D1 <- cbind(0, Delta_phi_D1)
    Delta_phi_S2 <- cbind(0, Delta_phi_S2)
    Delta_phi_D2 <- cbind(0, Delta_phi_D2)
    
    Delta_phi <- rbind(Delta_phi_S1,Delta_phi_D1,Delta_phi_S2,Delta_phi_D2)
    
    # Lambda
    
    Delta_lambda_S1 <- matrix(0, length(draw_S1), 24)
    Delta_lambda_D1 <- matrix(0, length(draw_D1), 24)
    Delta_lambda_S2 <- matrix(0, length(draw_S2), 24)
    Delta_lambda_D2 <- matrix(0, length(draw_D2), 24)
    
    for(i in 1:24){
      for(j in 1:length(draw_S1))
      Delta_lambda_S1[j, i] <- sm$delta_lambda[draw_S1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D1))
        Delta_lambda_D1[j, i] <- sm$delta_lambda[draw_D1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_S2))
        Delta_lambda_S2[j, i] <- sm$delta_lambda[draw_S2[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D2))
        Delta_lambda_D2[j, i] <- sm$delta_lambda[draw_D2[j],i]
    }
    
    Delta_lambda_S1 <- cbind(0, Delta_lambda_S1)
    Delta_lambda_D1 <- cbind(0, Delta_lambda_D1)
    Delta_lambda_S2 <- cbind(0, Delta_lambda_S2)
    Delta_lambda_D2 <- cbind(0, Delta_lambda_D2)
    
    Delta_lambda <- rbind(Delta_lambda_S1,Delta_lambda_D1,Delta_lambda_S2,Delta_lambda_D2)
    
    #Epsilon
    
    Delta_epsilon_S1 <- matrix(0, length(draw_S1), 24)
    Delta_epsilon_D1 <- matrix(0, length(draw_D1), 24)
    Delta_epsilon_S2 <- matrix(0, length(draw_S2), 24)
    Delta_epsilon_D2 <- matrix(0, length(draw_D2), 24)
    
    for(i in 1:24){
      for(j in 1:length(draw_S1))
        Delta_epsilon_S1[j, i] <- sm$delta_epsilon[draw_S1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D1))
        Delta_epsilon_D1[j, i] <- sm$delta_epsilon[draw_D1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_S2))
        Delta_epsilon_S2[j, i] <- sm$delta_epsilon[draw_S2[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D2))
        Delta_epsilon_D2[j, i] <- sm$delta_epsilon[draw_D2[j],i]
    }
    
    Delta_epsilon_S1 <- cbind(0, Delta_epsilon_S1)
    Delta_epsilon_D1 <- cbind(0, Delta_epsilon_D1)
    Delta_epsilon_S2 <- cbind(0, Delta_epsilon_S2)
    Delta_epsilon_D2 <- cbind(0, Delta_epsilon_D2)
    
    Delta_epsilon <- rbind(Delta_epsilon_S1,Delta_epsilon_D1,Delta_epsilon_S2,Delta_epsilon_D2)

    #Sigma
    
    Delta_sigma_S1 <- matrix(0, length(draw_S1), 24)
    Delta_sigma_D1 <- matrix(0, length(draw_D1), 24)
    Delta_sigma_S2 <- matrix(0, length(draw_S2), 24)
    Delta_sigma_D2 <- matrix(0, length(draw_D2), 24)
    
    for(i in 1:24){
      for(j in 1:length(draw_S1))
        Delta_sigma_S1[j, i] <- sm$delta_sigma[draw_S1[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D1))
        Delta_sigma_D1[j, i] <- sm$delta_sigma[draw_D1[j],i]
    }
    for(i in 1:24) {
      for(j in 1:length(draw_S2))
        Delta_sigma_S2[j, i] <- sm$delta_sigma[draw_S2[j],i]
    }
    for(i in 1:24){
      for(j in 1:length(draw_D2))
        Delta_sigma_D2[j, i] <- sm$delta_sigma[draw_D2[j],i]
    }
    
    Delta_sigma_S1 <- cbind(0, Delta_sigma_S1)
    Delta_sigma_D1 <- cbind(0, Delta_sigma_D1)
    Delta_sigma_S2 <- cbind(0, Delta_sigma_S2)
    Delta_sigma_D2 <- cbind(0, Delta_sigma_D2)
    
    Delta_sigma <- rbind(Delta_sigma_S1,Delta_sigma_D1,Delta_sigma_S2,Delta_sigma_D2)

    #Empty attraction matrix to hold initial attraction draws of social and non-social material
    
    Atx <- matrix(0, N_id, 2)
    
    #When making attraction draws, need to specify from which birds the estimates are to be pulled
    #This is because the attraction estimates are indexed in stan by choice and bird id
    #And if the forward simulation is not matched to study sample size, the attraction draws return 'subscript out of bounds' - for example, when trying to find bird 48, when we only have a maximum of 47 bird attractions estimated from true sample
    #So, we decide how attractions are to be drawn based on whether we're matched (or not) to study sample size
    
    #Combine our previously defined treatment draws into one vector, so that our for loops below, which loop through birds within treatments, draw attractions from the correct bird and treatment, because attractions in stan are indexed by choice and bird id
    
    draw_all <- c(draw_S1,draw_D1,draw_S2,draw_D2) 
    
    if(SS_matched == 1){ #If matched...
      
      #Draw attractions
      #We draw from estimated attraction scores at deposit 1 (versus estimated baseline attraction before deposit 1)
      #We do this because these draws will be closer to true attractions due to recursive nature of estimation 
      
      for(i in 1:N_S1){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, i], sm$Atx_non_soc[draw_all[i], 1, i])} #Satisfied-construction
      for(i in (max(N_S1) + 1):(max(N_S1) + max(N_D1))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, i], sm$Atx_non_soc[draw_all[i], 1, i])} #Dissatisfied-construction
      for(i in (max(N_S1) + max(N_D1) + 1):(max(N_S1) + max(N_D1) + max(N_S2))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, i], sm$Atx_non_soc[draw_all[i], 1, i])} #Satisfied-reproduction
      for(i in (max(N_S1) + max(N_D1) + max(N_S2) + 1):(max(N_S1) + max(N_D1) + max(N_S2) + max(N_D2))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, i], sm$Atx_non_soc[draw_all[i], 1, i])} #Dissatisfied-reproduction
      
    } else { #If not matched
      
      draw_S1_A <- sample(1:10, N_S1, replace = TRUE)
      draw_D1_A <- sample(11:18, N_D1, replace = TRUE)
      draw_S2_A <- sample(19:32, N_S2, replace = TRUE)
      draw_D2_A <- sample(33:47, N_D2, replace = TRUE)
      
      id_A <- c(draw_S1_A,draw_D1_A,draw_S2_A,draw_D2_A) 
      
      #Draw attractions
      #We draw from estimated attraction scores at deposit 1 (versus estimated baseline attraction before deposit 1)
      #We do this because these draws will be closer to true attractions due to recursive nature of estimation -
      
      for(i in 1:N_S1){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, id_A[i]], sm$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Satisfied-construction
      for(i in (max(N_S1) + 1):(max(N_S1) + max(N_D1))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, id_A[i]], sm$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Dissatisfied-construction
      for(i in (max(N_S1) + max(N_D1) + 1):(max(N_S1) + max(N_D1) + max(N_S2))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, id_A[i]], sm$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Satisfied-reproduction
      for(i in (max(N_S1) + max(N_D1) + max(N_S2) + 1):(max(N_S1) + max(N_D1) + max(N_S2) + max(N_D2))){Atx[i, ] <- c(sm$Atx_soc[draw_all[i], 1, id_A[i]], sm$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Dissatisfied-reproduction
      
    }
    
    #Define vectors specifying experiment and satisfaction level
    
    if(SS_matched == 0){
      exp_id <- rep(1:2, each = (N_choosers/2))
      sat_id <- rep(1:2, each = (N_choosers/4), times = 2)
    } else {
      exp_id <- rep(1:2, c(18,29))
      sat_id <- rep(c(1,2,1,2), c(10,8,14,15))
    }
    
    #Loop over choosers
    
    for(ind in 1:N_id){ 
      
      #print(ind) #To live-track simulation progress 
      
      #Set initial attractions to coloured nest material
      
      A <- c(Atx[ind,1], Atx[ind,2]) #(social option, non-social option)
      
      #Create output dataframe, which we want to hold:
      
      d <- data.frame(id = ind,
                      trial = 1:N_choices, 
                      experiment = exp_id[ind],
                      sat_level = sat_id[ind],
                      treat = NA,
                      choice = NA, 
                      social = NA,
                      copy = NA,
                      cum_soc = NA,
                      Atx_soc = NA,
                      sim = sim)
      
      #Empty matrix to fill with nestbox counts of each material-type; used when soc = 1 for probability determination
      
      NB <- matrix(0, 1, 2) 
      soc <- 0  
      
      #Loop over material-choice decisions i.e., trials
      
      for(trial in 1:N_choices){
        
        #First, before each choice, determine whether, and, if so, how much of the two material-types are within the nestbox
        
        if(trial > 1){ #Cannot have any material in nestbox if trial 1
          NB[1] <- length(which(d$copy[- trial] == 1)) #How many copy responses are from current trial to the previous trial
          NB[2] <- length(which(d$copy[- trial] == 0)) #How many no-copy responses are from current trial to the previous trial
        }
        
        #Decide if nestbox contains at least one social material, and, if so, turn social vector 'on'
        
        if(sum(NB[1]) >= 1) soc <- 1
        d$social[trial] <- soc
        
        #Lambda
        
        lambda <- lambda_f[ind] + (lambda_l[ind] - lambda_f[ind]) * sum(Delta_lambda[ind, 1:trial])
        lambda <- ifelse(lambda > 15, 15, lambda)  #Cap
        
        #Calculate decision probabilities 
        
        Prob <- c() #Empty vector to fill with decision probabilities
        if(soc == 0){
          Prob[1] <- exp(lambda * A[1]) / sum(exp(lambda * A)) 
          Prob[2] <- exp(lambda * A[2]) / sum(exp(lambda * A)) 
        } else {
          
          #Epsilon
          
          epsilon <- epsilon_f[ind] + (epsilon_l[ind] - epsilon_f[ind]) * sum(Delta_epsilon[ind, 1:trial])
          epsilon <- ifelse(epsilon > 10, 10, epsilon) #Cap
          
          #Sigma
          
          sigma <- sigma_f[ind] + (sigma_l[ind] - sigma_f[ind]) * sum(Delta_sigma[ind, 1:trial])
          
          #Prob
          
          Prob[1] <- (1 - sigma) * exp(lambda * A[1]) / sum(exp(lambda * A)) + (sigma * ((NB[1] ^ epsilon) / ( (NB[1] ^ epsilon) + (NB[2] ^ (1)) ))) 
          Prob[2] <- (1 - sigma) * exp(lambda * A[2]) / sum(exp(lambda * A)) + (sigma * ((NB[2] ^ (1)) / ( (NB[1] ^ epsilon) + (NB[2] ^ (1)) ))) 
          
        }
        
        #Make choice proportional to attraction scores
        
        d$choice[which(d$trial == trial)] <- sample(c(1:2), size = 1, prob = Prob) #Make the bird choose among the two choice-options, based on the choice probabilities defined above    
        
        #Determine copy based on choice
        
        if(d$choice[trial] == 1){ #If the bird chose the social option
          d$copy[trial] <- 1      #The bird is assigned to have copied 
        } else {                  #But if the bird chose the non-social option
          d$copy[trial] <- 0      #The bird is assigned to have not copied 
        }                         #End determination 
        
        #Phi
       
        phi <- phi_f[ind] + (phi_l[ind] - phi_f[ind]) * sum(Delta_phi[ind, 1:trial])
        
        #Update attractions
        
        if(d$choice[trial] == 1){                                                #If choose social
          A[1] <- (1 - phi) * A[d$choice[trial]] + phi * N_soc_pay     #Attraction for social increase
          A[2] <- (1 - phi) * A[d$choice[trial]] + phi * 0             #Attraction for non-social decrease
        } else {                                                                 #If choose non-social
          A[1] <- (1 - phi) * A[d$choice[trial]] + phi * 0             #Attraction for social decrease
          A[2] <- (1 - phi) * A[d$choice[trial]] + phi * N_non_soc_pay #Attraction for non-social increase
        }
        
        #Assign attraction scores to data frame 
        
        d$Atx_soc[trial] <- A[1]     #This is social
        d$Atx_non_soc[trial] <- A[2] #This is non-social
        
        #Assign treatment to data frame
        
        if(d$experiment[trial] == 1){  #If experimental treatment is construction
          if(d$sat_level[trial] == 1){ #And if satisfaction level is satisfied
            d$treat[trial] <- 1        #Treatment is 1 i.e., satisfied-construction
          } else {                     #But if satisfaction level is dissatisfied
            d$treat[trial] <- 2        #Treatment is 2 i.e., dissatisfied-construction
          }                            #End determination for treatment 1 & 2
        } else {                       #If experimental treatment is not construction but reproduction
          if(d$sat_level[trial] == 1){ #And if satisfaction level is satisfied
            d$treat[trial] <- 3        #Treatment is 3 i.e., Satisfied-reproduction
          } else {                     #But if satisfaction level is dissatisfied
            d$treat[trial] <- 4        #Treatment is 4 i.e., Dissatisfied-reproduction
          }                            #End determination for treatment 3 & 4
        }                              #End determination for all treatments
        
        #Also, assign cumulative social material count based on nestbox
        
        if(trial == 1){
          d$cum_soc[trial] <- d$copy[trial]
        } else {
          d$cum_soc[trial] <- length(which(d$copy[1:trial] == 1))
          
        }
        
      } #End looping over material-choice decisions
      
      #Add data from each bird to large output object        
      
      d_Overall <- rbind(d_Overall, d) 
      
    } #End looping over each bird
    
  } #End looping over simulations
  
  #Assign results to output object  
  
  return(d_Overall) 
  
} #End definition of function operations

