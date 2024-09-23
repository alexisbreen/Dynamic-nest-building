################################################################################################################################################################################

#Script to run post-study forward simulations for the manuscript

#Strategic social learning steers animal material technology

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de)

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#The following script must have already been run in order to fully execute the current script:

#SSL_EWA_Model_Execution.R

#
##
###Build simulation function
##
#

#To simulate new choosers from the estimated population, 
#first derive the full posterior distribution with random effects variance-covariance matrix 

#Extract posterior from EWA model

s <- extract.samples(m)

#1) Get standard deviation of learning parameters (variance) among individuals

sigma_id <- s$sigma_ID

#2) Compute correlation matrix from Cholesky factors

Correlations <- array(numeric(),c(4000,4,4))
for (i in 1:4000) {
  Correlations[i,,] <- s$Rho_ID[i,,] %*% t(s$Rho_ID[i,,])
}

#3) Compute variance-covariance matrix from variances and correlations

S <- array(numeric(),c(4000,4,4))
for (i in 1:4000) {
  S[i,,] <- diag(sigma_id[i,]) %*% Correlations[i,,]  %*% diag(sigma_id[i,])
}

#Extract posteriors for means

#Phi

phi_post <- array(numeric(), c(4000,2,2))
phi_post[,1,1] <- s$logit_phi[,1,1] #Satisfied-construction
phi_post[,1,2] <- s$logit_phi[,1,2] #Dissatisfied-construction
phi_post[,2,1] <- s$logit_phi[,2,1] #Satisfied-reproduction
phi_post[,2,2] <- s$logit_phi[,2,2] #Dissatisfied-reproduction

#Lambda

lambda_post <- array(numeric(), c(4000,2,2))
lambda_post[,1,1] <- s$log_lambda[,1,1] #Satisfied-construction
lambda_post[,1,2] <- s$log_lambda[,1,2] #Dissatisfied-construction
lambda_post[,2,1] <- s$log_lambda[,2,1] #Satisfied-reproduction
lambda_post[,2,2] <- s$log_lambda[,2,2] #Dissatisfied-reproduction
 
#Rho

rho_post <- array(numeric(), c(4000,2,2))
rho_post[,1,1] <- s$log_rho[,1,1] #Satisfied-construction
rho_post[,1,2] <- s$log_rho[,1,2] #Dissatisfied-construction
rho_post[,2,1] <- s$log_rho[,2,1] #Satisfied-reproduction
rho_post[,2,2] <- s$log_rho[,2,2] #Dissatisfied-reproduction

#Sigma

sigma_post <- array(numeric(), c(4000,2,2))
sigma_post[,1,1] <- s$logit_sigma[,1,1] #Satisfied-construction
sigma_post[,1,2] <- s$logit_sigma[,1,2] #Dissatisfied-construction
sigma_post[,2,1] <- s$logit_sigma[,2,1] #Satisfied-reproduction
sigma_post[,2,2] <- s$logit_sigma[,2,2] #Dissatisfied-reproduction

#Build post-hoc forward simulation function

PH_sim_fct <- function( #Begin function definition and parameter imputation
  N_sim = 1,            #Number of simulations to run
  N_choices = 25,       #Number of nest-material choices
  N_choosers = 500,     #Total number of choosers to draw from posterior - must be divisible by four! 
  N_soc_pay = 1,        #Reward payoff for choosing social material
  N_non_soc_pay = 1,    #Reward payoff for choosing non-social material
  SS_matched = 0)       #Whether to simulate within-treatment numbers matched to true sample size or use N_choosers value (0 = use N_choosers; 1 = use matched)

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
    
    Sim_S1 <- t(sapply(1:N_S1, function(i) rmvnorm(1, c(phi_post[draw_S1[i],1,1], lambda_post[draw_S1[i],1,1], rho_post[draw_S1[i], 1, 1], sigma_post[draw_S1[i], 1, 1]), S[draw_S1[i],,]))) #Satisfied-construction
    Sim_D1 <- t(sapply(1:N_D1, function(i) rmvnorm(1, c(phi_post[draw_D1[i],1,2], lambda_post[draw_D1[i],1,2], rho_post[draw_D1[i], 1, 2], sigma_post[draw_D1[i], 1, 2]), S[draw_D1[i],,]))) #Dissatisfied-construction
    Sim_S2 <- t(sapply(1:N_S2, function(i) rmvnorm(1, c(phi_post[draw_S2[i],2,1], lambda_post[draw_S2[i],2,1], rho_post[draw_S2[i], 2, 1], sigma_post[draw_S2[i], 2, 1]), S[draw_S2[i],,]))) #Satisfied-reproduction
    Sim_D2 <- t(sapply(1:N_D2, function(i) rmvnorm(1, c(phi_post[draw_D2[i],2,2], lambda_post[draw_D2[i],2,2], rho_post[draw_D2[i], 2, 2], sigma_post[draw_D2[i], 2, 2]), S[draw_D2[i],,]))) #Dissatisfied-reproduction
    
    #Combine draws
    
    Sim_all <- rbind(Sim_S1, Sim_D1, Sim_S2, Sim_D2)
      
    #Index draws to target parameters & transform back to original scale to be used later in the model
    
    phi <- inv_logit(Sim_all[, 1])
    lambda <- exp(Sim_all[, 2])
    rho <- exp(Sim_all[, 3])
    sigma <- inv_logit(Sim_all[, 4])
    
    #lambda and rho can 'run away' because past a certain number the model 'thinks' any draw is 'big'
    #This causes overflow when executing the equations below, which needs to be dampened by reducing draws past a reasonably high number - here, 200
    
    lambda <- ifelse(lambda > 200, lambda/100000000, lambda) #Huge reducer number b/c in simulation 2 the reward value i.e., multiplier of 2 in equation below quickly causes overflow 
    rho <- ifelse(rho > 200, rho/100000000, rho) #Huge reducer number b/c in simulation 2 the reward value i.e., multiplier of 2 in equation below quickly causes overflow 
      
    #Empty attraction matrix to hold initial attraction draws of social and non-social material
    
    Atx <- matrix(0, N_id, 2)
    
    #When making attraction draws, need to specify from which birds the estimates are to be pulled
    #This is because the attraction estimates are indexed in stan by choice and bird id
    #And if the forward simulation is not matched to study sample size, the attraction draws return 'subscript out of bounds' - for example, when trying to find bird 48, when we only have a maximum of 47 bird attractions estimated from true sample
    #So, we decide how attractions are to be drawn based on whether we're matched (or not) to study sample size
    
    #Combine our previously defined treatment draws into one vector, so that our for loops (below), which loop through birds within treatments, draw attractions from the correct bird and treatment, because attractions in stan are indexed by choice and bird id
    
    draw_all <- c(draw_S1,draw_D1,draw_S2,draw_D2) 
    
    if(SS_matched == 1){ #If matched...
      
      #Draw attractions
      #We draw from estimated attraction scores at deposit 1 (versus estimated baseline attraction before deposit 1)
      #We do this because these draws will be closer to true attractions due to recursive nature of estimation - see main text
      
      for(i in 1:N_S1){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, i], s$Atx_non_soc[draw_all[i], 1, i])} #Satisfied-construction
      for(i in (max(N_S1) + 1):(max(N_S1) + max(N_D1))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, i], s$Atx_non_soc[draw_all[i], 1, i])} #Dissatisfied-construction
      for(i in (max(N_S1) + max(N_D1) + 1):(max(N_S1) + max(N_D1) + max(N_S2))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, i], s$Atx_non_soc[draw_all[i], 1, i])} #Satisfied-reproduction
      for(i in (max(N_S1) + max(N_D1) + max(N_S2) + 1):(max(N_S1) + max(N_D1) + max(N_S2) + max(N_D2))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, i], s$Atx_non_soc[draw_all[i], 1, i])} #Dissatisfied-reproduction
      
    } else { #If not matched
      
      draw_S1_A <- sample(1:10, N_S1, replace = TRUE)
      draw_D1_A <- sample(11:18, N_D1, replace = TRUE)
      draw_S2_A <- sample(19:32, N_S2, replace = TRUE)
      draw_D2_A <- sample(33:47, N_D2, replace = TRUE)
      
      id_A <- c(draw_S1_A,draw_D1_A,draw_S2_A,draw_D2_A) 
      
      #Draw attractions
      #We draw from estimated attraction scores at deposit 1 (versus estimated baseline attraction before deposit 1)
      #We do this because these draws will be closer to true attractions due to recursive nature of estimation - see main 
      
      for(i in 1:N_S1){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, id_A[i]], s$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Satisfied-construction
      for(i in (max(N_S1) + 1):(max(N_S1) + max(N_D1))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, id_A[i]], s$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Dissatisfied-construction
      for(i in (max(N_S1) + max(N_D1) + 1):(max(N_S1) + max(N_D1) + max(N_S2))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, id_A[i]], s$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Satisfied-reproduction
      for(i in (max(N_S1) + max(N_D1) + max(N_S2) + 1):(max(N_S1) + max(N_D1) + max(N_S2) + max(N_D2))){Atx[i, ] <- c(s$Atx_soc[draw_all[i], 1, id_A[i]], s$Atx_non_soc[draw_all[i], 1, id_A[i]])} #Dissatisfied-reproduction
      
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
        
        #Calculate decision probabilities 
       
        Prob <- c() #Empty vector to fill with decision probabilities
        if(soc == 0){
          Prob[1] <- exp(lambda[ind] * A[1]) / sum(exp(lambda[ind] * A)) 
          Prob[2] <- exp(lambda[ind] * A[2]) / sum(exp(lambda[ind] * A)) 
        } else {
          Prob[1] <- (1 - sigma[ind]) * exp(lambda[ind] * A[1]) / sum(exp(lambda[ind] * A)) + (sigma[ind] * ((NB[1] ^ rho[ind]) / ( (NB[1] ^ rho[ind]) + (NB[2] ^ (1)) ))) 
          Prob[2] <- (1 - sigma[ind]) * exp(lambda[ind] * A[2]) / sum(exp(lambda[ind] * A)) + (sigma[ind] * ((NB[2] ^ (1)) / ( (NB[1] ^ rho[ind]) + (NB[2] ^ (1)) ))) 
          
        }
        
        #Make choice proportional to attraction scores
        
        d$choice[which(d$trial == trial)] <- sample(c(1:2), size = 1, prob = Prob) #Make the bird choose among the two choice-options, based on the choice probabilities defined above    
          
        #Determine copy based on choice
        
        if(d$choice[trial] == 1){ #If the bird chose the social option
          d$copy[trial] <- 1      #The bird is assigned to have copied 
        } else {                  #But if the bird chose the non-social option
          d$copy[trial] <- 0      #The bird is assigned to have not copied 
        }                         #End determination 
          
        #Update attractions
        
        if(d$choice[trial] == 1){                                                #If choose social
          A[1] <- (1 - phi[ind]) * A[d$choice[trial]] + phi[ind] * N_soc_pay     #Attraction for social increase
          A[2] <- (1 - phi[ind]) * A[d$choice[trial]] + phi[ind] * 0             #Attraction for non-social decrease
        } else {                                                                 #If choose non-social
          A[1] <- (1 - phi[ind]) * A[d$choice[trial]] + phi[ind] * 0             #Attraction for social decrease
          A[2] <- (1 - phi[ind]) * A[d$choice[trial]] + phi[ind] * N_non_soc_pay #Attraction for non-social increase
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
