################################################################################################################################################################################

#Script for the pre-study agent-based forward simulations for the manuscript

#Strategic social learning in avian nest construction and potentially beyond

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de)

################################################################################################################################################################################

#
##
###Build simulation function
##
#

PS_sim_fct <- function(                        #Begin function definition and parameter imputation
  N_sim = 2,                                   #Number of simulations: low and high social material use
  sim_N = 1,                                   #Which simulation; 1 = low social material use; 2 = high social material use
  N_choices = 25,                              #Number of nest-material choices
  sd_exp_phi = 0.1,                            #Experiment diffs for phi
  sd_exp_lambda = 0.1,                         #Experiment diffs for lambda
  sd_exp_rho = 0.1,                            #Experiment diffs for rho
  sd_exp_sigma = 0.1,                          #Experiment diffs for sigma
  sd_id_phi = 0.1,                             #Individual diffs for phi
  sd_id_lambda = 0.1,                          #Individual diffs for lambda
  sd_id_rho = 0.1,                             #Individual diffs for rho
  sd_id_sigma = 0.1,                           #Individual diffs for sigma
  real_phi = rbind(c(0.2, 0.2), c(0.2, 0.2)),  #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_lambda = rbind(c(5, 5), c(5, 5)),       #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_rho = rbind(c(3, 3), c(3, 3)),          #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_sigma = rbind(c(.6, .6), c(.6, .6)))    #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2)  

{ #Begin definition of function operations
  
  #Simulate subject birds corresponding to data set 
  
  N_exp <- 2                    #Experiment number
  N_per_exp <- matrix(0, 2, 2)  #Empty matrix to hold no. satisfied-dissatisfied birds per experiment
  N_per_exp[1, ] <- c(10, 8)    #(Satisfied, Dissatisfied)
  N_per_exp[2, ] <- c(14, 15)   #(Satisfied, Dissatisfied)
  
  #Define unique bird ID parameters based on above matrix 
  
  N_id <- sum(N_per_exp)    
  
  #Define scale of our phi and lambda values; transform to latent logit/log scale to add random offsets - see below
  
  Logit_phi <- logit(real_phi)
  Log_lambda <- log(real_lambda)
  Log_rho <- log(real_rho)
  Logit_sigma <- logit(real_sigma)
  
  #Create output object i.e., final dataset
  
  d_Overall <- c()
  
  for(sim in 1:N_sim){
    
    #In each experiment, begin following birds material-choices post-social-demonstration 
    #Loop over experiments
    
    for(exp in 1:N_exp){
      
      #Generate between and within-experiment individual-level variation for each latent parameter
      
      phi_offset <- rnorm(1, 0, sd_exp_phi) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_phi)
      lambda_offset <- rnorm(1, 0, sd_exp_lambda) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_lambda)
      rho_offset <- rnorm(1, 0, sd_exp_rho) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_rho)
      sigma_offset <- rnorm(1, 0, sd_exp_sigma) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_sigma)
      
      #Loop over birds
      
      for(ind in 1:sum(N_per_exp[exp, ])){ #For every bird in each experiment
        
        #Assign unique latent parameter value via offsets defined above
        
        phi <- inv_logit(Logit_phi + phi_offset[ind])
        lambda <- exp(Log_lambda + lambda_offset[ind])
        rho <- exp(Log_rho + rho_offset[ind])
        sigma <- inv_logit(Logit_sigma + sigma_offset[ind])
        
        #Generate baseline attraction to nest material, based on which pre-study simulation (1 or 2) is running
        
        if(sim_N == 1){
          A <- c(0.1, 0.7) #Low social material attraction and high non-social material attraction
        } else {
          A <- c(0.7, 0.1) #High social material attraction and low non-social material attraction
        } 
        
        #Create output matrix to record choices 
        
        d <- data.frame(sim = sim,
                        id = sum(N_per_exp[which(1:2 < exp), ]) + ind, 
                        sat_level = ifelse(ind <= N_per_exp[exp, 1], 1, 2),
                        experiment = exp, 
                        trial = 1:N_choices, 
                        cum_soc = NA,
                        cum_non_soc = NA,
                        choice = NA, 
                        copy = NA,
                        social = NA,
                        prob_soc = NA,
                        prob_non_soc = NA,
                        Atx_soc = NA,
                        Atx_non_soc =NA)
        
        #Define local variables used in material-choice decision loop      
        
        NB <- matrix(0, 1, 2) #Empty matrix to fill with nestbox counts of each material-type; used when social = 1 for Prob determination
        soc <- 0
        
        #Loop over material-choice decisions i.e., trials
        
        for(trial in 1:N_choices){
          
          #First, before each choice...
          #Determine whether two material-types in nestbox & add to data frame
          #If two material-types in nestbox, activate & calculate NB matrix to use in pA
          
          if(trial > 1){
            NB[1] <- length(which(d$copy[- trial] == 1))
            NB[2] <- length(which(d$copy[- trial] == 0))
          }
          
          #Decide if nestbox contains at least one social material, and, if so, turn social vector 'on'
          
          if(sum(NB[1]) >= 1) soc <- 1
          d$social[trial] <- soc
          
          #Calculate decision probabilities 
          
          Prob <- c() 
          if(soc == 0){
            Prob[1] <- exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[1]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) 
            Prob[2] <- exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[2]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) 
          } else {
            Prob[1] <- (1 - sigma[d$sat_level[trial], d$experiment[trial]]) * exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[1]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) + (sigma[d$sat_level[trial], d$experiment[trial]] * ((NB[1] ^ rho[d$sat_level[trial], d$experiment[trial]]) / ( (NB[1] ^ rho[d$sat_level[trial], d$experiment[trial]]) + (NB[2] ^ (1)) ))) 
            Prob[2] <- (1 - sigma[d$sat_level[trial], d$experiment[trial]]) * exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[2]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) + (sigma[d$sat_level[trial], d$experiment[trial]] * ((NB[2] ^ (1)) / ( (NB[1] ^ rho[d$sat_level[trial], d$experiment[trial]]) + (NB[2] ^ (1)) ))) 
            
          }
          
          #Store each probability separately for posterior prediction checks
          
          d$prob_soc[trial] <- Prob[1]
          d$prob_non_soc[trial] <- Prob[2]
          
          #Make choice proportional to attraction scores
          
          d$choice[which(d$trial == trial)] <- sample(c(1:2), size = 1, prob = Prob) #Make the bird choose among the two choice-options, based on the choice probabilities defined above    
          
          #Determine copy based on choice
          
          if(d$choice[trial] == 1){ #If the bird chose the social option
            d$copy[trial] <- 1      #The bird is assigned to have copied 
          } else {                  #But if the bird chose the non-social option
            d$copy[trial] <- 0      #The bird is assigned to have not copied 
          }                         #End determination 
          
          #Determine cumulative sum variables based on copying
          #We need these variables for the stan model to run
          
          if(trial > 1){                                                #If trial is greater than 1
            d$cum_soc[trial] <- length(which(d$copy[- trial] == 1))     #Count number of previous social material copies & assign to cum_soc
            d$cum_non_soc[trial] <- length(which(d$copy[- trial] == 0)) #Count number of previous non-social material copies & assign to cum_non_soc
          } else {                                                      #If the trial equals 1
            d$cum_soc[trial] <- 0                                       #Not possible to have social material previously deposited - assign 0
            d$cum_non_soc[trial] <- 0                                   #Not possible to have non-social material previously deposited - assign 0
          }                                                             #End determination 
          
          #Update...
          
          #Attraction for chosen option 
          
          A[d$choice[trial]] <- (1 - phi[d$sat_level[trial], d$experiment[trial]]) * A[d$choice[trial]] + phi[d$sat_level[trial], d$experiment[trial]] * 1
          
          #Attraction for non-chosen option
          
          if(d$choice[trial] == 1){
            A[2] <- (1 - phi[d$sat_level[trial], d$experiment[trial]]) * A[d$choice[trial]] + phi[d$sat_level[trial], d$experiment[trial]] * 0
          } else {
            A[1] <- (1 - phi[d$sat_level[trial], d$experiment[trial]]) * A[d$choice[trial]] + phi[d$sat_level[trial], d$experiment[trial]] * 0
          }
          
          ##If Simulation 1...because this is a reinforcement learning simulation where all choices reward (and all non-choices return 0),
          #the simulation will not approximate low levels of social material use, if the attraction scores aren't dampened,
          #no matter how much we set the initial latent learning parameters to be 'anti' social material
          
          if(sim_N == 1) A[1] <- A[1] / 10 #Dampen...
          
          #Assign attraction scores to data frame 
          
          d$Atx_soc[trial] <- A[1]
          d$Atx_non_soc[trial] <- A[2]
          
        } #End looping over material-choice decisions
        
        #Add data from each bird to large output object  
        
        d_Overall <- rbind(d_Overall, d) 
        
      } #End looping over each bird
      
    } #End looping over each experiment
  }
  
  #Assign results to output object  
  return(d_Overall) 
  
} #End definition of function operations
