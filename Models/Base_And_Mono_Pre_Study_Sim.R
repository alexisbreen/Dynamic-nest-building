################################################################################################################################################################################

#Script for the pre-study agent-based forward simulations for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) & Richard McElreath (richard_mcelreath@eva.mpg.de)

################################################################################################################################################################################

#Build simulation function

Pre_study_sim_fct <- function(                 #Begin function definition and parameter imputation
  Mono = 0,                                    #Whether to simulate baseline or monotonic structure (0 = baseline; 1 = monotonic)
  N_choices = 25,                              #Number of nest-material choices
  sd_exp_phi = 0.1,                            #Experiment diffs for phi
  sd_exp_lambda = 0.1,                         #Experiment diffs for lambda
  sd_exp_epsilon = 0.1,                        #Experiment diffs for epsilon
  sd_exp_sigma = 0.1,                          #Experiment diffs for sigma
  sd_id_phi = 0.1,                             #Individual diffs for phi
  sd_id_lambda = 0.1,                          #Individual diffs for lambda
  sd_id_epsilon = 0.1,                         #Individual diffs for epsilon
  sd_id_sigma = 0.1,                           #Individual diffs for sigma
  real_phi = rbind(c(0.5, 0.5), c(0.5, 0.5)),  #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_lambda = rbind(c(3, 3), c(3, 3)),       #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_epsilon = rbind(c(1, 1), c(1, 1)),      #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2) 
  real_sigma = rbind(c(.5, .5), c(.5, .5)))    #(Satisfied exp1, Satisfied exp2),(Dissatisfied exp1, Dissatisfied exp2)  

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
  Log_epsilon <- log(real_epsilon)
  Logit_sigma <- logit(real_sigma)
  
  #Create output object i.e., final data set
  
  d_Overall <- c()
  

  #In each experiment, begin following birds material-choices post-social-demonstration 
    
  for(exp in 1:N_exp){ #Loop over experiments
      
  #Generate between and within-experiment individual-level variation for each latent learning parameter
      
  phi_offset <- rnorm(1, 0, sd_exp_phi) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_phi)
  lambda_offset <- rnorm(1, 0, sd_exp_lambda) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_lambda)
  epsilon_offset <- rnorm(1, 0, sd_exp_epsilon) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_epsilon)
  sigma_offset <- rnorm(1, 0, sd_exp_sigma) + rnorm(sum(N_per_exp[exp, ]), 0, sd_id_sigma)
      
    #Loop over birds
        
    for(ind in 1:sum(N_per_exp[exp, ])){ #For every bird in each experiment
          
      #Assign unique latent parameter value via offsets defined above
          
      phi <- inv_logit(Logit_phi + phi_offset[ind])
      lambda <- exp(Log_lambda + lambda_offset[ind])
      epsilon <- exp(Log_epsilon + epsilon_offset[ind])
      sigma <- inv_logit(Logit_sigma + sigma_offset[ind])
          
      #Generate baseline attraction to nest material
          
      A <- c(0.1, 0.7) #Low social material attraction and high non-social material attraction
          
      #Create output matrix to record choices 
          
      d <- data.frame(id = sum(N_per_exp[which(1:2 < exp), ]) + ind, 
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
          
      #Loop over material-choice trials
          
      for(trial in 1:N_choices){
            
        #First, before each choice...
        #Determine whether two material-types in nestbox & add to data frame
        #If two material-types in nestbox, activate & calculate NB matrix to use in pA
            
        if(trial > 1){
          NB[1] <- length(which(d$copy[- trial] == 1))
          NB[2] <- length(which(d$copy[- trial] == 0))
        }
            
        #Decide if nestbox contains at least one social material, and if so, turn social vector 'on'
            
        if(sum(NB[1]) >= 1) soc <- 1
        d$social[trial] <- soc
            
        #Decide epsilon based on model structure
            
        e_t <- c() 
        if(Mono == 0){
          e_t <- epsilon[d$sat_level[trial], d$experiment[trial]]
        } else {
          if (d$sat_level[trial] == 1) {                             #Satisfied: epsilon decreases over time
            e_t <- epsilon[d$sat_level[trial], d$experiment[trial]] * (1 - 0.8 * ((trial - 1) / (N_choices - 1))) #Decrease by ca. 80%
          } else {                                                   #Dissatisfied: epsilon increases over time
            e_t <- epsilon[d$sat_level[trial], d$experiment[trial]] * (1 + 4 * ((trial - 1) / (N_choices - 1)))   #Increase by ca. 500%
          }
        }
            
        #Calculate decision probabilities 
            
        Prob <- c() #Empty vector to fill with decision probabilities
        if(soc == 0){
          Prob[1] <- exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[1]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) 
          Prob[2] <- exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[2]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) 
        } else {
          Prob[1] <- (1 - sigma[d$sat_level[trial], d$experiment[trial]]) * exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[1]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) + (sigma[d$sat_level[trial], d$experiment[trial]] * ((NB[1] ^ e_t) / ( (NB[1] ^ e_t) + (NB[2] ^ (1)) ))) 
          Prob[2] <- (1 - sigma[d$sat_level[trial], d$experiment[trial]]) * exp(lambda[d$sat_level[trial], d$experiment[trial]] * A[2]) / sum(exp(lambda[d$sat_level[trial], d$experiment[trial]] * A)) + (sigma[d$sat_level[trial], d$experiment[trial]] * ((NB[2] ^ (1)) / ( (NB[1] ^ e_t) + (NB[2] ^ (1)) ))) 
              
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
            
        #Assign attraction scores to data frame 
            
        d$Atx_soc[trial] <- A[1]
        d$Atx_non_soc[trial] <- A[2]
            
      } #End looping over material-choice decisions
          
      #Add data from each bird to large output object
          
      d_Overall <- rbind(d_Overall, d) 
          
    } #End looping over each bird
      
  } #End looping over each experiment
  
  #Assign results to output object  
  
  return(d_Overall) 
  
} #End definition of function operations
