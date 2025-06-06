///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Script to run the experience-weighted attraction learning model for the manuscript

//Dynamic strategic social learning in nest-building zebra finches and its generalisability

//Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) and Richard McElreath (richard_mcelreath@eva.mpg.de)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Begin stan decision model

//Data block: Define and name the size of each observed variable
  
data{ //Begin block
  
  //Indexing data
  
  int N;                         //Number of observations
  int N_id;                      //Number of birds
  
  //Initial preference data
  
  array[N_id] int touch_P;       //Total touches to pink 
  array[N_id] int touch_O;       //Total touches to orange
  array[N_id] real log_dur;      //Duration of test in seconds (log transformed)
  
  //Second-time nest-building data
  
  array[N] int id;              //Bird-specific individual identification 
  array[N] int manip;           //Type of experience manipulated: 1 = nest material; 2 = breeding outcome
  array[N] int sat;             //Level of manipulated experience: 1 = satisfied; 2 = dissatisfied
  array[N] int trial;           //Deposit number in-test (1 - 25)
  array[N] int choice;          //Choice of material: 1 = social; 2 = non-social; NOTE: CANNOT 0 & 1 DUMMY CODE CHOICE DATA
  array[N] int cum_soc;         //Cumulative number of social material in nestbox from previous trial to trial 1
  array[N] int cum_non_soc;     //Cumulative number of non-social in nestbox from previous trial to trial 1
  array[N] int social;          //Whether the first social material deposit has happened: 1 = yes; 0 = no
  
} //End block

//Parameter block: Define and name the size of each unobserved variable 
  
parameters{ //Begin block
    
  //Latent learning parameters
  
  matrix[2, 2] logit_phi;       //Matrix for latent phi values (asocial updating) - indexed by manipulation and satisfaction 
  matrix[2, 2] log_lambda;      //Matrix for latent lambda values (asocial weighting) - indexed by manipulation and satisfaction 
  matrix[2, 2] log_epsilon;     //Matrix for latent epsilon values (social updating) - indexed by manipulation and satisfaction 
  matrix[2, 2] logit_sigma;     //Matrix for latent sigma values (social weighting) - indexed by manipulation and satisfaction 

  //Varying effects clustered on individual; used non-centered approach to estimate individual-level offset as z-scores
  //These z-scores are later multiplied by vector of standard deviations of each parmeter and the cholesky factor to get right covariance structure among parameters
    
  matrix[4, N_id] z_ID;               //Matrix for our latent individual samples (z scores) - indexed by parameter (n = 4) & bird (n = 47)
  vector<lower = 0>[4] sigma_ID;      //Standard deviation of parameters among individuals 
  cholesky_factor_corr[4] Rho_ID; //Cholesky factor for covariance of parameters among individuals

  //Initial preference parameters
  
  matrix[N_id, 2] A_init; //Initial attractions - indexed by bird (n = 47) and material (n = 2)
  real<lower=0> wp;       //Coefficient weight of initial attractions in preference trial
  real ap;                //Intercept for preference trials - just rate of random touching
  
} //End block

//Transformed parameter block: Define and name additional parameters of interest to be saved in posterior output
  
transformed parameters{ //Begin block

  //Saved in output

  matrix[N_id, 4] v_ID;         //Indexed by bird (n = 47) and target parameters (n = 4)
  array[N] simplex[2] p;        //Choice probabilities 
  vector[N] soc_p;              //Social choice probability seperate
  vector[N] non_soc_p;          //Asocial choice probability seperate
  matrix[25, N_id] Atx_soc;     //Attraction score for social material - indexed by choices (n = 25) and bird (n = 47)
  matrix[25, N_id] Atx_non_soc; //Attraction score for non-social material - indexed by choices (n = 25) and bird (n = 47)
  matrix[25, N_id] soc_pp;      //Posterior prediction for social material - indexed by choices (n = 25) and bird (n = 47)
  
  v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)'; //Variance-covariance matrix clustered on individuals; based on z-scores, standard deviations and Cholesky factors - see: p. 467 in Rethinking (https://github.com/Booleans/statistical-rethinking/blob/master/Statistical%20Rethinking%202nd%20Edition.pdf)

  { //Begin local scope i.e, nothing below saved in output 

    //Initialise attraction scores via estimation & assign

    matrix[N_id, 2] A;      //Empty matrix to store in-test attraction scores for each material
    matrix[N_id, 2] A_base; //Empty matrix to hold baseline attraction scores for each material
    
    //Loop over individuals

    for(i in 1:N_id){  //For observation i across individuals
    
      A_base[i, 1] = inv_logit(A_init[i, 1]); //Assign baseline attractions as estimated below
      A_base[i, 2] = inv_logit(A_init[i, 2]); //Assign baseline attractions as estimated below

      //Need to make sure when attractions initialise that A[i,1] < A[i,2]
      //Because j is differently indexed between estimated baseline versus in-test nest-material attractions: 
      //baseline, j1 = pink and j2 = orange; in-test, j1 = social & j2 = non-social
      
      if(A_base[i,1] > A_base[i,2]){ //If pink bigger than orange
        A[i, 1] = A_base[i, 2];      //Pink indexed as non-social material
        A[i, 2] = A_base[i, 1];      //Orange indexed as social material
      } else {                       //If orange bigger than pink
        A[i, 1] = A_base[i, 1];      //Orange indexed as non-social material
        A[i, 2] = A_base[i, 2];      //Pink indexed as social material
          
      } //End ifelse
        
      //Safety check initial attractions assignment - looks good!
        
      //if(A[i,1] > A[i,2]){ 
          //print("A1 > A2");
          //}
            
    } //End loop over individuls

    //Loop over Choices
      
    for (i in 1:N){ //For choice i across choices
      
      //Define and name local variables that update across choices
      
      vector[2] pA;  //Vector of asocial choice probabilites
      vector[2] pS;  //Vector of social choice probabilites
      vector[2] pC;  //Vector of choice probabilities
      real phi;      //Asocial updating - on outcome scale as used in the model
      real lambda;   //Asocial weighting - on outcome scale as used in the model
      real epsilon;  //Social updating - on outcome scale as used in the model
      real sigma;    //Social weighting - on outcome scale as used in the model
      
      //Define asocial choice probability i.e., pA
      
      lambda = fmin(exp(log_lambda[manip[i], sat[i]] + v_ID[id[i], 1]), 15); //Cap at 15
      pA = softmax(lambda * A[id[i], 1:2]'); //Social choice probability using softmax function that normalizes attraction scores to sum to 1, so can interpret it as probability distribution; the ' symbol facilitates matrix multiplication via transposition                         
      
      //To determine the log-probability of choice...

      if(social[i] == 0){ //If no social material has yet been deposited in the nestbox
        
      pC = pA; //Prob. of choice
        
      } else { //If social material has been previously deposited at least once
      
      //Define social choice probability i.e., pS
      
      epsilon = fmin(exp(log_epsilon[manip[i], sat[i]] + v_ID[id[i], 2]/1000), 10); //Cap at 10
      pS[1] = (cum_soc[i]^epsilon) / (cum_soc[i]^epsilon + cum_non_soc[i]^(1));     //Prob. social material
      pS[2] = (cum_non_soc[i]^(1)) / (cum_soc[i]^epsilon + cum_non_soc[i]^(1));     //Prob. non-social material
        
      //Define weight of social cues on choice
          
      sigma = inv_logit(logit_sigma[manip[i], sat[i]] + v_ID[id[i], 3]); //Main and varying effects on sigma 
          
      //Combine above to compute relative choice probabilities i.e., influence of asocial & social cues on choice
    
      pC = (1 - sigma) * pA + sigma * pS; //Prob. choice is both social and asocial
        
      } //End determining social vector

      //And update attractions conditional on observed choice
        
      phi = inv_logit(logit_phi[manip[i], sat[i]] + v_ID[id[i], 4]/1000);  //Divided by 1000 to dampen overflow 
      if(choice[i] == 1){                                                  //If choice = social
        A[id[i], 1] = ((1-phi) * A[id[i], 1] + phi * 1);                   //Atx social gets reward i.e., 1
        A[id[i], 2] = ((1-phi) * A[id[i], 2] + phi * 0);                   //Atx non-social no reward i.e., 0  
      } else {                                                             //If choice = non-social
        A[id[i], 1] = ((1-phi) * A[id[i], 1] + phi * 0);                   //Atx social gets no reward i.e., 0
        A[id[i], 2] = ((1-phi) * A[id[i], 2] + phi * 1);                   //Atx non-social gets reward i.e., 1 
      }
      
      //Finally, assign target local parameters to global environment to save in output
      
      p[i] = pC;
      soc_p[i] = pC[1];
      non_soc_p[i] = pC[2];
      Atx_soc[trial[i], id[i]] = A[id[i], 1];
      Atx_non_soc[trial[i], id[i]] = A[id[i], 2];
      soc_pp[trial[i], id[i]] = pC[1];

    } //End looping over choices
  } //End local scope
} //End block
  
//Model block: Name and define the model
  
model{ //Begin block

  //Define (weakly-regularizing) priors for target parameters
  //Turn into vectors so that each can be used as a vectorized argument to the univariate normal density
  
  to_vector(logit_phi) ~  normal(0, 1);     //Standard normal priors...
  to_vector(log_lambda) ~  normal(0, 1);    
  to_vector(log_epsilon)  ~  normal(0, 1);   
  to_vector(logit_sigma) ~ normal(0, 1);
  
  //Define prior distribution of varying individual effects
    
  to_vector(z_ID) ~ normal(0, 1);     //Standard normal prior for z-scores
  sigma_ID ~ exponential(1);          //Exponential prior because variances are bound to be positive
  Rho_ID ~ lkj_corr_cholesky(4);      //Cholesky LKJ correlation distribution for correlation matrix; parameter value = 4 says that more prior probability is placed on small correlations

  //Priors for initial preference parameters
  
  wp ~ exponential(1);              //Standard normal prior for intercept
  ap ~ normal(0, 1);                //Exponential prior b/c probability of event is constant in time 
  to_vector(A_init) ~ normal(0, 1); //Logit scale

  //Initial preference observations
  
  for(i in 1:N_id){
    if(i > 18){ //Because birds 1 - 18 have NA for initial material preference data
      touch_P[i] ~ poisson_log(ap + log_dur[i] + wp * inv_logit(A_init[i, 1]));
      touch_O[i] ~ poisson_log(ap + log_dur[i] + wp * inv_logit(A_init[i, 2]));
    }
  }

  //Calculate multionmial likelihood of observed choice
  
  for(i in 1:N){
    choice[i] ~ categorical(p[i]);
  }
  
}