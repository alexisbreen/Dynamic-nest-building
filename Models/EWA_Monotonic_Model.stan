///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Script to run the experience-weighted attraction monotonic learning model for the manuscript

//Dynamic strategic social learning in nest-building zebra finches and its generalisability

//Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

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
  
  //Second-time nest construction data
  
  array[N] int id;               //Bird-specific individual identification 
  array[N] int manip;            //Type of experience manipulated: 1 = nest material; 2 = breeding outcome
  array[N] int sat;              //Level of manipulated experience: 1 = satisfied; 2 = dissatisfied
  array[N] int trial;            //Deposit number in-test (1 - 25)
  array[N] int choice;           //Choice of material: 1 = social; 2 = non-social; NOTE: CANNOT 0 & 1 DUMMY CODE CHOICE DATA
  array[N] int cum_soc;          //Cumulative number of social material in nestbox from previous trial to trial 1
  array[N] int cum_non_soc;      //Cumulative number of non-social in nestbox from previous trial to trial 1
  array[N] int social;           //Whether the first social material deposit has happened: 1 = yes; 0 = no
  
} //End block

//Parameter block: Define and name the size of each unobserved variable 
  
parameters{ //Begin block
    
  //Latent learning parameters
  
  matrix[2, 2] logit_phi_first;   //Matrix for latent phi first-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] logit_phi_last;    //Matrix for latent phi last-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] log_lambda_first;  //Matrix for latent lambda first-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] log_lambda_last;   //Matrix for latent lambda last-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] log_epsilon_first; //Matrix for latent epsilon first-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] log_epsilon_last;  //Matrix for latent epsilon last-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] logit_sigma_first; //Matrix for latent sigma first-choice values - indexed by manipulation and satisfaction 
  matrix[2, 2] logit_sigma_last;  //Matrix for latent sigma last-choice values - indexed by manipulation and satisfaction 
  
  // Monotonic effects
  
  simplex[24] delta_phi;      //24-length simplexes...
  simplex[24] delta_lambda;  
  simplex[24] delta_epsilon;       
  simplex[24] delta_sigma;   
  
  //Varying effects clustered on individual; used non-centered approach to estimate individual-level offset as z-scores
  //These z-scores are later multiplied by vector of standard deviations of each parmeter and the cholesky factor to get right covariance structure among parameters
    
  matrix[8, N_id] z_ID;               //Matrix for our latent individual samples (z scores) - indexed by parameter (n = 4) & bird (n = 47)
  vector<lower = 0>[8] sigma_ID;      //Standard deviation of parameters among individuals 
  cholesky_factor_corr[8] Rho_ID; //Cholesky factor for covariance of parameters among individuals

  //Initial preference parameters
  
  matrix[N_id, 2] A_init; //Initial attractions - indexed by bird (n = 47) and material (n = 2)
  real<lower=0> wp;       //Coefficient weight of initial attractions in preference trial
  real ap;                //Intercept for preference trials - just rate of random touching
  
} //End block

//Transformed parameter block: Define and name additional parameters of interest to be saved in posterior output
  
transformed parameters{ //Begin block

  //Saved in output

  matrix[N_id, 8] v_ID;         //Indexed by bird (n = 47) and target parameters (n = 4)
  array[N] vector[2] p;         //Choice probabilities 
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

    for(j in 1:N_id){  //For individual j

      A_base[j, 1] = inv_logit(A_init[j, 1]); //Assign baseline attractions as estimated below
      A_base[j, 2] = inv_logit(A_init[j, 2]); //Assign baseline attractions as estimated below

      //Need to make sure when attractions initialise that A[j,1] < A[j,2]
      //Because columns are differently indexed between estimated baseline versus in-test nest-material attractions: 
      //baseline, col 1 = pink and  col 2 = orange; in-test, col 1 = social & col 2 = non-social
      
      if(A_base[j,1] > A_base[j,2]){ //If pink bigger than orange
        A[j, 1] = A_base[j, 2];      //Pink indexed as non-social material
        A[j, 2] = A_base[j, 1];      //Orange indexed as social material
      } else {                       //If orange bigger than pink
        A[j, 1] = A_base[j, 1];      //Orange indexed as non-social material
        A[j, 2] = A_base[j, 2];      //Pink indexed as social material
          
      } //End ifelse
        
      //Safety check initial attractions assignment - looks good!
        
      //if(A[j,1] > A[j,2]){ 
          //print("A1 > A2");
          //}
            
    } //End j
    
    // Vectors for full 25-trial deltas
    
    vector[25] d_L;
    vector[25] d_P;
    vector[25] d_R;
    vector[25] d_S;
    
    // Fill delta vectors
      
    d_L  = append_row(0, delta_lambda);
    d_P  = append_row(0, delta_phi);
    d_R  = append_row(0, delta_epsilon);
    d_S  = append_row(0, delta_sigma);

    //Loop over Choices
      
    for (i in 1:N) { //For choice i across choices
      
      //Define and name local variables that update across choices
      
      vector[2] pAS; //Vector of asocial choice probabilites
      vector[2] pS;  //Vector of social choice probabilites
      vector[2] pC;  //Vector of choice probabilities
      real L_first;  //Reals for all parameter firsts, lasts and their trial-specific value (captilised)...
      real L_last;
      real LAMBDA;
      real P_first;
      real P_last;
      real PHI;
      real R_first;
      real R_last;
      real epsilon;
      real S_first;
      real S_last;
      real SIGMA;
      
      //Define asocial choice probability i.e., pAS
      
      L_first = exp(log_lambda_first[manip[i], sat[i]] + v_ID[id[i], 1]);
      L_last = exp(log_lambda_last[manip[i], sat[i]] + v_ID[id[i], 2]);

      // Define lambda - capped at 15 because bigger values result in identical choice behaviour
      
      LAMBDA = fmin(15, L_first + (L_last - L_first) * sum(d_L[1:trial[i]]));
      
      //Social choice probability using softmax function that normalizes attraction scores to sum to 1
      //Thus, can interpret it as probability distribution; the ' symbol facilitates matrix multiplication via transposition 
      
      pAS = softmax(LAMBDA * A[id[i], 1:2]'); 
      
      //To determine the log-probability of choice...

      if(social[i] == 0){ //If no social material has yet been deposited in the nestbox
        
      pC = pAS; //Prob. of choice

      } else { //If social material has been previously deposited at least once
      
      R_first = exp(log_epsilon_first[manip[i], sat[i]] + v_ID[id[i], 3] );
      R_last = exp(log_epsilon_last[manip[i], sat[i]] + v_ID[id[i], 4] );

      // Define epsilon - capped at 10 because bigger values result in identical choice behaviour; and bigger values can cause numerical overflow when exponentiated i.e., NaN issues
      
      epsilon = fmin(10, R_first + (R_last - R_first) * sum(d_R[1:trial[i]]));

      // Define social choice probability i.e., pS 
      // Nestbox-counts given very small non-zero for numerical stability 
      // Prevents issues with 0^epsilon, which can cause NaN or undefined behavior
      // Ensures that division by zero does not occur in the probability calculation     
      // Has no effect on true behavior: If cum_soc[i] or cum_non_soc[i] is truly 0, the resulting probability will still behave as if it were 0.
      
      real safe_cum_soc = fmax(cum_soc[i], 1e-6);
      real safe_cum_non_soc = fmax(cum_non_soc[i], 1e-6);

      pS[1] = (safe_cum_soc^epsilon) / (safe_cum_soc^epsilon + safe_cum_non_soc^(1)); // Prob social
      pS[2] = (safe_cum_non_soc^(1)) / (safe_cum_soc^epsilon + safe_cum_non_soc^(1)); // Prob non-social
      
      // Check for NaN issues
      
      if (is_nan(pS[1]) || is_nan(pS[2])) {
          print("Warning: NaN detected in pS computation for i=", i, 
                " epsilon=", epsilon, " cum_soc[i]=", cum_soc[i], " cum_non_soc[i]=", cum_non_soc[i]);
      }
      
      //Define weight of social cues on choice
      
      S_first = inv_logit(logit_sigma_first[manip[i], sat[i]] + v_ID[id[i], 5]);
      S_last = inv_logit(logit_sigma_last[manip[i], sat[i]] + v_ID[id[i], 6]);
      SIGMA = S_first + (S_last - S_first) * sum(d_S[1:trial[i]]);
      
      //Combine above to compute relative choice probabilities i.e., influence of asocial & social cues on choice
        
      pC = (1 - SIGMA) * pAS + SIGMA * pS; //Prob. choice is both social and asocial
      
      } //End determining social vector

      //And update attractions conditional on observed choice
      
      P_first = inv_logit(logit_phi_first[manip[i], sat[i]] + v_ID[id[i], 7]);
      P_last = inv_logit(logit_phi_last[manip[i], sat[i]] + v_ID[id[i], 8]);
      PHI = P_first + (P_last - P_first) * sum(d_P[1:trial[i]]);
      
      if(choice[i] == 1){                                                  //If choice = social
        A[id[i], 1] = ((1-PHI) * A[id[i], 1] + PHI * 1);                   //Atx social gets reward i.e., 1
        A[id[i], 2] = ((1-PHI) * A[id[i], 2] + PHI * 0);                   //Atx non-social no reward i.e., 0  
      } else {                                                             //If choice = non-social
        A[id[i], 1] = ((1-PHI) * A[id[i], 1] + PHI * 0);                   //Atx social gets no reward i.e., 0
        A[id[i], 2] = ((1-PHI) * A[id[i], 2] + PHI * 1);                   //Atx non-social gets reward i.e., 1 
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
  
  to_vector(logit_phi_first) ~  normal(0, 1); 
  to_vector(logit_phi_last) ~  normal(0, 1); 
  delta_phi ~ dirichlet(rep_vector(2.0, 24));

  to_vector(log_lambda_first) ~  normal(0, 1);    
  to_vector(log_lambda_last) ~  normal(0, 1);    
  delta_lambda ~ dirichlet(rep_vector(2.0, 24));

  to_vector(log_epsilon_first)  ~  normal(0, 1); 
  to_vector(log_epsilon_last)  ~  normal(0, 1); 
  delta_epsilon ~ dirichlet(rep_vector(2.0, 24));

  to_vector(logit_sigma_first) ~ normal(0, 1);
  to_vector(logit_sigma_last) ~ normal(0, 1);
  delta_sigma ~ dirichlet(rep_vector(2.0, 24));

  //Define prior distribution of varying individual effects
    
  to_vector(z_ID) ~ normal(0, 1);     //Standard normal prior for z-scores
  sigma_ID ~ exponential(1);          //Exponential prior because variances are bound to be positive
  Rho_ID ~ lkj_corr_cholesky(4);  //Cholesky LKJ correlation distribution for correlation matrix; parameter value = 4 says that more prior probability is placed on small correlations

  //Priors for initial preference parameters
  
  wp ~ exponential(1);              //Standard normal prior for intercept
  ap ~ normal(0, 1);                //Exponential prior b/c probability of event is constant in time - see p. 322 in Rethinking
  to_vector(A_init) ~ normal(0, 1); //Logit scale

  //Initial preference observations
  
  for(i in 1:N_id){
    if(i > 18){ //because birds 1 - 18 have NA for initial material preference data
      touch_P[i] ~ poisson_log(ap + log_dur[i] + wp * inv_logit(A_init[i, 1]));
      touch_O[i] ~ poisson_log(ap + log_dur[i] + wp * inv_logit(A_init[i, 2]));
    }
  }

  //Calculate multionmial likelihood of observed choice
  
  for(i in 1:N){
    choice[i] ~ categorical(p[i]);
  }
  
} 
  
