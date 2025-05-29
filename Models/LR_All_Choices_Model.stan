/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Script to run the hierarchical logistic regression all choices model for the manuscript

//Dynamic strategic social learning in nest-building zebra finches and its generalisability

//Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

data{
  int<lower=1> N;                             //Total number of observations (trials)
  int<lower=1> N_id;                          //Total number of unique individual birds
  int<lower=1> N_treat;                       //Number of treatment groups (4)
  array[N] int<lower=1, upper=N_id> id;       //Bird ID for each observation (1-47)
  array[N] int<lower=1, upper=N_treat> treat; //Treatment group index for each observation (1-4)
  array[N] int<lower=0, upper=1> choice;      //Binary choice outcome: 1 = chose social material, 0 = chose asocial material
  array[N] int<lower=1, upper=25> trial;      //Trial number (1-25)
}

parameters{
  vector[N_treat] a;                //Treatment-specific intercepts (baseline preference for social material per treatment)
  real beta_trial;                  //Global slope for trial number (average learning effect across all treatments)
  vector[N_treat] beta_treat_trial; //Treatment-specific deviations in learning slope (treatment × trial interaction)

  vector[N_id] b;                   //Random intercepts for individuals (captures individual baseline preference)
  vector[N_id] b_trial;             //Random slopes for individuals (captures variation in individual learning rates over trials)

  real<lower=0> sigma_b;            //Standard deviation of individual intercepts
  real<lower=0> sigma_b_trial;      //Standard deviation of individual slopes for trial
}

model{
  vector[N] logit_pi;               //Linear predictor for the log-odds of making a social choice

  //Priors on fixed effects 
  
  a ~ normal(0, 1);                //Weakly regularizing prior on treatment intercepts (centers preference at log-odds = 0)
  beta_trial ~ normal(0, 1);       //Prior on overall effect of trial number 
  beta_treat_trial ~ normal(0, 1); //Prior on interaction between treatment and trial effect

  //Priors on random effects
  
  b ~ normal(0, sigma_b);             //Hierarchical prior for individual baseline preference
  b_trial ~ normal(0, sigma_b_trial); //Hierarchical prior for individual learning slopes
  sigma_b ~ exponential(1);           //Prior for variability in individual intercepts
  sigma_b_trial ~ exponential(1);     //Prior for variability in individual trial slopes

  //Linear predictor: combines fixed and random effects
  
  for(n in 1:N){
    logit_pi[n] = a[treat[n]] +                           //Treatment-specific intercept
                  b[id[n]] +                              //Individual deviation from group mean
                  beta_trial * trial[n] +                 //Global learning effect
                  beta_treat_trial[treat[n]] * trial[n] + //Treatment-specific learning adjustment
                  b_trial[id[n]] * trial[n];              //Individual learning slope
  }

  //Likelihood: logistic regression for binary choice
  
  choice ~ bernoulli_logit(logit_pi);
}

generated quantities{
  vector[N] prob_social; //Posterior predictive probabilities of choosing social info

  //Compute posterior predicted probabilities for each observation
  
  for(n in 1:N){
    prob_social[n] = inv_logit(a[treat[n]] +                           //Treatment effect
                               b[id[n]] +                              //Individual baseline
                               beta_trial * trial[n] +                 //Overall learning effect
                               beta_treat_trial[treat[n]] * trial[n] + //Treatment-modulated learning
                               b_trial[id[n]] * trial[n]);             //Individual-specific learning
  }
}
