///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Script to run the hierarchical logistic regression first choice model for the manuscript

//Dynamic strategic social learning in nest-building zebra finches and its generalisability

//Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

data {
  int<lower=1> N;                                     //Number of observations (trials)
  array[N] int<lower=1, upper=4> treat;               //Treatment group index for each observation (1 through 4)
  array[N] int<lower=0, upper=1> choice;              //Binary outcome for each observation (1 = chose social, 0 = chose asocial)
}

parameters {
  vector[4] alpha;                                    // Log-odds intercepts for each treatment group
}

model {
  alpha ~ normal(0, 1);                               // Prior: moderately regularizing, centered on no preference (log-odds = 0)
  for (n in 1:N) {
    choice[n] ~ bernoulli_logit(alpha[treat[n]]);     // Likelihood: group-specific logistic regression on choice outcome
  }
} 
