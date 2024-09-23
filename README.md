# Dynamic-nest-construction

This GitHub repository hosts open materials for the manuscript

**Dynamic strategic social learning in avian nest construction and potentially beyond**

and all code was authored by Alexis J. Breen (alexis_breen@eva.mpg.de) & Richard McElreath (‍richard_mcelreath@eva.mpg.de)

**Data Processing folder contains:**

- SSL_Data_Processing.R script to wrangle the raw data
- The original, raw data sheet: SSL_Data_Original.csv

**Data folder contains:**
 
- SSL_IMP_Data_Processed.csv produced from SSL_Data_Processing.R script & used for all analyses/graphing related to initial material preference
- SSL_IMP_Test_Data_Processed.csv produced from SSL_Data_Processing.R script & used for all analyses/graphing related to final material preference

**Figures folder contains:**

- SSL_Fig2.R script to reproduce Figure 2 in the main text; this includes running the post-study agent-based forward simulations
- SSL_Fig2_Icon.png used in Figure 2
- SSL_FigS1.R script to reproduce supplementary Figure 1 showing baseline material attractions
- SSL_FigS2.R script to reproduce supplementary Figure 2; this includes running the pre-study agent-based forward simulations

**Models folder contains:**

- SSL_EWA_Model.stan script expressing the defined multi-level experience-weighted attraction model, examining asocial and social influence on material choice
- SSL_EWA_Model_Execution.stan script to prepare data for, run, and post-process (e.g., extract posteriors) the learning parameter estimates generated from the EWA model
- SSL_Post_Study_Simulation_Fxn.R script to build the agent-based forward simulation function to test for replication plus effect(s) of varying reward-payoffs on 'choosers' generally
- SSL_Post_Study_Simulation_Fxn.R script to build the agent-based forward simulation function to validate our EWA model fit a priori

**Software requirements:**

- Stan (for running the multi-level models): https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
- Rethinking (for processing fitted model outputs): https://github.com/rmcelreath/rethinking
- R (for running all code): https://www.rstudio.com/
