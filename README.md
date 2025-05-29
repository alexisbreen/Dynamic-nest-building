# Dynamic-nest-construction

This GitHub repository hosts open materials for the manuscript

**Dynamic strategic social learning in nest-building zebra finches and its generalisability**

All code was authored by Alexis J. Breen (alexis_breen@eva.mpg.de) & Richard McElreath (‍richard_mcelreath@eva.mpg.de)

**Data Processing folder contains:**

- SSL_Data_Processing.R script to wrangle the raw data
- The original, raw data sheet: SSL_Data_Original.csv

**Data folder contains:**
 
- SSL_IMP_Data_Processed.csv produced from SSL_Data_Processing.R script & used for all analyses/graphing related to initial material preference
- SSL_IMP_Test_Data_Processed.csv produced from SSL_Data_Processing.R script & used for all analyses/graphing related to final material preference

**Figures folder contains:**

- Figure_1.R script 
- Figure_1.R script
- Figure_3.R script
- Figure_4.R script
- Figure_S1.R script 
- Figure_S1.R script
- Figure_S3.R script
- Figure_S4.R script

**Models folder contains:**

- LR_FC_Model.stan script expressing the defined hierarchical logistic regression model, examining the effect of treatment on first choice
- LR_AC_Model.stan script expressing the defined hierarchical logistic regression model, examining the effect of treatment on all choices, with trial-number and treatment slopes

- EWA_Baseline_Model.stan script expressing the defined non-time-varying multi-level experience-weighted attraction model, examining asocial and social influence on material choice
- EWA_Montonic_Model.stan script expressing the defined time-structured multi-level experience-weighted attraction model, examining asocial and social influence on material choice
- EWA_Model_Summaries.R script to summarise the estimate of both EWA baseline and monotonic models in a table

- Base_And_Mono_Pre_Study_Sim.R script to simulate data for EWA model validation checks
- Baseline_Post_Study_Sim.R script to simulate data from the EWA baseline model posterior
- Monotonic_Post_Study_Sim.R script to simulate data from the EWA monotonic model posterior

- Model_Execution.R script to run all stan models

**Software requirements:**

- Stan (for running the multi-level models): https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
- Rethinking (for processing fitted model outputs): https://github.com/rmcelreath/rethinking
- R (for running all code): https://www.rstudio.com/
