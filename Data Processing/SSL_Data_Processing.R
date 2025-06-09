################################################################################################################################################################################

#Data processing script for the manuscript

#Dynamic strategic social learning in nest-building birds and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de)

################################################################################################################################################################################

#Library packages to load

library(tidyverse)
library(rethinking)

#Load raw data 

setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))
raw_d <- read.csv("Data Processing/SSL_Data_Original.csv")

#First, prepare initial material preference measures

IMP_d <- raw_d[ , c("ID", "IMP_Sec", "IMP_P_TT_Count", "IMP_O_TT_Count")] #Select relevant columns
IMP_d <- IMP_d[!(IMP_d$ID == "51" | IMP_d$ID == "52" | IMP_d$ID == "53" | IMP_d$ID == "12" | IMP_d$ID == "30" | IMP_d$ID == "43"), ] #Remove birds dropped from study
rownames(IMP_d) <- NULL #Reset row number after dropping
IMP_d$id <- c(1:nrow(IMP_d)) #Define new continuous id number
IMP_d <- IMP_d[,-1] #Remove old continuous ID column i.e., column 1
IMP_d$IMP_Sec <- log(IMP_d$IMP_Sec) #Transform duration to log scale
IMP_d$IMP_Sec[is.na(IMP_d$IMP_Sec)] <- (-999) #Set NA values to negative number that will make Stan complain if evaluated
IMP_d$IMP_P_TT_Count[is.na(IMP_d$IMP_P_TT_Count)] <- (-999) #Set NA values to negative number that will make Stan complain if evaluated
IMP_d$IMP_O_TT_Count[is.na(IMP_d$IMP_O_TT_Count)] <- (-999) #Set NA values to negative number that will make Stan complain if evaluated

#Next, prepare final material preference measures - a bit of wrangling, so let's pipe. 

FMP_d <- raw_d %>% 
  dplyr::filter(ID != "51" & ID != "52" & ID != "53" & ID != "12" & ID != "30" & ID != "43") %>% #Remove birds dropped from study
  dplyr::select(ID, Experiment, Satisfaction, Dems_Col, IMP_Sec, D1:D25) %>% #Select target columns
  gather("deposit_N", "material", D1:D25) %>% #Reformat from wide to long
  mutate(material = ifelse(material == "p", 1, 2)) %>% #Rename material colour levels
  mutate(Dems_Col = ifelse(Dems_Col == "p", 1, 2)) %>% #Rename demonstrated colour levels
  mutate(Satisfaction = ifelse(Satisfaction == "s", 1, 2)) %>% #Reorder and rename satisfaction levels
  mutate(choice = ifelse(Dems_Col == material, 1, 2)) %>% #Determine choice (1 = social; 2 = non-social) based on whether demonstrated colour matched deposited colour
  arrange(ID) %>% #Group by bird ID
  mutate(id = rep(1:47, each = 25)) %>% #Add new continuous id after dropping birds
  mutate(trial = rep(1:25, len = 1175)) %>% #Add trial number
  mutate(treat = rep(1:4, times = c((10*25),(8*25),(14*25),(15*25)))) %>% #Add treatment variable based on sample sizes (1 = Material/Satisfied; 2 = Material/Dissatisfied; 3 = Breeding/Satisfied; 4 = Breeding/Dissatisfied)
  dplyr::select(-ID, -Dems_Col, -deposit_N, -material, -IMP_Sec) #Remove non-relevant columns


#Build variables (cum_soc, cum_non_soc & social) needed for Stan model
#These are calculated based on existing variables (vars) about material in the nestbox (nb) in the raw data
#To calculate, define function that subsets within each bird - there's probably a more elegant way!

nb <- function(x){
  cum_soc <- c() 
  cum_non_soc <- c()
  social <- c()
  ss <- subset(FMP_d, id == x)
  for(i in 1:nrow(ss)){
    if(ss$trial[i] > 1){
      cum_soc[i] <- length(which(ss$choice[1:(i - 1)] == 1)) #Cumulative count of social material in nestbox
      cum_non_soc[i] <- length(which(ss$choice[1:(i - 1)] == 2)) #Cumulative count of non-social material in nestbox
      social[i] <- ifelse(1 %in% ss$choice[2:i - 1], 1, 0) #Has the first social material desposit occurred? 1 = Yes; 2 = No
    } else {
      cum_soc[i] <- 0
      cum_non_soc[i] <- 0
      social[i] <- 0
    }
  }
  vars <- cbind(cum_soc,cum_non_soc,social)
  return(vars)
}

#Apply function across all unique bird ID

vars_ID <-as.data.frame(lapply(unique(FMP_d$id),nb))

#Wrangle the newly calculated variable data

vars_cs <- vars_ID %>% 
  dplyr::select(starts_with("cum_soc")) %>% 
  gather("bird","cum_soc", 1:47) #47 b/c 47 birds
vars_cns <- vars_ID %>% 
  dplyr::select(starts_with("cum_non_soc")) %>% 
  gather("bird","cum_non_soc", 1:47)
vars_s <- vars_ID %>% 
  dplyr::select(starts_with("social")) %>% 
  gather("bird","social", 1:47)

#Combine only target variables into data frame 

vars_all <- cbind(vars_cs,vars_cns,vars_s) %>% dplyr::select(cum_soc,cum_non_soc,social)

#Combine target variable data frame with existing data frame containing the rest of the data

FMP_d <- cbind(FMP_d,vars_all)

#To download clean initial material preference data sheet to csv
#write.csv(IMP_d,"~/SSL/SSL_IMP_Data_Processed.csv", row.names = FALSE)
#To download clean test data sheet to csv
#write.csv(FMP_d,"~/SSL/SSL_Test_Data_Processed.csv", row.names = FALSE)

