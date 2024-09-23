#####################################################################################################################################################

#Script to run Figure S2 for the manuscript

#Strategic social learning in avian nest construction and potentially beyond

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

#####################################################################################################################################################

#
##
###Housekeeping
##
#

#The following script must have already been run in order to fully execute the current script:

#SL_Pre_Study_Simulation_Fxn.R 

#Required packages:

library(rethinking) 

#
##
###Figure pre-processing
##
#

#First, run pre-study simulation 1 and 2

#Pre-study simulation 1: low use of social material across treatments

PS_sim_1 <- PS_sim_fct( 
  
  N_sim = 10,
  sim_N = 1,
  real_phi = rbind(c(0.1, 0.1), c(0.1, 0.1)), #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_lambda = rbind(c(5, 5), c(5, 5)),      #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_rho = rbind(c(0, 0), c(0, 0)),         #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_sigma = rbind(c(.1, .1), c(.1, .1))    #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2)  
  
) 

#Pre-study simulation 2: high use of social material across treatments

PS_sim_2 <- PS_sim_fct( 
  
  N_sim = 10,
  sim_N = 2,
  real_phi = rbind(c(0.4, 0.4), c(0.4, 0.4)), #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_lambda = rbind(c(5, 5), c(5, 5)),      #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_rho = rbind(c(4, 4), c(4, 4)),         #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2) 
  real_sigma = rbind(c(.7, .7), c(.7, .7))    #(satisfied exp1, satisfied exp2),(dissatisfied exp1, dissatisfied exp2)  
  
) 

#Now, for each simulation, calculate proportion choosing social material across choices

#Empty matrices to hold Simulation 1 calculations...

pre_sim_1_S1 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-construction
pre_sim_1_D1 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-construction
pre_sim_1_S2 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-reproduction
pre_sim_1_D2 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-reproduction

pre_sim_2_S1 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-construction
pre_sim_2_D1 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-construction
pre_sim_2_S2 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-reproduction
pre_sim_2_D2 <- matrix(0, 25, 10) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-reproduction

#Perform calculations...

for(choice in 1:25){
  for(sim in 1:10){
    pre_sim_1_S1[choice, sim] <- mean(PS_sim_1$copy[PS_sim_1$trial == choice & PS_sim_1$sim == sim & PS_sim_1$experiment == 1 & PS_sim_1$sat_level == 1])
    pre_sim_1_D1[choice, sim] <- mean(PS_sim_1$copy[PS_sim_1$trial == choice & PS_sim_1$sim == sim & PS_sim_1$experiment == 1 & PS_sim_1$sat_level == 2])
    pre_sim_1_S2[choice, sim] <- mean(PS_sim_1$copy[PS_sim_1$trial == choice & PS_sim_1$sim == sim & PS_sim_1$experiment == 2 & PS_sim_1$sat_level == 1])
    pre_sim_1_D2[choice, sim] <- mean(PS_sim_1$copy[PS_sim_1$trial == choice & PS_sim_1$sim == sim & PS_sim_1$experiment == 2 & PS_sim_1$sat_level == 2])
  
    pre_sim_2_S1[choice, sim] <- mean(PS_sim_2$copy[PS_sim_2$trial == choice & PS_sim_2$sim == sim & PS_sim_2$experiment == 1 & PS_sim_2$sat_level == 1])
    pre_sim_2_D1[choice, sim] <- mean(PS_sim_2$copy[PS_sim_2$trial == choice & PS_sim_2$sim == sim & PS_sim_2$experiment == 1 & PS_sim_2$sat_level == 2])
    pre_sim_2_S2[choice, sim] <- mean(PS_sim_2$copy[PS_sim_2$trial == choice & PS_sim_2$sim == sim & PS_sim_2$experiment == 2 & PS_sim_2$sat_level == 1])
    pre_sim_2_D2[choice, sim] <- mean(PS_sim_2$copy[PS_sim_2$trial == choice & PS_sim_2$sim == sim & PS_sim_2$experiment == 2 & PS_sim_2$sat_level == 2])
  }
}

#Combine together into a data frame and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

pre_sim_1 <- as.data.frame(rbind(pre_sim_1_S1, pre_sim_1_D1, pre_sim_1_S2, pre_sim_1_D2))
pre_sim_1$count <- c(1:25, 27:51, 53:77, 79:103)
pre_sim_2 <- as.data.frame(rbind(pre_sim_2_S1, pre_sim_2_D1, pre_sim_2_S2, pre_sim_2_D2))
pre_sim_2$count <- c(1:25, 27:51, 53:77, 79:103)

#Next, write function that feeds all low and all high simulations one-by-one into our EWA model,
#runs our model, and extracts both attraction estimates and posterior predictions for the
#social material across all material choices per simulant

run_models <- function( sim_type, simulants ){
  
  for(sim in 1:max(simulants)){ 
    
    #Determine if simulation 1 or 2
    
    if(sim_type == "low"){
      sim_df <- PS_sim_1[(PS_sim_1$sim == sim),]
    } else {
      sim_df <- PS_sim_2[(PS_sim_2$sim == sim),]
    }
    
    #Make list for stan EWA model
    
    sim_list <- list(
      
      N = nrow(sim_df),
      N_id = length(unique(sim_df$id)),
      id = sim_df$id,
      manip = sim_df$experiment,
      sat = sim_df$sat_level,
      trial = sim_df$trial,
      choice = sim_df$choice,
      cum_soc = sim_df$cum_soc,
      cum_non_soc = sim_df$cum_non_soc,
      social = sim_df$social,
      touch_P = sample(10:100, 47, replace = TRUE),
      touch_O = sample(10:100, 47, replace = TRUE),
      log_dur = log(runif(47, 14400, (14400 * 4))) 
      
      )
    
    #Run EWA model
    
    mod <- cstan(file = "SSL_EWA_Model.stan", data = sim_list, cores = 4, chains = 4, refresh = 10, iter = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
    
    #Extract posterior from EWA model
    
    s_mod <- extract.samples(mod)
    
    #Empty matrices to hold...
    
    A_mu_mod <- matrix(0, 4, 25) #Mean attraction estimates for social material
    pp_mod <- matrix(0, 4, 25)   #Mean posterior predicted choice for social material
  
    #Fill empty matrices
    
    for(i in 1:25){
      
    A_mu_mod[1,i] <- mean(s_mod$Atx_soc[,i,1:10])  #Satisfied-construction
    A_mu_mod[2,i] <- mean(s_mod$Atx_soc[,i,11:18]) #Dissatisfied-construction
    A_mu_mod[3,i] <- mean(s_mod$Atx_soc[,i,19:32]) #Satisfied-reproduction
    A_mu_mod[4,i] <- mean(s_mod$Atx_soc[,i,33:47]) #Dissatisfied-reproduction
      
    pp_mod[1,i] <- mean(s_mod$soc_pp[,i,1:10])  #Satisfied-construction
    pp_mod[2,i] <- mean(s_mod$soc_pp[,i,11:18]) #Dissatisfied-construction
    pp_mod[3,i] <- mean(s_mod$soc_pp[,i,19:32]) #Satisfied-reproduction
    pp_mod[4,i] <- mean(s_mod$soc_pp[,i,33:47]) #Dissatisfied-reproduction
      
    }
    
    #Collapse matrices together into a matrix-specific data frame, and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis
    
    mod_A_mu <- as.data.frame(list(sim_type = ifelse(sim_type == "low", "low", "high"), simulant = sim, means = c(A_mu_mod[1,], A_mu_mod[2,], A_mu_mod[3,], A_mu_mod[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
    mod_pp <- as.data.frame(list(sim_type = ifelse(sim_type == "low", "low", "high"), simulant = sim, pred = c(pp_mod[1,], pp_mod[2,], pp_mod[3,], pp_mod[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
    
    #Saving into one big data frame
    
    if(sim == 1){
      
      mods_A_mu <- mod_A_mu
      mods_pp <- mod_pp
      
    } else {
      
      mods_A_mu <- rbind(mods_A_mu, mod_A_mu)
      mods_pp <- rbind(mods_pp, mod_pp)
      
    }
    
  }
  
  #Return multiple data frames in a list
  
  return(list(mods_A_mu, mods_pp))
  
}

#Execute the above function for all simulants in simulation 1 & 2
#NOTE: 10 SIMULANTS TAKES ALMOST 12 HOURS TO RUN ON EACH LOW AND HIGH SCRIPT!!!

low_sims <- run_models(sim_type = "low", simulants = 10)
high_sims <- run_models(sim_type = "high", simulants = 10)

#Put data frames in list into their own object

df_low_sims_A <- low_sims[[1]]
df_high_sims_A <- high_sims[[1]]
df_low_sims_pp <- low_sims[[2]]
df_high_sims_pp <- high_sims[[2]]

#Finally, write plotting function to examine each simulant & their associated EWA model output

plot_simulant <- function(df_A_low, df_A_high, df_pp_low, df_pp_high, df_sim_low, df_sim_high){
  
  #Plot set-up for shading
  
  plot(NULL, xlim = c(0,104), ylim = c(0,1), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
  corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  par(xpd = TRUE) #Draw outside plot area
  text(x = corners[2]+1, y = mean(corners[4] - .15), "High", srt = 270)
  text(x = corners[2]+1, y = mean(corners[3] + .15), "Low", srt = 270)
  par(xpd = FALSE) #End draw outside plot area
  rect(0,-1,26,1.2, col = "white", border = FALSE)
  rect(26,-1,52,1.2, col = "#f0f0f0", border = FALSE)
  rect(52,-1,78,1.2, col = "#d9d9d9", border = FALSE)
  rect(78,-1,104,1.2, col = "#bdbdbd", border = FALSE)
  
  #Plot set-up for data
  
  par(new = TRUE)
  plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
  if(i == 5) mtext("Mean percentage",  cex = 1, side = 2, line = 3)
  if(i == 10) mtext("Choice",  cex = 1, side = 1, line = 3)
  axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
  axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
  axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
  axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
  axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
  axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
  abline(v = c(26,52,78))
  
  #Underlay mean attraction estimates
  
  #Low sims
  
  points(df_A_low$count, df_A_low$means, col = "red", pch = 16)
  lines(df_A_low[1:25, 4], df_A_low[1:25, 3], col = "red", lty = 1)
  lines(df_A_low[26:50, 4], df_A_low[26:50, 3], col = "red", lty = 1)
  lines(df_A_low[51:75, 4], df_A_low[51:75, 3], col = "red", lty = 1)
  lines(df_A_low[76:100, 4], df_A_low[76:100, 3], col = "red", lty = 1)
  
  #High sims
  
  points(df_A_high$count, df_A_high$means, col = "red", pch = 16)
  lines(df_A_high[1:25, 4], df_A_high[1:25, 3], col = "red", lty = 1)
  lines(df_A_high[26:50, 4], df_A_high[26:50, 3], col = "red", lty = 1)
  lines(df_A_high[51:75, 4], df_A_high[51:75, 3], col = "red", lty = 1)
  lines(df_A_high[76:100, 4], df_A_high[76:100, 3], col = "red", lty = 1)
  
  #Overlay observed and predicted prop. social-material use
  
  #Predicted - low sims
  
  points(df_pp_low$count, df_pp_low$pred, col = rangi2, pch = 20)
  lines(df_pp_low[1:25, 4], df_pp_low[1:25, 3], col = rangi2)
  lines(df_pp_low[26:50, 4], df_pp_low[26:50, 3], col = rangi2)
  lines(df_pp_low[51:75, 4], df_pp_low[51:75, 3], col = rangi2)
  lines(df_pp_low[76:100, 4], df_pp_low[76:100, 3], col = rangi2)
  
  #Predicted - high sims
  
  points(df_pp_high$count, df_pp_high$pred, col = rangi2, pch = 20)
  lines(df_pp_high[1:25, 4], df_pp_high[1:25, 3], col = rangi2)
  lines(df_pp_high[26:50, 4], df_pp_high[26:50, 3], col = rangi2)
  lines(df_pp_high[51:75, 4], df_pp_high[51:75, 3], col = rangi2)
  lines(df_pp_high[76:100, 4], df_pp_high[76:100, 3], col = rangi2)
  
  #Plot mean social material-use trajectories per simulation 
  
  #Low sims
  
  points(pre_sim_1[,11], pre_sim_1[,i], col = "black", pch = 1)
  lines(pre_sim_1[1:25,(10+1)], pre_sim_1[1:25,i], col = "black")
  lines(pre_sim_1[26:50,(10+1)], pre_sim_1[26:50,i], col = "black")
  lines(pre_sim_1[51:75,(10+1)], pre_sim_1[51:75,i], col = "black")
  lines(pre_sim_1[76:100,(10+1)], pre_sim_1[76:100,i], col = "black")
  
  #High sims
  
  points(pre_sim_2[,11], pre_sim_2[,i], col = "black", pch = 1)
  lines(pre_sim_2[1:25,(10+1)], pre_sim_2[1:25,i], col = "black")
  lines(pre_sim_2[26:50,(10+1)], pre_sim_2[26:50,i], col = "black")
  lines(pre_sim_2[51:75,(10+1)], pre_sim_2[51:75,i], col = "black")
  lines(pre_sim_2[76:100,(10+1)], pre_sim_2[76:100,i], col = "black")
  
}

#
##
###PLOT Figure S2
##
#

#pdf(file = "Figure_S2.pdf", height = 12, width = 10) #If want pdf of plot

#General plot space setup

par(mfrow = c(11,1), mar = c(1,4,1,1), oma = c(3,1,2,1))

#Plot set-up for shading

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
text(x = corners[2]+1, y = mean(corners[4] - .15), "High", srt = 270)
text(x = corners[2]+1, y = mean(corners[3] + .15), "Low", srt = 270)
par(xpd = FALSE) #End draw outside plot area
rect(0,-1,26,1.2, col = "white", border = FALSE)
rect(26,-1,52,1.2, col = "#f0f0f0", border = FALSE)
rect(52,-1,78,1.2, col = "#d9d9d9", border = FALSE)
rect(78,-1,104,1.2, col = "#bdbdbd", border = FALSE)

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0,104), ylim = c(0,1), xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Satisfied-construction",  adj = .04, side = 3, line = 1)
mtext("Dissatisfied-construction",  adj = .34, side = 3, line = 1)
mtext("Satisfied-reproduction",  adj = .65, side = 3, line = 1)
mtext("Dissatisfied-reproduction",  adj = .975, side = 3, line = 1)
axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Simulation 1

#Draw mean social material-use trajectories per simulation 

for(j in 1:10) lines(pre_sim_1[1:25,(10+1)], pre_sim_1[1:25,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_1[26:50,(10+1)], pre_sim_1[26:50,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_1[51:75,(10+1)], pre_sim_1[51:75,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_1[76:100,(10+1)], pre_sim_1[76:100,j], col = col.alpha("black", alpha = .2))

#Draw mean social material-use trajectories across simulations

lines(1:25, apply(pre_sim_1[1:25, 1:10], 1 , mean), col = "black")
lines(27:51, apply(pre_sim_1[26:50, 1:10], 1 , mean), col = "black")
lines(53:77, apply(pre_sim_1[51:75, 1:10], 1 , mean), col = "black")
lines(79:103, apply(pre_sim_1[76:100, 1:10], 1 , mean), col = "black")

#Overlay mean social material-use values across simulation

for(i in 1:25) points(i, apply(pre_sim_1[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 26:50) points(i + 1, apply(pre_sim_1[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 51:75) points(i + 2, apply(pre_sim_1[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 76:100) points(i + 3, apply(pre_sim_1[i, 1:10], 1 , mean), col = "black", pch = 20)

#Simulation 2

#Draw mean social material-use trajectories per simulation 

for(j in 1:10) lines(pre_sim_2[1:25,(10+1)], pre_sim_2[1:25,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_2[26:50,(10+1)], pre_sim_2[26:50,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_2[51:75,(10+1)], pre_sim_2[51:75,j], col = col.alpha("black", alpha = .2))
for(j in 1:10) lines(pre_sim_2[76:100,(10+1)], pre_sim_2[76:100,j], col = col.alpha("black", alpha = .2))

#Draw mean social material-use trajectories across simulations

lines(1:25, apply(pre_sim_2[1:25, 1:10], 1 , mean), col = "black")
lines(27:51, apply(pre_sim_2[26:50, 1:10], 1 , mean), col = "black")
lines(53:77, apply(pre_sim_2[51:75, 1:10], 1 , mean), col = "black")
lines(79:103, apply(pre_sim_2[76:100, 1:10], 1 , mean), col = "black")

#Overlay mean social material-use values across simulation

for(i in 1:25) points(i, apply(pre_sim_2[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 26:50) points(i + 1, apply(pre_sim_2[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 51:75) points(i + 2, apply(pre_sim_2[i, 1:10], 1 , mean), col = "black", pch = 20)
for(i in 76:100) points(i + 3, apply(pre_sim_2[i, 1:10], 1 , mean), col = "black", pch = 20)

#Plot each simulant's choice for social material against the EWA model's estimated 
#social material attraction & posterior predicted choice for social material

for(i in 1:10){
  
  plot_simulant(df_A_low = df_low_sims_A[(df_low_sims_A$simulant == i),], 
              df_A_high = df_high_sims_A[(df_high_sims_A$simulant == i),], 
              df_pp_low = df_low_sims_pp[(df_low_sims_pp$simulant == i),], 
              df_pp_high = df_high_sims_pp[(df_high_sims_pp$simulant == i),], 
              df_sim_low = pre_sim_1[, c(i, 11)], #11 b/c it's the count column
              df_sim_high = pre_sim_2[, c(i, 11)]) #11 b/c it's the count column
}

#dev.off() #If want pdf of plot


