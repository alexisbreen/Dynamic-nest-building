#####################################################################################################################################################

#Script to run Figure 2 for the manuscript

#Dynamic strategic social learning in avian nest construction and potentially beyond

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

#####################################################################################################################################################

#
##
###Housekeeping
##
#

#The following scripts must have already been run in order to fully execute the current script:

#SSL_EWA_Model_Execution.R
#SL_Post_Study_Simulation_Fxn.R 

#Required packages:
library(rethinking) 
library(png)

#
##
###Pre-graph processing
##
#

#Load data if not already in global environment (file name: SSL_Test_Data_Processed.csv)
#SSL_d <- read.csv(file.choose(), header = T)

#Transform choice variable to calculate proportions b/c currently have 1 & 2 as choice-options due to STAN model dependent-variable coding requirements

SSL_d <- SSL_d %>% mutate(choice = ifelse(choice == 1, 1, 0))

#Empty matrices to hold calculations for each treatment across choices...

obs <- matrix(0, 4, 25)  #Mean observed choice for social material
A_mu <- matrix(0, 4, 25) #Mean attraction estimates for social material
pp <- matrix(0, 4, 25)   #Mean posterior predicted choice for social material

#Perform calculations

for(i in 1:25){
  
  obs[1,i] <- sapply(i, function(x) mean(SSL_d$choice[SSL_d$trial == x & SSL_d$Experiment == 1 & SSL_d$Satisfaction == 1])) #Satisfied-construction
  obs[2,i] <- sapply(i, function(x) mean(SSL_d$choice[SSL_d$trial == x & SSL_d$Experiment == 1 & SSL_d$Satisfaction == 2])) #Dissatisfied-construction
  obs[3,i] <- sapply(i, function(x) mean(SSL_d$choice[SSL_d$trial == x & SSL_d$Experiment == 2 & SSL_d$Satisfaction == 1])) #Satisfied-reproduction
  obs[4,i] <- sapply(i, function(x) mean(SSL_d$choice[SSL_d$trial == x & SSL_d$Experiment == 2 & SSL_d$Satisfaction == 2])) #Dissatisfied-reproduction
  
  A_mu[1,i] <- mean(s$Atx_soc[,i,1:10])  #Satisfied-construction
  A_mu[2,i] <- mean(s$Atx_soc[,i,11:18]) #Dissatisfied-construction
  A_mu[3,i] <- mean(s$Atx_soc[,i,19:32]) #Satisfied-reproduction
  A_mu[4,i] <- mean(s$Atx_soc[,i,33:47]) #Dissatisfied-reproduction
  
  pp[1,i] <- mean(s$soc_pp[,i,1:10])  #Satisfied-construction
  pp[2,i] <- mean(s$soc_pp[,i,11:18]) #Dissatisfied-construction
  pp[3,i] <- mean(s$soc_pp[,i,19:32]) #Satisfied-reproduction
  pp[4,i] <- mean(s$soc_pp[,i,33:47]) #Dissatisfied-reproduction
  
}

#Collapse matrices together into a matrix-specific data frame, and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

obs <- as.data.frame(list(real = c(obs[1,], obs[2,], obs[3,], obs[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
A_mu <- as.data.frame(list(means = c(A_mu[1,], A_mu[2,], A_mu[3,], A_mu[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
pp <- as.data.frame(list(pred = c(pp[1,], pp[2,], pp[3,], pp[4,]), count = c(1:25, 27:51, 53:77, 79:103)))

#Run 3 simulations... 

#Simulation 1: Can we approximate our data?

N <- 10 #Number of simulations matched to treatment sample size - 10 is not actually that many sims, so graph will vary more here, than in sims below; play around if you want!
d_sim_1 <- PH_sim_fct(N_sim = N, N_soc_pay = 1, N_non_soc_pay = 1, SS_matched = 1) 

#Simulation 2: What happens when only choose social rewards?

d_sim_2 <- PH_sim_fct(N_sim = 1, N_choosers = 10000, N_soc_pay = 1, N_non_soc_pay = 0, SS_matched = 0) 

#Simulation 3: What happens when choose social rewards 2x choose asocial?

d_sim_3 <- PH_sim_fct(N_sim = 1, N_choosers = 10000, N_soc_pay = 2, N_non_soc_pay = 1, SS_matched = 0) 

#Empty matrices to hold Simulation 1 calculations...

sim_1_S1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-construction
sim_1_D1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-construction
sim_1_S2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Satisfied-reproduction
sim_1_D2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - 10 in Dissatisfied-reproduction

#Perform calculations...

for(choice in 1:25){
  for(sim in 1:N){
      sim_1_S1[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 1 & d_sim_1$sat_level == 1])
      sim_1_D1[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 1 & d_sim_1$sat_level == 2])
      sim_1_S2[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 2 & d_sim_1$sat_level == 1])
      sim_1_D2[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 2 & d_sim_1$sat_level == 2])
  }
}

#Combined together into a data frame and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

sim_1 <- as.data.frame(rbind(sim_1_S1, sim_1_D1, sim_1_S2, sim_1_D2))
sim_1$count <- c(1:25, 27:51, 53:77, 79:103)

#Empty matrices to hold Simulation 2 & 3 calculations...

sim_2 <- matrix(0,4,25)
sim_3 <- matrix(0,4,25)

#Perform calculations

for(i in 1:25){
  
  sim_2[1,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 1])) #Satisfied-construction
  sim_2[2,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 2])) #Dissatisfied-construction
  sim_2[3,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 3])) #Satisfied-reproduction
  sim_2[4,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 4])) #Dissatisfied-reproduction
  
  sim_3[1,i] <- sapply(i, function(x) mean(d_sim_3$copy[d_sim_3$trial == x & d_sim_3$treat == 1])) #Satisfied-construction
  sim_3[2,i] <- sapply(i, function(x) mean(d_sim_3$copy[d_sim_3$trial == x & d_sim_3$treat == 2])) #Dissatisfied-construction
  sim_3[3,i] <- sapply(i, function(x) mean(d_sim_3$copy[d_sim_3$trial == x & d_sim_3$treat == 3])) #Satisfied-reproduction
  sim_3[4,i] <- sapply(i, function(x) mean(d_sim_3$copy[d_sim_3$trial == x & d_sim_3$treat == 4])) #Dissatisfied-reproduction
  
}

#Collapse matrices together into a matrix-specific data frame, and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

sim_2 <- as.data.frame(list(means = c(sim_2[1,], sim_2[2,], sim_2[3,], sim_2[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
sim_3 <- as.data.frame(list(means = c(sim_3[1,], sim_3[2,], sim_3[3,], sim_3[4,]), count = c(1:25, 27:51, 53:77, 79:103)))

#
##
###PLOT Figure 2
##
#

#pdf(file = "Figure2.pdf", height = 12, width = 10) #If want pdf of plot

#Overall plot spacing and layout set-up

par(mfrow = c(7,1), mar = c(1,4,1,2), oma = c(3,3,3,0))
layout(
  matrix(c(1,2,3,4,5,6,7), 
  nrow = 7, ncol = 1, byrow = TRUE),
  heights = c(2,1,1,1,1,2,2)
)

#
##
###Row 1 - BEHAVIOUR: MATERIAL CHOICE & EWA MODEL FIT
##
#

#Row-specific spacing

par(mar = c(7,5,0,4))

#Plot set-up for shading

plot(NULL, xlim = c(0,104), ylim = c(0,1), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
rect(0,-1,26,1.5, col = "white", border = FALSE)
rect(26,-1,52,1.5, col = "#f0f0f0", border = FALSE)
rect(52,-1,78,1.5, col = "#d9d9d9", border = FALSE)
rect(78,-1,104,1.5, col = "#bdbdbd", border = FALSE)

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Behaviour: material choice and EWA model fit",  cex = 1, side = 3, line = 1, font = 2)
mtext("Mean percentage",  cex = 1, side = 2, line = 3)
mtext("Choice",  cex = 1, side = 1, line = 3)
axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
mtext(expression(paste(bold(a))), side = 3, cex = 1, line = 1, adj = 0)
abline(v = c(26,52,78))
legend(x = 1, y = 1.05, legend = c("observed choose social", "predicted choose social", expression(paste("attraction ", italic(A), " for social"))), lty = 1, pch = c(1,20,16), col = c("black", rangi2, "red"),  bty = "n", cex = 1)

#Underlay mean attraction estimates

points(A_mu$count, A_mu$means, col = "red", pch = 16)
lines(A_mu[1:25, 2], A_mu[1:25, 1], col = "red", lty = 1)
lines(A_mu[26:50, 2], A_mu[26:50, 1], col = "red", lty = 1)
lines(A_mu[51:75, 2], A_mu[51:75, 1], col = "red", lty = 1)
lines(A_mu[76:100, 2], A_mu[76:100, 1], col = "red", lty = 1)

#Overlay observed and predicted prop. social-material use

#Predicted

points(pp$count, pp$pred, col = rangi2, pch = 20)
lines(pp[1:25, 2], pp[1:25, 1], col = rangi2)
lines(pp[26:50, 2], pp[26:50, 1], col = rangi2)
lines(pp[51:75, 2], pp[51:75, 1], col = rangi2)
lines(pp[76:100, 2], pp[76:100, 1], col = rangi2)

#Observed

points(obs$count, obs$real, col = "black", pch = 1)
lines(obs[1:25, 2], obs[1:25, 1], col = "black")
lines(obs[26:50, 2], obs[26:50, 1], col = "black")
lines(obs[51:75, 2], obs[51:75, 1], col = "black")
lines(obs[76:100, 2], obs[76:100, 1], col = "black")

#Add text for first-choice percentage choose social 

text(x = 1.6, y = .3 - .2, expression("30%"))
text(x = 27.6, y = .63 + .1, expression("63%"))
text(x = 53.6, y = .14 - .1, expression("14%"))
text(x = 79.6, y = .20 - .1, expression("20%"))

#
##
###Row 2 - COGNITION: EWA MODEL PARAMETER ESTIMATES
##
#

#PHI

#For x-axis

x <- (1:4)

#Panel-specific spacing

par(mar = c(1,5,0,30))

#Calculate mean & HPDI - also used in Table S1

mu_phi <- apply(s_phi, 2, mean)
HPDI_phi <- apply(s_phi, 2, HPDI)

#Plot set-up for shading

plot(NULL, xlim = c(0.5,4.5), ylim = c(0,1), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
rect(0,-1,1.5,1.5,col = "white")
rect(1.5,-1,2.5,1.5,col = "#f0f0f0")
rect(2.5,-1,3.5,1.5,col = "#d9d9d9")
rect(3.5,-1,5.5,1.5,col = "#bdbdbd")

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0.5,4.5), ylim = c(0,1), xaxs="i", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
axis(2, at = seq(0, .9, by = 0.1), labels = c(0,"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = 1, cex.axis = 1)
axis(side = 2, at = .5, labels = expression(paste("updating ", italic(phi))), col = "white", cex.axis = 1.5, line = 2)
for(i in 1:100){
  points(jitter(1:4, 1), s_phi[i, 1:4], col = rangi2) 
} 
arrows(x0 = x, y0 = HPDI_phi[1,], x1 = x, y1 = HPDI_phi[2,], length = 0, col = "black", lwd = 2)
points(x = x, y = mu_phi, col = "red", pch = 19)
mtext("Cognition: EWA model parameter estimates",  cex = 1, side = 3, line = 1, font = 2)
mtext(expression(paste(bold(b))), side = 3, cex = 1, line = 1, adj = 0)
mtext("Asocial-", at = -.1, side = 2, cex = 1, line = 5)
text(x = .6, y = .95, labels = expression(paste(italic("(i)"))))

#LAMBDA

#Panel-specific spacing

par(mar = c(1,5,0,30))

#Calculate mean & HPDI - also used in Table S1

mu_lambda <- apply(s_lambda, 2, mean)
HPDI_lambda <- apply(s_lambda, 2, HPDI)

#Plot set-up for shading

plot(NULL, xlim = c(0.5,4.5), ylim = c(0,7), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
rect(0,-1,1.5,7.5,col = "white")
rect(1.5,-1,2.5,7.5,col = "#f0f0f0")
rect(2.5,-1,3.5,7.5,col = "#d9d9d9")
rect(3.5,-1,5.5,7.5,col = "#bdbdbd")

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0.5,4.5), ylim = c(0,7), xaxs="i", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
axis(2, at = seq(0, 6, by = 1), labels = c(0,"","","","","",""), cex.axis = 1)
axis(2, at = 7, labels = 7, cex.axis = 1)
axis(side = 2, at = 3, labels = expression(paste("sensitivity ", italic(lambda))), col = "white", cex.axis = 1.5, line = 2)
for(i in 1:100){
  points(jitter(1:4, 1), s_lambda[i, 1:4], col = rangi2) 
} 
arrows(x0 = x, y0 = HPDI_lambda[1,], x1 = x, y1 = HPDI_lambda[2,], length = 0, col = "black", lwd = 2)
points(x = x, y = mu_lambda, col = "red", pch = 19)
text(x = .6, y = 6.5, labels = expression(paste(italic("(ii)"))))

#RHO

#Panel-specific spacing

par(mar = c(1,5,0,30))

#Calculate mean & HPDI - also used in Table S1

mu_rho <- apply(s_rho, 2, mean)
HPDI_rho <- apply(s_rho, 2, HPDI)

#Plot set-up for shading

plot(NULL, xlim = c(0.5,4.5), ylim = c(0,10), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
rect(0,-1,1.5,10.5,col = "white")
rect(1.5,-1,2.5,10.5,col = "#f0f0f0")
rect(2.5,-1,3.5,10.5,col = "#d9d9d9")
rect(3.5,-1,5.5,10.5,col = "#bdbdbd")

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0.5,4.5), ylim = c(0,10), xaxs="i", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
axis(2, at = seq(0, 9, by = 1), labels = c(0,"","","","","","","","",""), cex.axis = 1)
axis(2, at = 10, labels = 10, cex.axis = 1)
axis(side = 2, at = 5, labels = expression(paste("reactivity ", italic(rho))), col = "white", cex.axis = 1.5, line = 2)
for(i in 1:100){
  points(jitter(1:4, 1), s_rho[i, 1:4], col = rangi2) 
} 
abline(h = 1, lty = 2)
arrows(x0 = x, y0 = HPDI_rho[1,], x1 = x, y1 = HPDI_rho[2,], length = 0, col = "black", lwd = 2)
points(x = x, y = mu_rho, col = "red", pch = 19)
mtext("Social-", at = -.1, side = 2, cex = 1, line = 5)
text(x = .6, y = 9.5, labels = expression(paste(italic("(iii)"))))

#SIGMA

#Panel-specific spacing

par(mar = c(1,5,0,30))

#Calculate mean & HPDI - also used in Table S1

mu_sig <- apply(s_sig, 2, mean)
HPDI_sig <- apply(s_sig, 2, HPDI)

#Plot set-up for shading

plot(NULL, xlim = c(0.5,4.5), ylim = c(0,1), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
rect(0,-1,1.5,1.5,col = "white")
rect(1.5,-1,2.5,1.5,col = "#f0f0f0")
rect(2.5,-1,3.5,1.5,col = "#d9d9d9")
rect(3.5,-1,5.5,1.5,col = "#bdbdbd")

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0.5,4.5), ylim = c(0,1), xaxs="i", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
axis(2, at = seq(0, .9, by = 0.1), labels = c(0,"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = 1, cex.axis = 1)
axis(side = 2, at = .5, labels = expression(paste("bias ", italic(sigma))), col = "white", cex.axis = 1.5, line = 2)
for(i in 1:100){
  points(jitter(1:4, 1), s_sig[i, 1:4], col = rangi2) 
} 
arrows(x0 = x, y0 = HPDI_sig[1,], x1 = x, y1 = HPDI_sig[2,], length = 0, col = "black", lwd = 2)
points(x = x, y = mu_sig, col = "red", pch = 19)
text(x = .6, y = .95, labels = expression(paste(italic("(iv)"))))

#
##
###Row 3
##
#

#Row-specific spacing

par(mar = c(4,5,3,4))

#Plot set-up for shading

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
rect(0,-1,26,1.5, col = "white", border = FALSE)
rect(26,-1,52,1.5, col = "#f0f0f0", border = FALSE)
rect(52,-1,78,1.5, col = "#d9d9d9", border = FALSE)
rect(78,-1,104,1.5, col = "#bdbdbd", border = FALSE)

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0,104), ylim = c(0,1), xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Replication: agent-based forward simulations",  cex = 1, side = 3, line = 1, font = 2)
mtext("Mean percentage",  cex = 1, side = 2, line = 3)
mtext("Choice",  cex = 1, side = 1, line = 3)
axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))
legend(x = 1, y = 1.05, legend = c("observed choose social", "simulated choose social"), lty = c(1,1), pch = c(1,20), col = c("black", rangi2),  bty = "n", cex = 1)
mtext(expression(paste(bold(c))), side = 3, cex = 1, line = 1, adj = 0)

#Simulation 1 versus observed data

#Draw mean social material-use trajectories per simulation 

for(j in 1:N) lines(sim_1[1:25,(N+1)], sim_1[1:25,j], col = col.alpha(rangi2, alpha = .2))
for(j in 1:N) lines(sim_1[26:50,(N+1)], sim_1[26:50,j], col = col.alpha(rangi2, alpha = .2))
for(j in 1:N) lines(sim_1[51:75,(N+1)], sim_1[51:75,j], col = col.alpha(rangi2, alpha = .2))
for(j in 1:N) lines(sim_1[76:100,(N+1)], sim_1[76:100,j], col = col.alpha(rangi2, alpha = .2))

#Draw mean social material-use trajectories across simulations

lines(1:25, apply(sim_1[1:25, 1:N], 1 , mean), col = rangi2)
lines(27:51, apply(sim_1[26:50, 1:N], 1 , mean), col = rangi2)
lines(53:77, apply(sim_1[51:75, 1:N], 1 , mean), col = rangi2)
lines(79:103, apply(sim_1[76:100, 1:N], 1 , mean), col = rangi2)

#Overlay mean social material-use values across simulation

for(i in 1:25) points(i, apply(sim_1[i, 1:N], 1 , mean), col = rangi2, pch = 20)
for(i in 26:50) points(i + 1, apply(sim_1[i, 1:N], 1 , mean), col = rangi2, pch = 20)
for(i in 51:75) points(i + 2, apply(sim_1[i, 1:N], 1 , mean), col = rangi2, pch = 20)
for(i in 76:100) points(i + 3, apply(sim_1[i, 1:N], 1 , mean), col = rangi2, pch = 20)

#Observed prop. use social material

points(obs$count, obs$real, col = "black", pch = 1)
lines(obs[1:25, 2], obs[1:25, 1], col = "black")
lines(obs[26:50, 2], obs[26:50, 1], col = "black")
lines(obs[51:75, 2], obs[51:75, 1], col = "black")
lines(obs[76:100, 2], obs[76:100, 1], col = "black")

#
##
###Row 4
##
#

#Row-specific spacing

par(mar = c(4,5,3,4))

#Plot set-up for shading

plot(NULL, xlim = c(0,104), ylim = c(0,1), axes = FALSE, xaxs="i", bty = "n", xlab = "", ylab = "") 
mtext("Implications: agent-based forward simulations",  cex = 1, side = 3, line = 1, font = 2)
rect(0,-1,26,1.5, col = "white", border = FALSE)
rect(26,-1,52,1.5, col = "#f0f0f0", border = FALSE)
rect(52,-1,78,1.5, col = "#d9d9d9", border = FALSE)
rect(78,-1,104,1.5, col = "#bdbdbd", border = FALSE)

#Plot set-up for data

par(new = TRUE)
plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Mean percentage",  cex = 1, side = 2, line = 3)
mtext("Choice",  cex = 1, side = 1, line = 3)
axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))
mtext(expression(paste(bold(d))), side = 3, cex = 1, line = 1, adj = 0)
legend(x = 1, y = 1.05, legend = c(expression(paste(italic("social pays"), " choose social")), expression(paste(italic("social pays double"), " choose social"))), lty = 1, pch = 20, col = c("black", rangi2),  bty = "n", cex = 1)

#Simulation 2 versus Simulation 3

#Draw mean social material-use for Simulation 2

points(sim_2$count, sim_2$means, col = "black", pch = 20)
lines(sim_2[1:25, 2], sim_2[1:25, 1], col = "black")
lines(sim_2[26:50, 2], sim_2[26:50, 1], col = "black")
lines(sim_2[51:75, 2], sim_2[51:75, 1], col = "black")
lines(sim_2[76:100, 2], sim_2[76:100, 1], col = "black")

#Draw mean social material-use for Simulation 3

points(sim_3$count, sim_3$means, col = rangi2, pch = 20)
lines(sim_3[1:25, 2], sim_3[1:25, 1], col = rangi2)
lines(sim_3[26:50, 2], sim_3[26:50, 1], col = rangi2)
lines(sim_3[51:75, 2], sim_3[51:75, 1], col = rangi2)
lines(sim_3[76:100, 2], sim_3[76:100, 1], col = rangi2)

#Legend

#Reset spacing for legend & render a new plot

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')

#Import png & align on plot (achieved via trial-and-error)

icons <- readPNG("SSL_icon.png")
rasterImage(icons, xleft = .3, ybottom = -.16, xright = .55, ytop = .58)

#Add legend text (alignment achieved via trial-and-error)

text(.63,.64,"Treatment", cex = 1.5, font = 2)
text(.75,.49,"Satisfied-construction", cex = 1.5)
text(.75,.31,"Dissatisfied-construction", cex = 1.5)
text(.75,.12,"Satisfied-reproduction", cex = 1.5)
text(.75,-.06,"Dissatisfied-reproduction", cex = 1.5)

#Add dashed box around legend (alignment achieved via trial-and-error)

rect(xleft = .27, ybottom = -.1775, xright = .97, ytop = .605, density = NULL,
col = NA, border = TRUE, lty = 2)

#dev.off() #If want pdf of plot
