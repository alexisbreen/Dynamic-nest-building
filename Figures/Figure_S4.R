#####################################################################################################################################################

#Script to run Figure S4 for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

#####################################################################################################################################################

#
##
###Housekeeping
##
#

#Required packages

library(rethinking) 
library(dplyr)

#Load data

SSL_d <- read.csv(file.choose(), header = TRUE)

#Convert choice to binary

SSL_d <- SSL_d %>% mutate(choice = ifelse(choice == 1, 1, 0))

# Mean observed choice for each treatment at each trial

obs <- matrix(0, 4, 25)
for(i in 1:25){
  obs[1,i] <- mean(SSL_d$choice[SSL_d$trial == i & SSL_d$Experiment == 1 & SSL_d$Satisfaction == 1]) # MS
  obs[2,i] <- mean(SSL_d$choice[SSL_d$trial == i & SSL_d$Experiment == 1 & SSL_d$Satisfaction == 2]) # MD
  obs[3,i] <- mean(SSL_d$choice[SSL_d$trial == i & SSL_d$Experiment == 2 & SSL_d$Satisfaction == 1]) # IS
  obs[4,i] <- mean(SSL_d$choice[SSL_d$trial == i & SSL_d$Experiment == 2 & SSL_d$Satisfaction == 2]) # ID
}

#Flatten for plotting

obs <- as.data.frame(list(
  real = c(obs[1,], obs[2,], obs[3,], obs[4,]),
  count = c(1:25, 27:51, 53:77, 79:103)
))

#Mean and SE observed choice for social for each treatment across trials

#First, compute individual-level mean choice across trials

ind_means <- aggregate(choice ~ id + treat, data = SSL_d, FUN = mean)

#Convert to percentage

ind_means$choice <- ind_means$choice * 100

#Run 2 simulations... 
#Note this will overwrite objects from Figure_3.R, if already run

N <- 100

#EWA Baseline model

#First, run the Baseline_Post_Study_Sim.R code 

d_sim_1 <- PH_base_sim_fct(N_sim = N, N_soc_pay = 1, N_non_soc_pay = 1, SS_matched = 1) 

#EWA Monotonic model

#Now, run the Monotonic_Post_Study_Sim.R code

d_sim_2 <- PH_sim_mono_fct(N_sim = N, N_soc_pay = 1, N_non_soc_pay = 1, SS_matched = 1) 

#Empty matrices to hold Simulation 1 calculations...

sim_1_S1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Satisfied-construction
sim_1_D1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Dissatisfied-construction
sim_1_S2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Satisfied-reproduction
sim_1_D2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Dissatisfied-reproduction

sim_2_S1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Satisfied-construction
sim_2_D1 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Dissatisfied-construction
sim_2_S2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Satisfied-reproduction
sim_2_D2 <- matrix(0, 25, N) #Proportion choosing social across choices 1 - 25 and simulations 1 - N in Dissatisfied-reproduction

#Perform calculations...

for(choice in 1:25){
  for(sim in 1:N){
    sim_1_S1[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 1 & d_sim_1$sat_level == 1])
    sim_1_D1[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 1 & d_sim_1$sat_level == 2])
    sim_1_S2[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 2 & d_sim_1$sat_level == 1])
    sim_1_D2[choice, sim] <- mean(d_sim_1$copy[d_sim_1$trial == choice & d_sim_1$sim == sim & d_sim_1$experiment == 2 & d_sim_1$sat_level == 2])
    
    sim_2_S1[choice, sim] <- mean(d_sim_2$copy[d_sim_2$trial == choice & d_sim_2$sim == sim & d_sim_2$experiment == 1 & d_sim_2$sat_level == 1])
    sim_2_D1[choice, sim] <- mean(d_sim_2$copy[d_sim_2$trial == choice & d_sim_2$sim == sim & d_sim_2$experiment == 1 & d_sim_2$sat_level == 2])
    sim_2_S2[choice, sim] <- mean(d_sim_2$copy[d_sim_2$trial == choice & d_sim_2$sim == sim & d_sim_2$experiment == 2 & d_sim_2$sat_level == 1])
    sim_2_D2[choice, sim] <- mean(d_sim_2$copy[d_sim_2$trial == choice & d_sim_2$sim == sim & d_sim_2$experiment == 2 & d_sim_2$sat_level == 2])
  }
}

#Combined together into a data frame and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

sim_1 <- as.data.frame(rbind(sim_1_S1, sim_1_D1, sim_1_S2, sim_1_D2))
sim_1$count <- c(1:25, 27:51, 53:77, 79:103)

sim_2 <- as.data.frame(rbind(sim_2_S1, sim_2_D1, sim_2_S2, sim_2_D2))
sim_2$count <- c(1:25, 27:51, 53:77, 79:103)

#
##
###Plot Figure S4
##
#

# Start outputting to a PDF file for the plot grid
pdf(file = "Figure_S4.pdf", height = 7, width = 10)

#Overall plot spacing and layout set-up

layout(matrix(c(1,1,1,1,
                2,2,2,2,
                3,3,3,3), nrow = 3, byrow = TRUE))
par(mar = c(5, 5, 4, 2), oma = c(4, .5, 4, 0))

#
##
###Row 1: Behavioural data
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Choose social",  cex = 1, side = 2, line = 3)
axis(1, at = 1:25, labels = c("1", rep("", 23), "25"))
axis(1, at = 27:51, labels = c("1", rep("", 23), "25"))
axis(1, at = 53:77, labels = c("1", rep("", 23), "25"))
axis(1, at = 79:103, labels = c("1", rep("", 23), "25"))
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Observed social material choice

lines(obs[1:25, 2], obs[1:25, 1], col = rangi2, lwd = 2)              #Mat-Sat
lines(obs[26:50, 2], obs[26:50, 1], col = rangi2, lty = 2, lwd = 2)   #Mat-Diss
lines(obs[51:75, 2], obs[51:75, 1], col = "gold", lwd = 2)            #Inc-Sat
lines(obs[76:100, 2], obs[76:100, 1], col = "gold", lty = 2, lwd = 2) #Inc-Diss

#Fill area under each treatment line

polygon(c(obs[1:25, 2], rev(obs[1:25, 2])), #Mat-Sat
        c(rep(0, 25), rev(obs[1:25, 1])),
        col = col.alpha(rangi2, 0.2), border = NA)

polygon(c(obs[26:50, 2], rev(obs[26:50, 2])), #Mat-Diss
        c(rep(0, 25), rev(obs[26:50, 1])),
        col = col.alpha(rangi2, 0.2), border = NA)

polygon(c(obs[51:75, 2], rev(obs[51:75, 2])), #Inc-Sat
        c(rep(0, 25), rev(obs[51:75, 1])),
        col = col.alpha("gold", 0.2), border = NA)

polygon(c(obs[76:100, 2], rev(obs[76:100, 2])), #Inc-Diss
        c(rep(0, 25), rev(obs[76:100, 1])),
        col = col.alpha("gold", 0.2), border = NA)

#Add treatment labels
#Center x positions of each block

label_xs <- c(13, 39, 65, 91)
label_names <- c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss")

# Plot text above the box

for(i in 1:4){
  text(x = label_xs[i], y = 1.15, labels = label_names[i], cex = 1.5, xpd = NA)
}

#
##
###Row 2: Social pays
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Choose social",  cex = 1, side = 2, line = 3)
axis(1, at = 1:25, labels = c("1", rep("", 23), "25"))
axis(1, at = 27:51, labels = c("1", rep("", 23), "25"))
axis(1, at = 53:77, labels = c("1", rep("", 23), "25"))
axis(1, at = 79:103, labels = c("1", rep("", 23), "25"))
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Draw mean social material-use trajectories per simulation (of 100)

for(j in 1:N) lines(sim_1[1:25,(N+1)], sim_1[1:25,j], col = col.alpha(rangi2, alpha = .05))
for(j in 1:N) lines(sim_1[26:50,(N+1)], sim_1[26:50,j], col = col.alpha(rangi2, alpha = .05))
for(j in 1:N) lines(sim_1[51:75,(N+1)], sim_1[51:75,j], col = col.alpha("gold", alpha = .05))
for(j in 1:N) lines(sim_1[76:100,(N+1)], sim_1[76:100,j], col = col.alpha("gold", alpha = .05))

#Draw mean social material-use trajectories across 100 simulations

lines(1:25, apply(sim_1[1:25, 1:N], 1 , mean), col = rangi2, lwd = 2)
lines(27:51, apply(sim_1[26:50, 1:N], 1 , mean), col = rangi2, lty = 2, lwd = 2)
lines(53:77, apply(sim_1[51:75, 1:N], 1 , mean), col = "gold", lwd = 2)
lines(79:103, apply(sim_1[76:100, 1:N], 1 , mean), col = 'gold', lty = 2, lwd = 2)

#
##
###Row 3: Social pays more
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Choose social",  cex = 1, side = 2, line = 3)
axis(1, at = 1:25, labels = c("1", rep("", 23), "25"))
axis(1, at = 27:51, labels = c("1", rep("", 23), "25"))
axis(1, at = 53:77, labels = c("1", rep("", 23), "25"))
axis(1, at = 79:103, labels = c("1", rep("", 23), "25"))
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Draw mean social material-use trajectories per simulation (of 100)

for(j in 1:N) lines(sim_2[1:25,(N+1)], sim_2[1:25,j], col = col.alpha(rangi2, alpha = .05))
for(j in 1:N) lines(sim_2[26:50,(N+1)], sim_2[26:50,j], col = col.alpha(rangi2, alpha = .05))
for(j in 1:N) lines(sim_2[51:75,(N+1)], sim_2[51:75,j], col = col.alpha("gold", alpha = .05))
for(j in 1:N) lines(sim_2[76:100,(N+1)], sim_2[76:100,j], col = col.alpha("gold", alpha = .05))

#Draw mean social material-use trajectories across 100 simulations

lines(1:25, apply(sim_2[1:25, 1:N], 1 , mean), col = rangi2, lwd = 2)
lines(27:51, apply(sim_2[26:50, 1:N], 1 , mean), col = rangi2, lty = 2, lwd = 2)
lines(53:77, apply(sim_2[51:75, 1:N], 1 , mean), col = "gold", lwd = 2)
lines(79:103, apply(sim_2[76:100, 1:N], 1 , mean), col = "gold", lty = 2, lwd = 2)

#Add common outer labels

mtext("Choice", side = 1, line = 1.75, outer = TRUE, cex = 1)
mtext("Approximation", side = 3, line = -1, outer = TRUE, font = 2, cex = 1)

dev.off() #Close PDF file
