################################################################################################################################################################################

#Script for plotting Figure 3 for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Simulation 1: What happens when only choose social rewards?

#First, run Baseline_Post_Study_Sim.R

d_sim_1 <- PH_base_sim_fct(N_sim = 1, N_choosers = 10000, N_soc_pay = 1, N_non_soc_pay = 0, SS_matched = 0) 

#Simulation 2: What happens when choose social rewards 2x choose asocial?

d_sim_2 <- PH_base_sim_fct(N_sim = 1, N_choosers = 10000, N_soc_pay = 2, N_non_soc_pay = 1, SS_matched = 0) 

#Empty matrices to hold Simulation 2 & 3 calculations...

sim_1 <- matrix(0,4,25)
sim_2 <- matrix(0,4,25)

#Perform calculations

for(i in 1:25){
  
  sim_1[1,i] <- sapply(i, function(x) mean(d_sim_1$copy[d_sim_1$trial == x & d_sim_1$treat == 1])) #Mat-Sat
  sim_1[2,i] <- sapply(i, function(x) mean(d_sim_1$copy[d_sim_1$trial == x & d_sim_1$treat == 2])) #Mat-Diss
  sim_1[3,i] <- sapply(i, function(x) mean(d_sim_1$copy[d_sim_1$trial == x & d_sim_1$treat == 3])) #Inc-Sat
  sim_1[4,i] <- sapply(i, function(x) mean(d_sim_1$copy[d_sim_1$trial == x & d_sim_1$treat == 4])) #Inc-Diss
  
  sim_2[1,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 1])) #Mat-Sat
  sim_2[2,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 2])) #Mat-Diss
  sim_2[3,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 3])) #Inc-Sat
  sim_2[4,i] <- sapply(i, function(x) mean(d_sim_2$copy[d_sim_2$trial == x & d_sim_2$treat == 4])) #Inc-Diss
  
}

#Collapse matrices together into a matrix-specific data frame, and add custom trial count that skips 26, 52, & 78 to allow for correct spacing in graph along x-axis

sim_1 <- as.data.frame(list(means = c(sim_1[1,], sim_1[2,], sim_1[3,], sim_1[4,]), count = c(1:25, 27:51, 53:77, 79:103)))
sim_2 <- as.data.frame(list(means = c(sim_2[1,], sim_2[2,], sim_2[3,], sim_2[4,]), count = c(1:25, 27:51, 53:77, 79:103)))

#
##
###Plot Figure 3
##
#

pdf(file = "Figure3.pdf", height = 5, width = 10)

#General plot layout

layout(matrix(c(1,1,1,1,
                2,2,2,2), nrow = 2, byrow = TRUE))
par(mar = c(3, 5, 1, 2), oma = c(4, .5, 4, 0))


#
##
###Row 1: Social pays
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext(expression("Choose " * italic("social pays")), cex = 1, side = 2, line = 3)
mtext("Simulations",  cex = 1, side = 3, line = 2, font = 2)
axis(1, at = 1:25, labels = c("1", rep("", 23), "25"))
axis(1, at = 27:51, labels = c("1", rep("", 23), "25"))
axis(1, at = 53:77, labels = c("1", rep("", 23), "25"))
axis(1, at = 79:103, labels = c("1", rep("", 23), "25"))
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Draw mean social material-use for Simulation 1

lines(sim_1[1:25, 2], sim_1[1:25, 1], col = rangi2, lwd = 2)
lines(sim_1[26:50, 2], sim_1[26:50, 1], col = rangi2, lty = 2, lwd = 2)
lines(sim_1[51:75, 2], sim_1[51:75, 1], col = "gold", lwd = 2)
lines(sim_1[76:100, 2], sim_1[76:100, 1], col = 'gold', lty = 2, lwd = 2)

#Fill area under each treatment line

polygon( #Mat-Sat
  x = c(sim_1[1:25, 2], rev(sim_1[1:25, 2])),
  y = c(rep(0, 25), rev(sim_1[1:25, 1])),
  col = adjustcolor(rangi2, alpha.f = 0.3),
  border = NA
)

polygon( #Mat-Diss
  x = c(sim_1[26:50, 2], rev(sim_1[26:50, 2])),
  y = c(rep(0, 25), rev(sim_1[26:50, 1])),
  col = adjustcolor(rangi2, alpha.f = 0.3),
  border = NA
)

polygon( #Inc-Sat
  x = c(sim_1[51:75, 2], rev(sim_1[51:75, 2])),
  y = c(rep(0, 25), rev(sim_1[51:75, 1])),
  col = adjustcolor("gold", alpha.f = 0.3),
  border = NA
)

polygon( #Inc-Diss
  x = c(sim_1[76:100, 2], rev(sim_1[76:100, 2])),
  y = c(rep(0, 25), rev(sim_1[76:100, 1])),
  col = adjustcolor("gold", alpha.f = 0.3),
  border = NA
)

#Add treatment labels

label_xs <- c(13, 39, 65, 91)
label_names <- c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss")

#Plot text above 

for (i in 1:4) {
  text(x = label_xs[i], y = 1.15, labels = label_names[i], cex = 1.5, xpd = NA)
}

#
##
###Row 2: Social pays more
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext(expression("Choose " * italic("social pays more")), cex = 1, side = 2, line = 3)
mtext("Choice",  cex = 1, side = 1, line = 3)
axis(1, at = 1:25, labels = c("1", rep("", 23), "25"))
axis(1, at = 27:51, labels = c("1", rep("", 23), "25"))
axis(1, at = 53:77, labels = c("1", rep("", 23), "25"))
axis(1, at = 79:103, labels = c("1", rep("", 23), "25"))
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
abline(v = c(26,52,78))

#Draw mean social material-use for Simulation 2

lines(sim_2[1:25, 2], sim_2[1:25, 1], col = rangi2, lwd = 2)
lines(sim_2[26:50, 2], sim_2[26:50, 1], col = rangi2, lty = 2, lwd = 2)
lines(sim_2[51:75, 2], sim_2[51:75, 1], col = "gold", lwd = 2)
lines(sim_2[76:100, 2], sim_2[76:100, 1], col = 'gold', lty = 2, lwd = 2)

#Fill area under each treatment line

polygon( #Mat-Sat
  x = c(sim_2[1:25, 2], rev(sim_2[1:25, 2])),
  y = c(rep(0, 25), rev(sim_2[1:25, 1])),
  col = adjustcolor(rangi2, alpha.f = 0.3),
  border = NA
)

polygon( #Mat-Diss
  x = c(sim_2[26:50, 2], rev(sim_2[26:50, 2])),
  y = c(rep(0, 25), rev(sim_2[26:50, 1])),
  col = adjustcolor(rangi2, alpha.f = 0.3),
  border = NA
)

polygon( #Inc-Sat
  x = c(sim_2[51:75, 2], rev(sim_2[51:75, 2])),
  y = c(rep(0, 25), rev(sim_2[51:75, 1])),
  col = adjustcolor("gold", alpha.f = 0.3),
  border = NA
)

polygon( #Inc-Diss
  x = c(sim_2[76:100, 2], rev(sim_2[76:100, 2])),
  y = c(rep(0, 25), rev(sim_2[76:100, 1])),
  col = adjustcolor("gold", alpha.f = 0.3),
  border = NA
)

dev.off() #Close PDF file
