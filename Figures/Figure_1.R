################################################################################################################################################################################

#Script for plotting Figure 1 for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Required packages

library(rethinking) 
library(dplyr)

#Load data

SSL_d <- read.csv(file.choose(), header = TRUE) #File: SSL_Data_Test_Processed.csv

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

#Define standard error function

se <- function(x) sd(x) / sqrt(length(x))

#Now compute group-level mean and SE per treatment

mean_df <- aggregate(choice ~ treat, data = ind_means, FUN = mean)
se_df <- aggregate(choice ~ treat, data = ind_means, FUN = se)

#Merge

mean_se_result <- merge(mean_df, se_df, by = "treat", suffixes = c("_mean", "_se"))

#Extract samples

s_FC <- extract.samples(FC_m)
s_AC <- extract.samples(AC_m)

#Probabilities

p <- plogis(s_FC$alpha)
s_AC$alpha <- s_AC$a
p2 <- plogis(s_AC$alpha)

#Contrasts

contrasts1 <- list(
  `Mat-Sat` = p[,2] - p[,1],
  `Inc-Sat` = p[,2] - p[,3],
  `Inc-Diss` = p[,2] - p[,4]
)

contrasts2 <- list(
  `Mat-Sat` = p2[,2] - p2[,1],
  `Inc-Sat` = p2[,2] - p2[,3],
  `Inc-Diss` = p2[,2] - p2[,4]
)

#Index and labels

n <- length(contrasts1)
labels <- names(contrasts1)

#
##
###Plot Figure 1
##
#

pdf(file = "Figure1.pdf", height = 7, width = 10) #If want pdf of plot

#General layout

layout(matrix(c(1,1,1,1,2,3,4,5), nrow = 2, byrow = TRUE), heights = c(1.25, 2))
par(mar = c(4, 5, 5, 2))

#
##
###Row 1: Behavioural data
##
#

#General plot setup

plot(NULL, xlim = c(0,104), ylim = c(0,1),  xaxs="i", xaxt="n", yaxt="n", xlab = NA, ylab = NA)
mtext("Behaviour",  cex = 1, side = 3, line = 3, font = 2)
mtext("Choose social",  cex = 1, side = 2, line = 3)
mtext("Choice",  cex = 1, side = 1, line = 3)
axis(1, at = c(1:25), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(27:51), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(53:77), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(1, at = c(79:103), labels = c("1","","","","","","","","","","","","","","","","","","","","","","","","25"), cex = 1)
axis(2, at = seq(0, .9, by = 0.1), labels = c(expression("0%"),"","","","","","","","",""), cex.axis = 1)
axis(2, at = 1, labels = expression("100%"), cex.axis = 1)
mtext(expression(paste(bold(a))), side = 3, cex = 1, line = 2, adj = 0)
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

#Add text for first-choice percentage choose social 

text(x = 1.6, y = .3 + .1, expression("30%"))    #Mat-Sat
text(x = 27.6, y = .63 + .1, expression("63%"))  #Mat-Diss
text(x = 53.6, y = .14 + .15, expression("14%")) #Inc-Sat
text(x = 79.6, y = .20 + .16, expression("20%")) #Inc-Diss

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
###Row 1, Plot 1: First choice estimates
##
#

#General plot setup

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 4.7), xlab = "Prob. choose social", cex.lab = 1.5, ylab = "", yaxt = "n")
axis(2, at = 4:1, labels = c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss"), cex = 1.5, las = 1)
mtext(expression(paste(bold(b))), side = 3, cex = 1, line = 2, adj = 0)
text(x = 0.5, y = 5.1, labels = "First choice", xpd = NA, cex = 1.5)

#Plot densities

for(i in 1:4){
  
  #Required objects
  
  row <- 5 - i
  d <- density(p[,i])        
  x <- d$x
  y <- d$y / max(d$y) * 0.6 #Re-scaling so no overlap
  
  #Set colours & line types
  
  fill_col <- if (i %in% c(3, 4)) col.alpha("gold", 0.2) else col.alpha(rangi2, 0.2)
  border_col <- if (i %in% c(3, 4)) "gold" else rangi2
  border_lty <- if (i %in% c(2, 4)) 2 else 1  
  
  #Density plots + means + HPDI
  
  polygon(c(x, rev(x)), c(rep(row, length(x)), row + rev(y)),
          col = fill_col,
          border = border_col,
          lty = border_lty, lwd = 2)
  h <- HPDI(p[,i], prob = 0.89)
  segments(h[1], row, h[2], row, lwd = 2, col = "black")
  points(mean(p[,i]), row, pch = 19, col = "red")
  
}

#
##
###Row 2, Plot 2: First choice contrasts
##
#

#General plot setup

plot(NULL, xlim = c(-0.5, 1), ylim = c(0.5, n + 0.75), xlab = "Posterior difference", cex.lab = 1.5, ylab = "", yaxt = "n")
axis(2, at = n:1, labels = labels, las = 1)
text(x = .25, y = 4.1, labels = "Contrasts", xpd = NA, cex = 1.5)
text(x = 0.95, y = n + 0.75, labels = "Mat-Diss vs.", adj = 1, cex = 1)

#Plot densities

for(i in seq_along(contrasts1)){ 
  
  #Required objects
  contrast <- contrasts1[[i]]
  d <- density(contrast)
  x <- d$x
  y <- d$y / max(d$y) * 0.6
  row <- n - i + 1
  
  #Density plots + means + HPDI
  
  polygon(c(x, rev(x)), c(rep(row, length(x)), row + rev(y)),
          col = col.alpha("black", 0.2), border = FALSE)
  h <- HPDI(contrast, prob = 0.89)
  segments(h[1], row, h[2], row, lwd = 2, col = "black")
  points(mean(contrast), row, pch = 19, col = "red")
  
}

#Line at 0

abline(v = 0, lty = 3)

#
##
###Row 2, Plot 3: All choices estimates
##
#

#General plot setup

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 4.7), xlab = "Prob. choose social", cex.lab = 1.5, ylab = "", yaxt = "n")
axis(2, at = 4:1, labels = c("Mat-Sat", "Mat-Diss", "Inc-Sat", "Inc-Diss"), cex = 1.5, las = 1)
text(x = 0.5, y = 5.1, labels = "All choices", xpd = NA, cex = 1.5)

#Plot densities

for(i in 1:4){
  
  #Required objects
  
  row <- 5 - i
  d <- density(p2[,i])        
  x <- d$x
  y <- d$y / max(d$y) * 0.6
  
  #Set colours & line types
  
  fill_col <- if (i %in% c(3, 4)) col.alpha("gold", 0.2) else col.alpha(rangi2, 0.2)
  border_col <- if (i %in% c(3, 4)) "gold" else rangi2
  border_lty <- if (i %in% c(2, 4)) 2 else 1  
  
  #Density plots + means + HPDI
  
  polygon(c(x, rev(x)), c(rep(row, length(x)), row + rev(y)),
          col = fill_col,
          border = border_col,
          lty = border_lty, lwd = 2)
  h <- HPDI(p[,i], prob = 0.89)
  segments(h[1], row, h[2], row, lwd = 2, col = "black")
  points(mean(p[,i]), row, pch = 19, col = "red")
  
}

#
##
###Row 2, Plot 4: All choices contrasts
##
#

#General plot setup

plot(NULL, xlim = c(-0.5, 1), ylim = c(0.5, n + 0.75), xlab = "Posterior difference", cex.lab = 1.5, ylab = "", yaxt = "n")
axis(2, at = n:1, labels = labels, las = 1)
text(x = .25, y = 4.1, labels = "Contrasts", xpd = NA, cex = 1.5)
text(x = 0.95, y = n + 0.75, labels = "Mat-Diss vs.", adj = 1, cex = 1)

#Plot densities

for(i in seq_along(contrasts2)){
  
  #Required objects
  
  contrast <- contrasts2[[i]]
  d <- density(contrast)
  x <- d$x
  y <- d$y / max(d$y) * 0.6
  row <- n - i + 1
  
  #Density plots + means + HPDI
  
  polygon(c(x, rev(x)), c(rep(row, length(x)), row + rev(y)),
          col = col.alpha("black", 0.2), border = FALSE)
  h <- HPDI(contrast, prob = 0.89)
  segments(h[1], row, h[2], row, lwd = 2, col = "black")
  points(mean(contrast), row, pch = 19, col = "red")
}

#Line at 0

abline(v = 0, lty = 3)

dev.off() #Close PDF file

