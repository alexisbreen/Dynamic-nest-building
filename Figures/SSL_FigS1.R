#####################################################################################################################################################

#Script to run Figure S1 for the manuscript

#Strategic social learning in avian nest construction and potentially beyond

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

#####################################################################################################################################################

#
##
###Housekeeping
##
#

#The following script must have already been run in order to fully execute the current script:

#SSL_EWA_Model_Execution.R

#Required packages:

library(rethinking) 

#
##
###Figure pre-processing
##
#

#Empty vectors to fill with initial attraction to pink and orange string per bird

IP <- vector("list", 47)
IO <- vector("list", 47)

#Fill vectors

for(i in 1:47) {
  IP[[i]] <- density(s$A_init[,i,1])
  IO[[i]] <- density(s$A_init[,i,2])
}

#Function for stacked density plot, based on: https://stackoverflow.com/questions/25328533/overlapping-stacked-density-plots

stacked.density <- function(data, fac = 3, xlim, col = 'black', 
                             show.xaxis = T, xlab = '', ylab = ''){
  
  xvals = unlist(lapply(data, function(d) d$x))
  if(missing(xlim)) xlim=c(min(xvals), max(xvals))
  
  col = sapply(col, alpha)
  if(length(col) == 1) col = rep(col, length(data))
  
  plot(1, type = "n", xlim = xlim, ylim = c(1, length(data) + 1.5), cex.axis = .7,
       yaxt='n', bty='n', xaxt=ifelse(show.xaxis, 'l', 'n'), xlab = xlab, ylab = ylab)
  
  z = length(data):1
  for(i in 1:length(data)){
    d = data[[ z[i] ]]
    lines(d$x, fac*d$y + i, col = alpha(col[i], alpha = 1), lwd = 1)
    polygon(d$x, fac*d$y+ i, col = alpha(col[i], alpha = .3), border = NA)
    abline(h = i, lwd = 1)
    abline(v = 0, lty = 3, lwd = .25, col = "black")
  }
}

#
##
###Plot Figure S1
##
#

#pdf(file = "Figure_S1.pdf", height = 7, width = 5) #If want pdf of plot

#Setup plot space

par(mfrow = c(1,2), mar = c(2,2,1,1), oma = c(2,2,2,1))

#Pink string

stacked.density(rev(IP[1:47]), fac = .8, col=c("#ff007f"), show.xaxis = T)
mtext("pink string",  cex = 1, side = 3, line = 0.5, font = 1)
axis(2, at = seq(1, 47, by = 1), lwd = 0, labels = TRUE, cex.axis = .7)

#General labels

mtext("Initial attraction",  cex = 1, side = 3, line = 1.5, font = 2, adj = 3.5)
mtext("Material attraction [logit scale]",  cex = 1, side = 1, line = 2, adj = -1.15)
par(xpd = TRUE) #Draw outside plot area
mtext("Posterior density per bird",  cex = 1, side = 2, line = 2)
par(xpd = FALSE) #Draw outside plot area

#Orange string

stacked.density(rev(IO[1:47]), fac = .8, col=c("#FF5F1F"), show.xaxis = T)
mtext("orange string",  cex = 1, side = 3, line = 0.5, font = 1)
axis(2, at = seq(1, 47, by = 1), lwd = 0, labels = TRUE, cex.axis = .7)

#dev.off() #If want pdf of plot

