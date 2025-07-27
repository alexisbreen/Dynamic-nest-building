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

#Prep results output and add a total row
elpd_results <- rbind(
  elpd_results,
  data.frame(
    bird_id = "TOTAL",                                                                  #Label for total row
    elpd_baseline = sum(elpd_results$elpd_baseline, na.rm = TRUE),                      #Sum of baseline ELPD across birds
    elpd_monotonic = sum(elpd_results$elpd_monotonic, na.rm = TRUE),                    #Sum of monotonic ELPD
    delta_elpd = sum(elpd_results$delta_elpd, na.rm = TRUE)                             #Sum of delta ELPD (monotonic - baseline)
  )
)

#Subset k-values for target trials
k_baseline <- subset(k_df, model == "baseline" & trial >= L+1 & trial <= 24)            #Extract Pareto k for baseline model
k_monotonic <- subset(k_df, model == "monotonic" & trial >= L+1 & trial <= 24)          #Extract Pareto k for monotonic model

#Insert a blank row before TOTAL for visual gap
gap_row <- data.frame(
  bird_id = "GAP",                                                                      #Dummy label for spacing
  elpd_baseline = NA,                                                                   #Empty cells
  elpd_monotonic = NA,
  delta_elpd = NA
)

#Insert the gap row 
elpd_with_gap <- rbind(
  elpd_results[1:47, ],                                                                  #First 47 birds
  gap_row,                                                                               #Gap row
  elpd_results[48, ]                                                                     #Total row
)

#Bar labels & colour scheme
bar_labels <- as.character(elpd_with_gap$bird_id)                                        #Extract bird IDs for x-axis labels
bar_labels[bar_labels == "GAP"] <- ""                                                    #Blank label for gap
bar_labels[bar_labels == "TOTAL"] <- "ALL"                                               #Rename total to ALL

bar_colors <- ifelse(
  elpd_with_gap$delta_elpd > 0, "red",                                                   #Red if monotonic better
  ifelse(elpd_with_gap$delta_elpd < 0, "blue", NA)                                       #Blue if baseline better
)

#
##
###Plot Figure S4
##
#

# Start outputting to a PDF file for the plot grid

pdf(file = "Figure_S4.pdf", height = 8, width = 10)                                    

#Overall plot spacing and layout set-up

layout(
  matrix(c(1, 2,
           3, 3), nrow = 2, byrow = TRUE),                                            
  heights = c(1, 1)                                                                     
)

par(mar = c(5, 5, 3, 1))                                                                

#
##
###Row 1, top left: Baseline model predictive performance
#
#

plot(
  NA, xlim = c((L+1)-.5, 24.5), ylim = c(-2.5, 2.5),                                     
  xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "s",
  xlab = "Trial", ylab = expression("Pareto-" * italic(k)),
  main = "", frame.plot = TRUE
)
title("Predictive performance (baseline)", line = 1.5)                                  
abline(h = 0.7, lty = 2, col = "black")                                               
text(x = 6, y = 0.8, labels = expression(italic("k") == 0.7), adj = 0, cex = 0.7)     

for (b in unique(k_baseline$bird)) {
  points(subset(k_baseline, bird == b)[, c("trial", "k")], pch = 3, col = "blue")     

axis(1, at = 6:24, labels = ifelse(6:24 %in% c(6, 24), 6:24, ""), tick = TRUE)           
box()                                                                                   

#
##
###Row 1, top right: Montonic model predictive performance
#
#

plot(
  NA, xlim = c((L+1)-.5, 24.5), ylim = c(-2.5, 2.5),                                     
  xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "s",
  xlab = "Trial", ylab = expression("Pareto-" * italic(k)),
  main = "", frame.plot = TRUE
)
title("Predictive performance (monotonic)", line = 1.5)                                 
abline(h = 0.7, lty = 2, col = "black")                                               
text(x = 6, y = 0.8, labels = expression(italic("k") == 0.7), adj = 0, cex = 0.7)     

for (b in unique(k_monotonic$bird)) {
  points(subset(k_monotonic, bird == b)[, c("trial", "k")], pch = 3, col = "red")     

axis(1, at = 6:24, labels = ifelse(6:24 %in% c(6, 24), 6:24, ""), tick = TRUE)
box()

#
##
###Row 2: Baseline vs Montonic model predictive advantage
#
#

par(mar = c(5, 5, 5, 1))                                                                

x_positions <- 1:49                                                                      #X-axis position for each bar
x_lefts <- x_positions - 0.4                                                             #Left edge of bars
x_rights <- x_positions + 0.4                                                            #Right edge of bars
bar_vals <- elpd_with_gap$delta_elpd                                                     #Bar heights

plot(
  NULL,
  xlim = c(0, 50), ylim = c(-4.5, 4.5),
  xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "s",
  xlab = "Bird ID",
  ylab = expression("Difference in expected log point density (" * Delta * "ELPD)"),
  main = "", frame.plot = TRUE
)
title("Predictive advantage (monotonic - baseline)", line = 1.5)                         

for (i in seq_len(length(bar_vals))) {
  if (!is.na(bar_vals[i])) {
    rect(
      xleft = x_lefts[i],
      xright = x_rights[i],
      ybottom = 0,
      ytop = bar_vals[i],
      col = bar_colors[i],
      border = NA
    )
  }
}

abline(h = 0, lty = 2)                                                                   
box()

axis(1, at = c(1, 47, 49), labels = c(bar_labels[1], bar_labels[47], "ALL"), las = 1, tick = FALSE)  

#Add tick marks without GAP position
tick_positions <- x_positions[x_positions != 48]                                         
axis(1, at = tick_positions, labels = FALSE, tick = TRUE)                                

#Vertical gap divider
abline(v = 48, lty = 2)                                                                  

#Legend
legend(
  "top",
  inset = c(0, 0),
  horiz = TRUE,
  legend = c("monotonic > baseline", "baseline > monotonic"),
  fill = c("red", "blue"),
  border = NA,
  bty = "n",
  xpd = NA
)

dev.off()                                                                                 #Close PDF file