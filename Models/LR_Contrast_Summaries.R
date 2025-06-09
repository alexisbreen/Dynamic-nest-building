################################################################################################################################################################################

#Script for first choice and all choices model contrast summaries for the manuscript

#Dynamic strategic social learning in nest-building zebra finches and its generalisability

#Code authored by Alexis J Breen (alexis_breen@eva.mpg.de) 

################################################################################################################################################################################

#
##
###Housekeeping
##
#

#Load required libraries

library(rethinking)

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

# Define a helper function to compute mean and HPDI
get_summary <- function(samples, prob = 0.89) {
  h <- HPDI(samples, prob = prob)
  data.frame(
    Mean = mean(samples),
    HPDI_Lower = h[1],
    HPDI_Upper = h[2]
  )
}

#Apply to each contrast in first choice model
summary_contrasts1 <- lapply(contrasts1, get_summary)
summary_contrasts1 <- do.call(rbind, summary_contrasts1)
summary_contrasts1$Contrast <- rownames(summary_contrasts1)
summary_contrasts1$Model <- "First Choice"

#Apply to each contrast in all choices model
summary_contrasts2 <- lapply(contrasts2, get_summary)
summary_contrasts2 <- do.call(rbind, summary_contrasts2)
summary_contrasts2$Contrast <- rownames(summary_contrasts2)
summary_contrasts2$Model <- "All Choices"

#Combine into one table
contrast_summary_table <- rbind(summary_contrasts1, summary_contrasts2)

#Reorder columns
contrast_summary_table <- contrast_summary_table[, c("Model", "Contrast", "Mean", "HPDI_Lower", "HPDI_Upper")]

#Round for readability 
contrast_summary_table <- contrast_summary_table %>%
  mutate(across(c(Mean, HPDI_Lower, HPDI_Upper), \(x) round(x, 2)))

#View the table
print(contrast_summary_table)
