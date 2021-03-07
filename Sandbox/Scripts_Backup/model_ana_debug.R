library(ggplot2)
library(dplyr)
rm(list = ls())

setwd("~/Documents/IC/CMEECourseWork/MiniProject/Code")

# Loading dataset
quad_info <- read.csv("../Results/respirationrate/fitInfos/quad_info.csv")
cubic_info <- read.csv("../Results/respirationrate/fitInfos/cubic_info.csv")
briere_info <- read.csv("../Results/respirationrate/fitInfos/Briere_info.csv")
school_info <- read.csv("../Results/respirationrate/fitInfos/Schoolfield_info.csv")

id <- quad_info$quad.id
quad_aicc <- quad_info$quad.aic
cubic_aicc <- cubic_info$cubic.aic
briere_aicc <- briere_info$bri.aic
school_aicc <- school_info$sch.aic

aic_info <- data.frame(cbind(id, quad_aicc, cubic_aicc, briere_aicc, school_aicc))

# To define a funtion of calculating numbers of best-fit occurence.
# ------------------------------------------------------------------------------------
best_fit_number <- function(x = aic_info[1,], num = 4, threshold = 4){
  minIdx <- which.min(x[1:num])
  minSort <- sort(x[1:num], decreasing = F)
  x[num + minIdx] <- x[num + minIdx] + 1
  for (i in 2:num){
    if (abs(minSort[1] - minSort[i]) < threshold){
      x[match(minSort[i],x) + num] <- x[match(minSort[i],x) + num] + 1
    }
  }
  return(c(x[(num+1):(num+4)]))
}

colnames(aic_info) <- c('id', "Quadratic", "Cubic", 'Briere', "Schoolfield_High_Temp")
aic_info$quad_best <- 0
aic_info$cubic_best <- 0
aic_info$bri_best <- 0
aic_info$sch_best <- 0
aic_info[is.na(aic_info)] <- 9999.9999

# Applying the algorithm above
# ------------------------------------------------------------------------------------
aic_info[6:9] <- t(apply(aic_info[2:9], 1, FUN = best_fit_number, num = 4, threshold = 4))
aic_info$best_fit_aic <- colnames(aic_info[2:9])[apply(aic_info[2:5], 1, which.min)]

bestFitOccur <- data.frame(Model_Name=factor(colnames(aic_info)[2:5], levels = colnames(aic_info)[2:5]),
                           Best_Fit_Count=c(sum(aic_info$quad_best), 
                                            sum(aic_info$cubic_best), 
                                            sum(aic_info$bri_best),
                                            sum(aic_info$sch_best)))


# # Create pie chart
# value <- c(perc_quad_aic, perc_cubic_aic, perc_briere_aic, perc_school_aic)
# models <- c("Quadratic Model", "Cubic Model", "Briere Model", "Schoolfield with High Temperature")
# pie_d <- data.frame(Models = models, Percentage = value)
# 
# # Compute the position of labels
# pie_d <- pie_d %>% 
#   arrange(desc(Models)) %>%
#   mutate(prop = Percentage / sum(pie_d$Percentage) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )
# 
# # Plot pie chart
# pdf("../Results/Fitting_info/bestfit_pie_aic.pdf")
# ggplot(pie_d, aes(x = "", y = Percentage, fill = Models)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y", start = 0) +
#   theme_void() +
#   geom_text(aes(y = ypos, label = paste(Percentage, "%", sep = "")), color = "white", size=4) +
#   scale_fill_brewer(palette="Set1")
# dev.off()
#################################################################################################