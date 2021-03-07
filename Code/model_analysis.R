## Script Name: model_analysis.R
## Author: Jingkai Sun
## Email Address: jingkai.sun20@imperial.ac.uk

## Things to be done
# - Read tables of fitting info from results directory and transform it as DataFrame
# - Merge ID, aic, bic, aicc, rsq columns as one table
# - Add best_fit column to find the model with the lowest aicc value
# - Calculate percentage of each model with lowest aicc value during fitting all IDs
# - Plot bar chart of the percentage

options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
rm(list = ls())

setwd("~/Documents/IC/ALL_Courses/CMEECourseWork/MiniProject/Code")

# Define functions
best_fit_count <- function(x, num, threshold){
  for (i in 1:num){
    if (x[i] <= threshold) {x[i + num] <- x[i + num] + 1}
  }
  return(c(x[(num+1):(num+4)]))
}

Delta <- function(x = aic_info[1,4:11], num = 4){
  minIdx <- which.min(x[1:num])
  minSort <- sort(x[1:num], decreasing = F)
  # x[num + minIdx] <- 0
  for (i in 2:num){
    delta <- minSort[i] - minSort[1]
    x[match(minSort[i],x) + num] <- x[match(minSort[i],x) + num] + delta
  }
  return(c(x[(num+1):(num+4)]))
}

fileName = "netphotosynthesisrate"
quad_info1 <- read.csv(paste("../Results/", fileName ,"/fitInfos/quad_info.csv", sep = ""))
cubic_info1 <- read.csv(paste("../Results/", fileName ,"/fitInfos/cubic_info.csv", sep = ""))
briere_info1 <- read.csv(paste("../Results/", fileName ,"/fitInfos/Briere_info.csv", sep = ""))
school_info1 <- read.csv(paste("../Results/", fileName ,"/fitInfos/Schoolfield_info.csv", sep = ""))
quad_info1$traits <- cubic_info1$traits <- "photosynthesis"
briere_info1$traits <- school_info1$traits <- "photosynthesis"

fileName = "respirationrate"
quad_info2 <- read.csv(paste("../Results/", fileName ,"/fitInfos/quad_info.csv", sep = ""))
cubic_info2 <- read.csv(paste("../Results/", fileName ,"/fitInfos/cubic_info.csv", sep = ""))
briere_info2 <- read.csv(paste("../Results/", fileName ,"/fitInfos/Briere_info.csv", sep = ""))
school_info2 <- read.csv(paste("../Results/", fileName ,"/fitInfos/Schoolfield_info.csv", sep = ""))
quad_info2$traits <- cubic_info2$traits <- "respiration"
briere_info2$traits <- school_info2$traits <- "respiration"

quad_info <- rbind(quad_info1, quad_info2)
cubic_info <- rbind(cubic_info1, cubic_info2)
briere_info <- rbind(briere_info1, briere_info2)
school_info <- rbind(school_info1, school_info2)

aic_info <- data.frame(cbind(quad_info$id, quad_info$n, quad_info$traits,
                             quad_info$quad.aic, cubic_info$cubic.aic, 
                             briere_info$bri.aic, school_info$sch.aic))

bic_info <- data.frame(cbind(quad_info$id, quad_info$n, quad_info$traits,
                             quad_info$quad.bic, cubic_info$cubic.bic, 
                             briere_info$bri.bic, school_info$sch.bic))

params_info <- data.frame(cbind(quad_info[1:5], cubic_info[3:6], 
                                briere_info[3:5], school_info[3:6]))

remove(quad_info, cubic_info, briere_info, school_info, quad_info1, cubic_info1, briere_info1,
       school_info1, quad_info2, cubic_info2, briere_info2, school_info2)

colnames(aic_info) <- c('id', "n", "traits", "Quadratic", "Cubic", 'Briere', "Schoolfield")
colnames(bic_info) <- c('id', "n", "traits", "Quadratic", "Cubic", 'Briere', "Schoolfield")

aic_info$Quadratic <- as.double(aic_info$Quadratic)
aic_info$Cubic <- as.double(aic_info$Cubic)
aic_info$Briere <- as.double(aic_info$Briere)
aic_info$Schoolfield <- as.double(aic_info$Schoolfield)
bic_info$Quadratic <- as.double(bic_info$Quadratic)
bic_info$Cubic <- as.double(bic_info$Cubic)
bic_info$Briere <- as.double(bic_info$Briere)
bic_info$Schoolfield <- as.double(bic_info$Schoolfield)

aic_info$quad_delta <- bic_info$quad_delta <- 0
aic_info$cubic_delta <- bic_info$cubic_delta <- 0
aic_info$bri_delta <- bic_info$bri_delta <- 0
aic_info$sch_delta <- bic_info$sch_delta <- 0
aic_info$quad_best <- bic_info$quad_best <- 0
aic_info$cubic_best <- bic_info$cubic_best <- 0
aic_info$bri_best <- bic_info$bri_best <- 0
aic_info$sch_best <- bic_info$sch_best <- 0
aic_info <- na.omit(aic_info)
bic_info <- na.omit(bic_info)

# Applying the algorithm above
# ------------------------------------------------------------------------------------
aic_info[8:11] <- t(apply(aic_info[4:11], 1, FUN = Delta, num = 4))
aic_info[12:15] <- t(apply(aic_info[8:15], 1, FUN = best_fit_count, num = 4, threshold = 2))
aic_info$best_fit_aic <- colnames(aic_info[4:11])[apply(aic_info[4:7], 1, which.min)]

bic_info[8:11] <- t(apply(bic_info[4:11], 1, FUN = Delta, num = 4))
bic_info[12:15] <- t(apply(bic_info[8:15], 1, FUN = best_fit_count, num = 4, threshold = 2))
bic_info$best_fit_aic <- colnames(bic_info[4:11])[apply(bic_info[4:7], 1, which.min)]

count_all <- c(sum(aic_info$quad_best), sum(aic_info$cubic_best), sum(aic_info$bri_best),
               sum(aic_info$sch_best), sum(bic_info$quad_best), sum(bic_info$cubic_best),
               sum(bic_info$bri_best), sum(bic_info$sch_best))

aic_info1 <- subset(aic_info, traits == "photosynthesis")
bic_info1 <- subset(bic_info, traits == "photosynthesis")
aic_info2 <- subset(aic_info, traits == "respiration")
bic_info2 <- subset(bic_info, traits == "respiration")

count1 <- c(sum(aic_info1$quad_best), sum(aic_info1$cubic_best), sum(aic_info1$bri_best),
               sum(aic_info1$sch_best), sum(bic_info1$quad_best), sum(bic_info1$cubic_best),
               sum(bic_info1$bri_best), sum(bic_info1$sch_best))

count2 <- c(sum(aic_info2$quad_best), sum(aic_info2$cubic_best), sum(aic_info2$bri_best),
            sum(aic_info2$sch_best), sum(bic_info2$quad_best), sum(bic_info2$cubic_best),
            sum(bic_info2$bri_best), sum(bic_info2$sch_best))

count_basic1 <- c(nrow(subset(aic_info1, best_fit_aic == "Quadratic")), 
                  nrow(subset(aic_info1, best_fit_aic == "Cubic")),
                  nrow(subset(aic_info1, best_fit_aic == "Briere")),
                  nrow(subset(aic_info1, best_fit_aic == "Schoolfield")),
                  nrow(subset(bic_info1, best_fit_aic == "Quadratic")), 
                  nrow(subset(bic_info1, best_fit_aic == "Cubic")),
                  nrow(subset(bic_info1, best_fit_aic == "Briere")),
                  nrow(subset(bic_info1, best_fit_aic == "Schoolfield")))

count_basic2 <- c(nrow(subset(aic_info2, best_fit_aic == "Quadratic")), 
                  nrow(subset(aic_info2, best_fit_aic == "Cubic")),
                  nrow(subset(aic_info2, best_fit_aic == "Briere")),
                  nrow(subset(aic_info2, best_fit_aic == "Schoolfield")),
                  nrow(subset(bic_info2, best_fit_aic == "Quadratic")), 
                  nrow(subset(bic_info2, best_fit_aic == "Cubic")),
                  nrow(subset(bic_info2, best_fit_aic == "Briere")),
                  nrow(subset(bic_info2, best_fit_aic == "Schoolfield")))

bestFitOccur <- data.frame(Model_Name=factor(rep(colnames(aic_info)[4:7], 2), level = colnames(aic_info)[4:7]),
                          Best_Fit_Count=count_all,
                          total=rep(nrow(aic_info), 8),
                          Method=(rep(c("AIC", "BIC"), each = 4)))

bestFitOccur12 <- data.frame(Model_Name=factor(rep(colnames(aic_info1)[4:7], 4), level = colnames(aic_info1)[4:7]),
                           Best_Fit_Count=c(count1, count2),
                           total=rep(c(nrow(aic_info1), nrow(aic_info2)), each = 8),
                           Method=(rep(c("AIC", "BIC"), each = 4, 2)),
                           Trits=rep(c("Photosynthesis", "Respiration"), each = 8))

bestFitall12 <- data.frame(Model_Name=factor(rep(colnames(aic_info1)[4:7], 8), level = colnames(aic_info1)[4:7]),
                             Best_Fit_Count=c(count1, count2, count_basic1, count_basic2),
                             total = rep(c(nrow(aic_info1), nrow(aic_info2)), each = 8, 2),
                             Assessment=(rep(c("AIC", "BIC"), each = 4, 4)),
                             Trits=rep(c("Photosynthesis", "Respiration"), each = 8, 2),
                             Method=rep(c("Rule of AIC/BIC difference", "Rule of single best fit"), each = 16))

bestFitall12$Percentage = round(bestFitall12$Best_Fit_Count/bestFitall12$total, 2)
bestFitall12$Label = paste(bestFitall12$Best_Fit_Count,"(", bestFitall12$Percentage * 100, "%", ")", sep = "")
bestFitOccur12$Percentage = round(bestFitOccur12$Best_Fit_Count/bestFitOccur12$total, 2)
bestFitOccur12$Label = paste(bestFitOccur12$Best_Fit_Count,"(", bestFitOccur12$Percentage * 100, "%", ")", sep = "")
bestFitOccur$Percentage = round(bestFitOccur$Best_Fit_Count/bestFitOccur$total, 2)
bestFitOccur$Label = paste(bestFitOccur$Best_Fit_Count,"(", bestFitOccur$Percentage * 100, "%", ")", sep = "")

print("The number of datasets in which cubic and Schoolfield both are the best-fit model by AIC difference rule:")
cat(nrow(subset(aic_info, cubic_best == 1 & sch_best == 1)), '\n')

print("The number of datasets in which cubic and Schoolfield both are the best-fit model by BIC difference rule:")
cat(nrow(subset(bic_info, cubic_best == 1 & sch_best == 1)), '\n')

aic_info <- cbind(aic_info[1:11], aic_info[15])
bic_info <- cbind(bic_info[1:11], bic_info[15])

write.csv(aic_info, "../Results/aic_info_all.csv", row.names = F)
print("The AIC values of all models have been saved into Results Directory! ")
write.csv(bic_info, "../Results/bic_info_all.csv", row.names = F)
print("The BIC values of all models have been saved into Results Directory! ")

# Saving the plot of best-fit occurence in to local directory
# ------------------------------------------------------------------------------------
pdf("../Results/bestfit_delta_all.pdf")
ggplot(data = bestFitOccur, aes(x = Model_Name, y = Best_Fit_Count, fill = Method)) +
  geom_histogram(stat = "identity", position = "dodge") +
  geom_text(aes(label=Label), position=position_dodge(1), vjust = -0.5, colour="black", size = 3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab("Model Names") +
  ylab("Counts") +
  theme(axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "bottom")
dev.off()

print("The figure of the best-fit model occurrence for the whole datasets has been saved into the results directory!")

pdf("../Results/bestfit_delta_grouped.pdf")
ggplot(data = bestFitall12, aes(x = Model_Name, y = Best_Fit_Count, fill = Assessment),) +
  geom_histogram(stat = "identity", position = "dodge") +
  facet_grid(Method ~ Trits) +
  geom_text(aes(label=Label), position=position_dodge(1), vjust = -0.5, colour="black", size = 2) +
  theme_classic() +
  xlab("Model Names") +
  ylab("Counts") +
  theme(axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "bottom")
dev.off()

print("The figure of the best-fit model occurrence for the grouped datasets has been saved into the results directory!")




