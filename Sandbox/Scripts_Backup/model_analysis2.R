## Things to be done
# - Read tables of fitting info from results directory and transform it as DataFrame
# - Merge ID, aic, bic, aicc, rsq columns as one table
# - Add best_fit column to find the model with the lowest aicc value
# - Calculate percentage of each model with lowest aicc value during fitting all IDs
# - Plot bar chart of the percentage
library(ggplot2)
library(dplyr)

rm(list = ls())

# Loading dataset
quad_info <- read.csv("../Results/Fitting_info/quad_info.csv")
cubic_info <- read.csv("../Results/Fitting_info/cubic_info.csv")
# briere_info <- read.csv("../Results/Fitting_info/Briere_info.csv")
school_info <- read.csv("../Results/Fitting_info/Schoolfield_info.csv")

quad_info <- quad_info %>% select(id, 7:10) %>%
  rename(quad.rsq = Rsq, quad.aic = aic, quad.bic = bic, quad.aicc = aicc)
cubic_info <- cubic_info %>% select(8:11) %>%
  rename(cubic.rsq = Rsq, cubic.aic = aic, cubic.bic = bic, cubic.aicc = aicc)
# briere_info <- briere_info %>% select(8:11) %>%
#   rename(briere.rsq = Rsq, briere.aic = aic, briere.bic = bic, briere.aicc = aicc)
school_info <- school_info %>% select(9:12) %>%
  rename(school.rsq = Rsq, school.aic = aic, school.bic = bic, school.aicc = aicc)

fitting_info <- cbind(quad_info, cubic_info, school_info)

id = fitting_info$id
aicc_info <- fitting_info %>% select(quad.aicc, cubic.aicc, school.aicc)
row.names(aicc_info) <- id
aicc_info$best_fit_aicc <- colnames(aicc_info)[apply(aicc_info, 1, which.min)]

aic_info <- fitting_info %>% select(quad.aic, cubic.aic,school.aic)
row.names(aic_info) <- id
aic_info$best_fit_aic <- colnames(aic_info)[apply(aic_info, 1, which.min)]

bic_info <- fitting_info %>% select(quad.bic, cubic.bic,  school.bic)
row.names(bic_info) <- id
bic_info$best_fit_bic <- colnames(bic_info)[apply(bic_info, 1, which.min)]

model_evaluation <- cbind(id, aicc_info, aic_info, bic_info)
write.csv(model_evaluation, "../Results/Fitting_info/model_evaluation.csv", row.names = FALSE)

# Save all plots
###########################################################################
pdf("../Results/Fitting_info/best_fit_aicc.pdf")
ggplot(data = aicc_info) +
  geom_histogram(aes(x = best_fit_aicc), stat = "count") +
  labs(title = "Best Fit of Each Model Besed on AICc Values") +
  xlab("Model Names") +
  ylab("Counts")
dev.off()

pdf("../Results/Fitting_info/best_fit_aic.pdf")
ggplot(data = aic_info) +
  geom_histogram(aes(x = best_fit_aic), stat = "count") +
  labs(title = "Best Fit of Each Model Besed on AIC Values") +
  xlab("Model Names") +
  ylab("Counts")
dev.off()

pdf("../Results/Fitting_info/best_fit_bic.pdf")
ggplot(data = bic_info) +
  geom_histogram(aes(x = best_fit_bic), stat = "count") +
  labs(title = "Best Fit of Each Model Besed on BIC Values") +
  xlab("Model Names") +
  ylab("Counts")
dev.off()
###########################################################################

# Calculating percentage of beconming best fit for each model
# based on AICc
#################################################################################################
# perc_briere_aicc <- nrow(aicc_info %>% subset(best_fit_aicc == "briere.aicc")) / nrow(aicc_info)
perc_school_aicc <- nrow(aicc_info %>% subset(best_fit_aicc == "school.aicc")) / nrow(aicc_info)
perc_cubic_aicc <- nrow(aicc_info %>% subset(best_fit_aicc == "cubic.aicc")) / nrow(aicc_info)
perc_quad_aicc <- nrow(aicc_info %>% subset(best_fit_aicc == "quad.aicc")) / nrow(aicc_info)
cat("Percentage of Best Fit for Each Model basd on AICc Values", "\n")
cat("--------------------------------------------------------------", "\n")
cat("Quadratic Model: ", round(perc_quad_aicc * 100, 3), "%", "\n")
cat("Cubic Model: ", round(perc_cubic_aicc * 100, 3), "%", "\n")
# cat("Briere Model: ", round(perc_briere_aicc * 100, 3), "%", "\n")
cat("Simplified Schoolfield Model: ", round(perc_school_aicc * 100, 3), "%", "\n")
cat("--------------------------------------------------------------", "\n")

# based on AIC
# perc_briere_aic <- nrow(aic_info %>% subset(best_fit_aic == "briere.aic")) / nrow(aic_info)
perc_school_aic <- nrow(aic_info %>% subset(best_fit_aic == "school.aic")) / nrow(aic_info)
perc_cubic_aic <- nrow(aic_info %>% subset(best_fit_aic == "cubic.aic")) / nrow(aic_info)
perc_quad_aic <- nrow(aic_info %>% subset(best_fit_aic == "quad.aic")) / nrow(aic_info)
cat("\n", "Percentage of Best Fit for Each Model basd on AIC Values", "\n")
cat("--------------------------------------------------------------", "\n")
cat("Quadratic Model: ", round(perc_quad_aic * 100, 3), "%", "\n")
cat("Cubic Model: ", round(perc_cubic_aic * 100, 3), "%", "\n")
# cat("Briere Model: ", round(perc_briere_aic * 100, 3), "%", "\n")
cat("Simplified Schoolfield Model: ", round(perc_school_aic * 100, 3), "%", "\n")
cat("--------------------------------------------------------------", "\n")

# based on BIC
# perc_briere_bic <- nrow(bic_info %>% subset(best_fit_bic == "briere.bic")) / nrow(bic_info)
perc_school_bic <- nrow(bic_info %>% subset(best_fit_bic == "school.bic")) / nrow(bic_info)
perc_cubic_bic <- nrow(bic_info %>% subset(best_fit_bic == "cubic.bic")) / nrow(bic_info)
perc_quad_bic <- nrow(bic_info %>% subset(best_fit_bic == "quad.bic")) / nrow(bic_info)
cat("\n", "Percentage of Best Fit for Each Model basd on BIC Values", "\n")
cat("--------------------------------------------------------------", "\n")
cat("Quadratic Model: ", round(perc_quad_bic * 100, 3), "%", "\n")
cat("Cubic Model: ", round(perc_cubic_bic * 100, 3), "%", "\n")
# cat("Briere Model: ", round(perc_briere_bic * 100, 3), "%", "\n")
cat("Simplified Schoolfield Model: ", round(perc_school_bic * 100, 3), "%", "\n")
cat("--------------------------------------------------------------", "\n")
#################################################################################################
