# Loading packages
require(minpack.lm)
require(tidyverse)
require(ggplot2)

# Clear all the variables
rm(list = ls())

# Loading data
d <- read.csv("../Data/ThermRespData.csv")
dim(d)
length(colnames(d))
# 78 Columns

# Delete null values in OriginalTraitValue and ConTemp
d <- subset(d, !is.null(OriginalTraitValue))
d <- subset(d, !is.null(ConTemp))
dim(d)
dplyr::glimpse(d)

# See the value of these three variables
unique(d$OriginalTraitUnit)
## units of the response variable
unique(d$ConTempUnit)
## Units of the independent variable
unique(d$ID)

# To see the range of variables
range(d$OriginalTraitValue)
range(scale(d$OriginalTraitValue))
range(d$ConTemp)
range(scale(d$ConTemp))

# Linear Model Fitting for different IDs

for (i in 1:6){
  d_sub <- subset(d, ID == i)
  
  ## Define functions
  briere_func <- function(T, B0, T0, Tm) {
    return(B0 * T * (T - T0) * (abs(Tm-T)^(1/2)) * as.numeric(T < Tm) * as.numeric(T > T0))
  }
  
  ## Start points of parameters in nonlinear functions
  B0_start <- 0.01
  T0_start <- 10
  Tm_start <- 40
  
  ## Modeling
  lin_model <- lm(OriginalTraitValue ~ ConTemp, data = d_sub)
  quad_model <- lm(OriginalTraitValue ~ poly(ConTemp, 2), data = d_sub)
  poly_model <- lm(OriginalTraitValue ~ poly(ConTemp, 3), data = d_sub)
  briere_model <- nlsLM(OriginalTraitValue ~ briere_func(T = ConTemp, B0, T0, Tm), data = d_sub,
                        start = list(B0 = B0_start, T0 = T0_start, Tm = Tm_start))
  
  ## Summary of models
  summ_lin <- summary(lin_model)
  summ_quad <- summary(quad_model)
  summ_poly <- summary(poly_model)
  summ_briere <- summary(briere_model)
  
  ## create x values for model prediction
  x_values <- seq(min(d_sub$ConTemp), max(d_sub$ConTemp), len = 300)
  
  ## Model prediction
  lin_pred <- predict.lm(lin_model, data.frame(ConTemp = x_values))
  quad_pred <- predict.lm(quad_model, data.frame(ConTemp = x_values))
  poly_pred <- predict.lm(poly_model, data.frame(ConTemp = x_values))
  
  briere_pred <- briere_func(
    T = x_values,
    B0 = coef(briere_model)["B0"],
    T0 = coef(briere_model)["T0"],
    Tm = coef(briere_model)["Tm"]
  )
  
  # poly_pred <- poly_func(
  #   x = x_values,
  #   B0 = coef(poly_model)["B0"],
  #   B1 = coef(poly_model)["B1"],
  #   B2 = coef(poly_model)["B2"], 
  #   B3 = coef(poly_model)["B3"]
  # )
  
  ## Create model frame for plotting
  df_lin <- data.frame(x_values, lin_pred)
  df_lin$model <- "Linear Regression"
  names(df_lin) <- c("ConTemp", "OriginalTraitValue", "Model")
  
  df_quad <- data.frame(x_values, quad_pred)
  df_quad$model <- "Quadratic Regression"
  names(df_quad) <- c("ConTemp", "OriginalTraitValue", "Model")
  
  df_poly <- data.frame(x_values, poly_pred)
  df_poly$model <- "Polynomial Regression"
  names(df_poly) <- c("ConTemp", "OriginalTraitValue", "Model")
  
  df_briere <- data.frame(x_values, briere_pred)
  df_briere$model <- "Briere Model"
  names(df_briere) <- c("ConTemp", "OriginalTraitValue", "Model")
  
  model_frame <- rbind(df_lin, df_quad, df_poly, df_briere)
  
  ## plotting models
  pdfFilePath <- paste("../Results/model_fitting", i, ".pdf", sep = "")
  pdf(pdfFilePath)
  p <- ggplot(d_sub, aes(x = ConTemp, y = OriginalTraitValue)) +
    geom_point() + 
    geom_line(data = model_frame, aes(x = ConTemp, y = OriginalTraitValue, col = Model), size = 1)
    # geom_text(aes(x = mean(ConTemp), y = max(OriginalTraitValue),
    #               label = paste("R**2: ", signif(summ_lin$r.squared, 3))),
    #           size = 4, parse = T,
    #           colour = "red")
  print(p)
  dev.off()
}


