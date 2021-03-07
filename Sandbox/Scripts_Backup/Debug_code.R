#####################################################################################################
############################## The linear model example of the given ID #############################
#####################################################################################################

require(minpack.lm)
require(tidyverse)
require(ggplot2)

d <- read.csv("../Data/ThermRespData.csv")

## Find subset
d_sub <- subset(d, ID == 2)
dim(d_sub)
ggplot(d_sub, aes(x = ConTemp, y = OriginalTraitValue))+
  geom_point(color = "blue")

## Linear Model Fit
lin_model <- lm(OriginalTraitValue ~ ConTemp, data = d_sub)
summ_lin <- summary(lin_model)
summ_lin$r.squared
AIC(lin_model)
BIC(lin_model)
## Model Prediction
Lengths <- seq(min(d_sub$ConTemp), max(d$ConTemp), len = 200)
lin_pred <- predict.lm(lin_model, data.frame(ConTemp = Lengths))

## Model Plot
plot(d_sub$ConTemp, d_sub$OriginalTraitValue)
lines(Lengths, lin_pred, col = "blue", lwd = 2)

## Model ggplot
ggplot(data = d_sub, aes(x = ConTemp, y = OriginalTraitValue)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_text(aes(x = median(ConTemp), y = max(OriginalTraitValue) + 5,
                label = paste("R**2: ", signif(summ_lin$r.squared, 3))),
            size = 4, parse = T,
            colour = "red")

ggplot(data = d_sub, aes(x = ConTemp, y = OriginalTraitValue)) +
  geom_point() +
  geom_abline(
    intercept = summ_lin$coefficients[1][1],
    slope = summ_lin$coefficients[2][1],
    colour = "blue") +
  geom_text(aes(x = median(ConTemp), y = max(OriginalTraitValue),
                label = paste("R**2: ", signif(summ_lin$r.squared, 3))),
            size = 4, parse = T,
            colour = "red")

#####################################################################################################
############################## The linear model example of the given ID #############################
#####################################################################################################

#---------------------------------------------------------------------------------------------------#

#####################################################################################################
############################## The non-linear model example of the given ID #########################
#####################################################################################################

## Find subset
rm(list = ls())
d <- read.csv("MiniProject/Data/ThermRespData.csv")
d_sub <- subset(d, ID == 1)
dim(d_sub)
ggplot(d_sub, aes(x = ConTemp, y = OriginalTraitValue))+
  geom_point(color = "blue")

quad_func <- function(B0, B1, B2, x){
  return(B0 + B1 * x + B2 * x^2)
}

B0_start <- .1
B1_start <- .1
B2_start <- .1

lin_model <- lm(OriginalTraitValue ~ ConTemp, data = d_sub)
quad_model <- nlsLM(OriginalTraitValue ~ quad_func(B0, B1, B2, ConTemp), data = d_sub, 
                    start = list(B0 = B0_start, B1 = B1_start, B2 = B2_start))

summ_lin <- summary(lin_model)
summ_quad <- summary(quad_model)

x_values <- seq(min(d_sub$ConTemp), max(d_sub$ConTemp), len = 200)

lin_pred <- predict.lm(lin_model, data.frame(ConTemp = x_values))
quad_pred <- quad_func(
  x = x_values,
  B0 = coef(quad_model)["B0"],
  B1 = coef(quad_model)["B1"],
  B2 = coef(quad_model)["B2"]
)

df_lin <- data.frame(x_values, lin_pred)
df_lin$model <- "Linear Regression"
names(df_lin) <- c("ConTemp", "OriginalTraitValue", "Model")

df_quad <- data.frame(x_values, quad_pred)
df_quad$model <- "Quadratic Regression"
names(df_quad) <- c("ConTemp", "OriginalTraitValue", "Model")

model_frame <- rbind(df_lin, df_quad)

ggplot(d_sub, aes(x = ConTemp, y = OriginalTraitValue)) +
  geom_point(size = 3) + 
  geom_line(data = model_frame, aes(x = ConTemp, y = OriginalTraitValue, col = Model), size = 1) +
  theme() +
  labs(x = "ConTemp", y = "OriginalTraitValue")


#####################################################################################################
############################## The non-linear model example of the given ID #########################
#####################################################################################################