library(ggplot2)
library(dplyr)
library(tidyverse)

rm(list = ls())

# Loading dataset
quad_info <- read.csv("../Sandbox/Results Backup/Fitting_info_best/quad_info.csv")
cubic_info <- read.csv("../Sandbox/Results Backup/Fitting_info_best/cubic_info.csv")
briere_info <- read.csv("../Sandbox/Results Backup/Fitting_info_best/Briere_info.csv")
school_info <- read.csv("../Sandbox/Results Backup/Fitting_info_best/Schoolfield_info.csv")

# Delete NAs
school_info <- school_info %>% subset(is.na(aic) == F)

# Choose columns
quad_info <- quad_info %>% select(1:4) %>%
  rename(quad.a = param1, quad.b1 = param2, quad.b2 = param3)

ggplot(quad_info, aes(x = "Intercept", y = quad.a)) +
  geom_boxplot()

cubic_info <- cubic_info %>% select(1:5) %>%
  rename(cubic.a = param1, cubic.b1 = param2, cubic.b2 = param3, cubic.b3 = param4)

briere_info <- briere_info %>% select(1:4) %>%
  rename(bri.B0 = B0, bri.T0 = T0, bri.Tm = Tm)

school_info <- school_info %>% select(1:5) %>%
  rename(sch.B0 = B0, sch.E = E, sch.Th = Th, sch.Eh = Eh)

# Plot figures
ggplot2(quad_info, aes(x = )) +
  geom_boxplot(data = quad_info, ()