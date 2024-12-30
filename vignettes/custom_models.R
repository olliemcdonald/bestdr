## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bestdr)
library(tidyverse)

## ----birthdeath-noDR----------------------------------------------------------
b <- 1
d <- 0.1

N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(100) # Number of ancestors

set.seed(90210)
params <- c(b, d)
model <- estipop::process_model(
  estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = c(2)), # birth
  estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = c(1)) # death
)
dat <- estipop::branch(model, params, init_pop, times, N, silent = T)

# Plot the longitudinal cell growth data
dat %>% ggplot(aes(x = time, y = type1, group = rep)) +
  geom_line() +
  scale_y_log10()

