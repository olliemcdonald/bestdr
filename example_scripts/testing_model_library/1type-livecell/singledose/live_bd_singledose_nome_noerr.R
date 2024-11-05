### Single Dose Estimation of Birth and Death Rates
### No Error
### No Mixed Effects

################################################################################
# Set Project Directory Location---------------------------------
library(tidyverse)
library(estipop)
library(cmdstanr)

################################################################################
# Generate data -------------------------------------------------
# Generate Parameters:
br <- 1 # birth rate
dr <- 0.5 # death rate
lambda <- br - dr
N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors

# Simulate process with estipop
set.seed(89273459734)
params <- c(br, dr)
model <- process_model(
  transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
  transition(rate = rate(params[2]), parent = 1, offspring = 0)
) # death
dat <- branch(model, params, init_pop, times, N, silent = T)

dat <- dat %>%
  rename(sample = rep, live = type1) %>%
  group_by(sample) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    live_prev = lag(live)
  )

# dat %>% ggplot(aes(x = time, y = live, group = sample)) + geom_line()
dat <- dat %>% filter(!is.na(dt))

################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/singledose/birthdeath.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)


data.list <- list(
  N = nrow(dat),
  count = dat$live,
  count_prev = dat$live_prev,
  dt = dat$dt
)
prior.list <- list(
  prior_mu_b = 1,
  prior_mu_d = 0.5,
  prior_s_b = 1,
  prior_s_d = 1
)
initfun <- function() {
  list(
    b = runif(1, 0, 1),
    d = runif(1, 0, 1),
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 4,
  parallel_chains=4,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000)
)

params_posterior <- model_posterior$draws(variables = c("b", "d"), format = "df")

################################################################################
# Compare to ground truth
ground_truth <- data.frame(stat = c("b", "d"), true_value = c(br, dr))

params_posterior_l <- params_posterior %>% select(.chain, b, d) %>%
  pivot_longer(!c(".chain"), names_to="stat")

params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = ground_truth, aes(xintercept = true_value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()


params_posterior_l %>%
  left_join(ground_truth) %>%
  mutate(diff = value - true_value) %>%
  group_by(.chain, stat) %>%
  summarize(mean_diff = mean(diff)) %>%
  pivot_wider(names_from = stat, values_from = mean_diff)
