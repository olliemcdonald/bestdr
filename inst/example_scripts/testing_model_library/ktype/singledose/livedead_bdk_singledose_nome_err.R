### Single Dose Estimation of Birth, Death, Clearance Rates
### 
### Error
### No Mixed Effects

################################################################################
library(tidyverse)
library(estipop)
library(cmdstanr)

################################################################################
# Generate data -------------------------------------------------
# Generate Parameters:
br <- 1 # birth rate
dr <- 0.5 # death rate
kr <- 0.02
obs_err <- 0.05


lambda <- br - dr
N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(200, 0) # Number of ancestors

# Simulate process with estipop
set.seed(032224)
params <- c(br, dr, kr)
# time-homogeneous two-type model
model <- process_model(
  transition(rate = rate(params[1]), parent = 1, offspring = c(2, 0)), # birth
  transition(rate = rate(params[2]), parent = 1, offspring = c(0, 1)), # death
  transition(rate = rate(params[3]), parent = 2, offspring = c(0, 0)) # clearance
) 
dat <- branch(model, params, init_pop, times, N, silent = T)

# Plot to look at data
dat %>%
  pivot_longer(!c(rep, time), names_to = "type") %>%
  ggplot(aes(x = time, group = interaction(type, as.factor(rep)), color = type)) +
  geom_line(aes(y = value)) +
  theme_bw() +
  scale_y_log10()

# Add errors
dat <- dat %>%
  mutate(
    type1_noerror = type1,
    type2_noerror = type2,
    er1 = round(type1 * rnorm(nrow(dat), 0, obs_err)),
    er2 = round(type2 * rnorm(nrow(dat), 0, obs_err)),
    type1 = type1 + er1,
    type2 = type2 + er2
  )

dat %>%
  select(c(rep, time, type1, type2)) %>%
  pivot_longer(!c(rep, time), names_to = "type") %>%
  ggplot(aes(x = time, group = interaction(type, as.factor(rep)), color = type)) +
  geom_line(aes(y = value)) +
  theme_bw() +
  scale_y_log10()


dat <- dat %>%
  rename(sample = rep) %>%
  group_by(sample) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    type1_prev = lag(type1),
    type2_prev = lag(type2)
  )


dat <- dat %>% filter(!is.na(dt))

################################################################################
# Set up inputs for sampler

dt <- array(unique(dat$dt))
ntypes <- 2
nevents <- 3

data.list <- list(
  ntypes = ntypes,
  nevents = nevents,
  n_dt_unique = length(unique(dat$dt)),
  N = nrow(dat),
  count = cbind(dat$type1, dat$type2),
  count_prev = cbind(dat$type1_prev, dat$type2_prev),
  times_idx = rep(1, nrow(dat)),
  dt = dt
)
model.list <- list(
  e_mat = matrix(c(2, 0, 0, 1, 0, 0), nrow = 3, ncol = 2, byrow = T),
  p_vec = c(1, 1, 2)
)
prior.list <- list(
  prior_mu = c(1, 0.5, 0.1),
  prior_s = c(1, 1, 0.05),
  prior_s_obs_err = 0.1
)
initfun <- function() {
  list(
    theta = c(runif(1, 0.5, 1), runif(1, 0, 0.5), runif(1, 0., 0.1)),
    count_err = matrix(runif(data.list$N * ntypes, -0.5, 0.5), ncol = ntypes)
  )
}

################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/ktype/singledose/ktype_error.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

model_posterior <- stanmodel$sample(
  data = c(data.list, model.list, prior.list),
  chains = 1,
  parallel_chains=1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)

params_posterior <- model_posterior$draws(variables = c("theta", "obs_err"), format = "df")
names(params_posterior)[1:3] <- c("b", "d", "k")
################################################################################
# Compare to ground truth
ground_truth <- data.frame(stat = c("b", "d", "k", "obs_err"), true_value = c(br, dr, kr, obs_err))

params_posterior_l <- params_posterior %>% select(.chain, b, d, k, obs_err) %>%
  pivot_longer(!c(".chain"), names_to="stat")

params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = ground_truth, aes(xintercept = true_value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()
