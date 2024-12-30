### Single Dose Estimation of Birth, Death, Clearance Rates
### No Error
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
lambda <- br - dr
N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(200, 0) # Number of ancestors

# Simulate process with estipop
set.seed(102938)
params <- c(br, dr, kr)
# time-homogeneous two-type model
model <- process_model(
  transition(rate = rate(params[1]), parent = 1, offspring = c(2, 0)), # birth
  transition(rate = rate(params[2]), parent = 1, offspring = c(0, 1)), # death
  transition(rate = rate(params[3]), parent = 2, offspring = c(0, 0)) # clearance
) 
dat <- branch(model, params, init_pop, times, N, silent = T)

dat <- dat %>%
  rename(sample = rep, live = type1, dead = type2) %>%
  group_by(sample) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    live_prev = lag(live),
    dead_prev = lag(dead)
  )

# dat %>% ggplot(aes(x = time, y = live, group = sample)) + geom_line()
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
  count = cbind(dat$live, dat$dead),
  count_prev = cbind(dat$live_prev, dat$dead_prev),
  times_idx = rep(1, nrow(dat)),
  dt = dt
)
model.list <- list(
  e_mat = matrix(c(2, 0, 0, 1, 0, 0), nrow = 3, ncol = 2, byrow = T),
  p_vec = c(1, 1, 2)
)
prior.list <- list(
  prior_mu_b = 1,
  prior_mu_d = 0.5,
  prior_mu_k = 0.1,
  prior_s_b = 1,
  prior_s_d = 1,
  prior_s_k = 0.05
)
initfun <- function() {
  list(
    b = runif(1, 0.5, 1),
    d = runif(1, 0, 0.5),
    k = runif(1, 0., 0.1),
    count_err = matrix(runif(data.list$N * ntypes, -0.5, 0.5), ncol = ntypes)
  )
}

################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/2type-livedead/singledose/birthdeathclearance.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)


model_posterior <- stanmodel$sample(
  data = c(data.list, model.list, prior.list),
  chains = 1,
  parallel_chains=1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000)
)

params_posterior <- model_posterior$draws(variables = c("b", "d", "k"), format = "df")

################################################################################
# Compare to ground truth
ground_truth <- data.frame(stat = c("b", "d", "k"), true_value = c(br, dr, kr))

params_posterior_l <- params_posterior %>% select(.chain, b, d, k) %>%
  pivot_longer(!c(".chain"), names_to="stat")

params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = ground_truth, aes(xintercept = true_value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()


################################################################################
# Exponential Matrix Version Solver - not really working
# stanmodel_quad <- cmdstan_model("./stan/livedead_singledose/livedead_bdk_singledose_nome_noerr_quad.stan", force_recompile=T)


# zz <- file("~/Desktop/output.txt", open = "wt")
# sink(zz, type = "output")
# sink(zz, type = "message")
# model_posterior_quad <- stanmodel_quad$sample(
#   data = c(data.list, model.list, prior.list),
#   chains = 1,
#   parallel_chains=1,
#   refresh = 1,
#   init = initfun,
#   iter_warmup = 1,
#   iter_sampling = 1)
# sink()
# sink(type = "message")
# close(zz)
