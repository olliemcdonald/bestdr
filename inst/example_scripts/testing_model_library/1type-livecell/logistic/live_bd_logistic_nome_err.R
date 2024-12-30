### Logistic Concentration of Birth and Death Rates
### w/ Additive Error
### No Mixed Effects

################################################################################
library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters:
# Generate curve from 4-parameter logistic function for birth and death
b0 <- 1 / 16 # left asymptote
bi <- 0.01 # right asymptote
b50 <- -2 # log-dose at midpoint of birthrate (should switch to log scale)
bh <- 2 # Hill Coefficient
db <- b0 - bi # Difference in asymptotes (for parameterization)

d0 <- 0.005
di <- 0.08
d50 <- 0
dh <- 1
dd <- di - d0

obs_err <- 0.02


param_truth <- data.frame(
  bi = bi, b50 = b50, bh = bh, db = db,
  d0 = d0, d50 = d50, dh = dh, dd = dd,
  obs_err = obs_err
)


concentrations <- c(0, 10^seq(-3, 2, 0.5))
log_conc <- log(concentrations)
b_curve <- function(lc, bi, db, bh, b50) bi + db / (1 + exp(bh * (lc - b50)))
d_curve <- function(lc, d0, dd, dh, d50) d0 + dd - dd / (1 + exp(dh * (lc - d50)))

b_truth <- b_curve(log_conc, bi, db, bh, b50)
d_truth <- d_curve(log_conc, d0, dd, dh, d50)
g_truth <- b_truth - d_truth
curve_truth <- data.frame(concentration = concentrations, log_conc = log_conc, b = b_truth, d = d_truth)
curve_truth %>% ggplot(aes(x = concentrations)) +
  geom_point(aes(y = b_truth), color = "blue") +
  geom_line(aes(y = b_truth), color = "blue") +
  geom_point(aes(y = d_truth), color = "red") +
  geom_line(aes(y = d_truth), color = "red") +
  geom_point(aes(y = g_truth), color = "magenta") +
  geom_line(aes(y = g_truth), color = "magenta") +
  scale_x_log10() +
  labs(x = "Concentration", y = "True Values") +
  theme_bw()


N <- 10 # Number of samples
tot_time <- 72
dt <- 2 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 2000 # Number of ancestors


# Simulate process with estipop
set.seed(8475)
dat <- data.frame()
for (i in seq_along(concentrations)) {
  params <- c(b_truth[i], d_truth[i])
  # time-homogeneous two-type model
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
    transition(rate = rate(params[2]), parent = 1, offspring = 0) # death
  ) # clearance

  dat_temp <- branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- concentrations[i]
  dat_temp$log_conc <- log_conc[i]
  dat_temp$b = b_truth[i]
  dat_temp$d = d_truth[i]
  dat <- bind_rows(dat, dat_temp)
}

dat <- dat %>%
  rename(sample = rep) %>%
  group_by(sample, concentration) %>%
  mutate(
    type1_prev = lag(type1),
    dt = signif(time - lag(time), 6),
    type1_mean = type1_prev * exp((b - d) * dt),
    error = rnorm(n(), mean = 0, sd = obs_err * type1_mean),
    count = ifelse(is.na(type1_prev), type1, round(type1 + error)),
    count_prev = lag(count)
)


# Plot to look at data
dat %>% ggplot(aes(x = time, group = as.factor(sample))) +
  geom_line(aes(y = count), col = "black") + # geom_line(aes(y = dead), col = "red") +
  facet_wrap(concentration ~ .) +
  theme_bw() +
  scale_y_log10()

dat <- dat %>% filter(!is.na(dt))


################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/logistic/birthdeath_error.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)


# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)

data.list <- list(
  N = nrow(dat),
  count = dat$count,
  count_prev = dat$count_prev,
  dt = dat$dt,
  c_idx = dat$log_conc_idx,
  nc = num_conc,
  conc = log_conc
)
prior.list <- list(
  prior_mu_bi = 0.05,
  prior_s_bi = 0.05,
  prior_mu_db = 0.1,
  prior_s_db = 0.05,
  prior_mu_b50 = 0,
  prior_s_b50 = 1,
  prior_mu_bh = 1,
  prior_s_bh = 0.5,
  prior_mu_d0 = 0.01,
  prior_s_d0 = 0.01,
  prior_mu_dd = 0.1,
  prior_s_dd = 0.05,
  prior_mu_d50 = 0,
  prior_s_d50 = 1,
  prior_mu_dh = 1,
  prior_s_dh = 1,
  prior_mu_obs_err = 0.05,
  prior_s_obs_err = 0.02
)
initfun <- function() {
  list(
    bi = runif(1, 0, 0.05),
    db = runif(1, 0, 0.1),
    b50 = runif(1, -2, 2),
    bh = runif(1, 0, 2),
    d0 = runif(1, 0, 0.05),
    dd = runif(1, 0, 0.1),
    d50 = runif(1, -2, 2),
    dh = runif(1, 0, 2),
    obs_err = runif(1, 0, 0.05),
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

# model_posterior <- stanmodel$pathfinder(
#   data = c(data.list, prior.list),
#   init = initfun,
#   refresh = 10
# )
model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  init = initfun,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=10)


print(model_posterior$summary(), n = 15)
params_posterior <- model_posterior$draws(variables = c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh", "obs_err"), format = "df")

################################################################################
# Compare to ground truth
param_truth <- param_truth %>% pivot_longer(everything(), names_to = "stat")

params_posterior_l <- params_posterior %>%
  select(.chain, bi, db, b50, bh, d0, dd, d50, dh, obs_err) %>%
  pivot_longer(!c(".chain"), names_to = "stat")

# Marginal Posteriors
params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = param_truth, aes(xintercept = value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()


# Samples of posterior curves
n_posteriors <- 300
crange <- c(10^seq(-6, 4, 0.1))
lcrange <- log(crange)
sample_post <- sample(nrow(params_posterior), n_posteriors, replace = FALSE)
b_post <- data.frame(
  crange,
  sapply(sample_post, function(x) {
    b_curve(
      lcrange,
      params_posterior$bi[x],
      params_posterior$db[x],
      params_posterior$bh[x],
      params_posterior$b50[x]
    )
  })
) %>%
  pivot_longer(cols = !crange, names_to = "sample", values_to = "b")
d_post <- data.frame(
  crange,
  sapply(sample_post, function(x) {
    d_curve(
      lcrange,
      params_posterior$d0[x],
      params_posterior$dd[x],
      params_posterior$dh[x],
      params_posterior$d50[x]
    )
  })
) %>%
  pivot_longer(cols = !crange, names_to = "sample", values_to = "d")
posterior_samples <- b_post %>%
  inner_join(d_post) %>%
  mutate(crange = ifelse(crange == 0, 0.000001, crange))
posterior_samples %>%
  ggplot(aes(x = crange)) +
  geom_line(aes(y = b, group = sample), color = "darkblue", alpha = 0.1) +
  geom_line(aes(y = d, group = sample), color = "darkred", alpha = 0.1) +
  geom_point(data = curve_truth, aes(x = concentration, y = b), color = "blue", size = 2) +
  geom_point(data = curve_truth, aes(x = concentration, y = d), color = "red", size = 2) +
  scale_x_log10() +
  theme_bw()




################################################################################
# Compare to single dose cases
singledose_model <- cmdstan_model("./stan/live_singledose/live_bd_singledose_nome_err.stan")

alldoses_draws <- data.frame()
for (i in unique(dat$log_conc_idx)) {
  dat_temp <- dat %>% filter(log_conc_idx == i)

  data.list <- list(
    N = nrow(dat_temp),
    count = dat_temp$count,
    count_prev = dat_temp$count_prev,
    dt = dat_temp$dt
  )
  prior.list <- list(
    prior_mu_b = 0.1,
    prior_mu_d = 0.1,
    prior_s_b = 0.05,
    prior_s_d = 0.05,
    prior_mu_obs_err = 0.05,
    prior_s_obs_err = 0.05
  )
  initfun <- function() {
    list(
      b = runif(1, 0, 0.2),
      d = runif(1, 0, 0.2),
      obs_err = runif(1, 0, 0.1) # ,
      # count_err = runif(data.list$N, -0.5, 0.5)
    )
  }
  
  singledose_posterior <- singledose_model$sample(
    data = c(data.list, prior.list),
    chains = 2,
    parallel_chains = 2,
    init = initfun,
    iter_warmup = 1000,
    iter_sampling = 500,
    refresh = 100
  )

  singledose_draws <- singledose_posterior$draws(variables = c("b", "d", "obs_err"), format = "df")
  singledose_draws$concentration <- concentrations[i]
  
  alldoses_draws <- rbind(alldoses_draws, singledose_draws)

}


posterior_samples %>%
  ggplot(aes(x = crange)) +
  geom_line(aes(y = b, group = sample), color = "darkblue", alpha = 0.1) +
  geom_line(aes(y = d, group = sample), color = "darkred", alpha = 0.1) +
  geom_violin(data = alldoses_draws, aes(x = concentration, y = b, group = concentration), color = "blue") +
  geom_violin(data = alldoses_draws, aes(x = concentration, y = d, group = concentration), color = "red") +
  geom_point(data = curve_truth, aes(x = concentration, y = b), color = "blue", size = 2) +
  geom_point(data = curve_truth, aes(x = concentration, y = d), color = "red", size = 2) +
  scale_x_log10() +
  theme_bw()


