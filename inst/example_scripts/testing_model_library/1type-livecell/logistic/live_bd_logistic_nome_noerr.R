### Logistic Concentration of Birth and Death Rates
### No Error
### No Mixed Effects

################################################################################
library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters:
# Generate curve from 4-parameter logistic function for birth and death
b0 <- 1 # left asymptote
bi <- 0.2 # right asymptote
b50 <- -3 # log-dose at midpoint of birthrate (should switch to log scale)
bh <- 0.5 # Hill Coefficient
db <- b0 - bi # Difference in asymptotes (for parameterization)

d0 <- 0.01
di <- 0.5
d50 <- 0
dh <- 1
dd <- di - d0

param_truth <- data.frame(
  bi = bi, b50 = b50, bh = bh, db = db,
  d0 = d0, d50 = d50, dh = dh, dd = dd
)


concentrations <- c(0, 10^seq(-4, 4, 1))
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
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors


# Simulate process with estipop
set.seed(45387687634)
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
  dat <- bind_rows(dat, dat_temp)
}

dat <- dat %>%
  rename(sample = rep, live = type1) %>%
  group_by(sample, concentration) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    live_prev = lag(live)
  )

# Plot to look at data
dat %>% ggplot(aes(x = time, group = as.factor(sample))) +
  geom_line(aes(y = live), col = "black") + # geom_line(aes(y = dead), col = "red") +
  facet_wrap(concentration ~ .) +
  theme_bw()
dat <- dat %>% filter(!is.na(dt))


################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/logistic/birthdeath.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)

data.list <- list(
  N = nrow(dat),
  count = dat$live,
  count_prev = dat$live_prev,
  dt = dat$dt,
  c_idx = dat$log_conc_idx,
  nc = num_conc,
  conc = log_conc
)
prior.list <- list(
  prior_mu_bi = 1,
  prior_s_bi = 0.3,
  prior_mu_db = 0.5,
  prior_s_db = 0.5,
  prior_mu_b50 = 0,
  prior_s_b50 = 1,
  prior_mu_bh = 1,
  prior_s_bh = 0.5,
  prior_mu_d0 = 0.3,
  prior_s_d0 = 0.3,
  prior_mu_dd = 0.5,
  prior_s_dd = 0.5,
  prior_mu_d50 = 0,
  prior_s_d50 = 1,
  prior_mu_dh = 1,
  prior_s_dh = 0.5
)
initfun <- function() {
  list(
    bi = runif(1, 0, 2),
    db = runif(1, 0, 2),
    b50 = runif(1, -2, 2),
    bh = runif(1, 0, 2),
    d0 = runif(1, 0, 2),
    dd = runif(1, 0, 2),
    d50 = runif(1, -2, 2),
    dh = runif(1, 0, 2),
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 4,
  parallel_chains = 4,
  refresh = 10,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)


params_posterior <- model_posterior$draws(variables = c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh"), format = "df")

################################################################################
# Compare to ground truth
param_truth <- param_truth %>% pivot_longer(everything(), names_to = "stat")

params_posterior_l <- params_posterior %>%
  select(.chain, bi, db, b50, bh, d0, dd, d50, dh) %>%
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

