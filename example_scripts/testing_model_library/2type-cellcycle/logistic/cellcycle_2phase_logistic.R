### Logistic Concentration of Birth, Death, Clearance Rates 
### No Error
### No Mixed Effects

################################################################################

library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters:
# Model:
# G1 -> S/G2/M   (u_g1)
# G1 -> X        (d_g1)
# S/G2/M -> 2 G1 (u_g2)
# S/G2/M -> X    (d_g2)
# Generate Parameters:
# Rates are given as color transitions or color deaths

# Constants
u_g1 <- 1 / 8       
d_g1 <- 0.01
# Generate curve from 4-parameter logistic function for SG2M -> 2G1 and SG2M death
u_g2_0 <- 1 / 16
u_g2_i <- 1 / 100
u_g2_50 <- -3
u_g2_h <- 0.5
u_g2_delta <- u_g2_0 - u_g2_i

d_g2_0 <- 0.01  
d_g2_i <- 0.03
d_g2_50 <- -1 
d_g2_h <- 1
d_g2_delta <- d_g2_i - d_g2_0

param_truth <- data.frame(
  u_g1 = u_g1, d_g1 = d_g1,
  u_g2_i = u_g2_i, u_g2_delta = u_g2_delta, u_g2_50 = u_g2_50, u_g2_delta = u_g2_delta,
  d_g2_0 = d_g2_0, d_g2_delta = d_g2_delta, d_g2_50 = d_g2_50, d_g2_delta = d_g2_delta
)


concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)
b_curve <- function(lc, bi, db, bh, b50) bi + db / (1 + exp(bh * (lc - b50)))
d_curve <- function(lc, d0, dd, dh, d50) d0 + dd - dd / (1 + exp(dh * (lc - d50)))

u_g2_truth <- b_curve(log_conc, u_g2_i, u_g2_delta, u_g2_h, u_g2_50)
d_g2_truth <- d_curve(log_conc, d_g2_0, d_g2_delta, d_g2_h, d_g2_50)
curve_truth <- data.frame(concentration = concentrations, log_conc = log_conc, u_g2 = u_g2_truth, d_g2 = d_g2_truth)
curve_truth %>% ggplot(aes(x = concentrations)) +
  geom_point(aes(y = u_g2_truth), color = "blue") +
  geom_line(aes(y = u_g2_truth), color = "blue") +
  geom_point(aes(y = d_g2_truth), color = "red") +
  geom_line(aes(y = d_g2_truth), color = "red") +
  geom_hline(aes(yintercept = param_truth$u_g1), color = "darkgreen") +
  geom_hline(aes(yintercept = param_truth$d_g1), color = "darkorange") +
  scale_x_log10() +
  labs(x = "Concentration", y = "True Values") +
  theme_bw()


N <- 10 # Number of samples
tot_time <- 72
dt <- 3 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(1000, 0) # Number of ancestors


# Simulate process with estipop
set.seed(90210)
dat <- data.frame()
for (i in seq_along(concentrations)) {
  params <- c(u_g1, d_g1, u_g2_truth[i], d_g2_truth[i])
  # time-homogeneous two-type model
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = c(0, 1)),
    transition(rate = rate(params[2]), parent = 1, offspring = c(0, 0)),
    transition(rate = rate(params[3]), parent = 2, offspring = c(2, 0)),
    transition(rate = rate(params[4]), parent = 2, offspring = c(0, 0))
  )

  dat_temp <- branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- concentrations[i]
  dat_temp$log_conc <- log_conc[i]
  dat <- bind_rows(dat, dat_temp)
}

# Plot to look at data
dat %>%
  pivot_longer(!c(rep, time, concentration, log_conc), names_to = "type") %>%
  ggplot(aes(x = time, group = interaction(type, as.factor(rep)), color = type)) +
  geom_line(aes(y = value)) +
  facet_wrap(concentration ~ .) +
  theme_bw() +
  scale_y_log10()

dat <- dat %>%
  rename(sample = rep) %>%
  group_by(sample, concentration) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    type1_prev = lag(type1),
    type2_prev = lag(type2)
  )


dat <- dat %>% filter(!is.na(dt))


################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/2type-cellcycle/logistic/cellcycle-2phase.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
ntypes <- 2
nevents <- 4
dt <- array(unique(dat$dt))

# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)

data.list <- list(
  N = nrow(dat),
  ntypes = ntypes,
  nevents = nevents,
  n_dt_unique = length(unique(dat$dt)),
  count = cbind(dat$type1, dat$type2),
  count_prev = cbind(dat$type1_prev, dat$type2_prev),
  times_idx = rep(1, nrow(dat)), # when different times need to fix this
  dt = dt,
  c_idx = dat$log_conc_idx,
  nc = num_conc,
  conc = log_conc
)
model.list <- list(
  e_mat = matrix(c(0, 1, 0, 0, 2, 0, 0, 0), nrow = 4, ncol = 2, byrow = T),
  p_vec = c(1, 1, 2, 2)
)
prior.list <- list(
  prior_mu_u_g2_i = 1,
  prior_s_u_g2_i = 0.3,
  prior_mu_u_g2_delta = 1,
  prior_s_u_g2_delta = 0.3,
  prior_mu_u_g2_50 = 0,
  prior_s_u_g2_50 = 1,
  prior_mu_u_g2_h = 1,
  prior_s_u_g2_h = 0.5,
  prior_mu_d_g2_0 = 0.3,
  prior_s_d_g2_0 = 0.3,
  prior_mu_d_g2_delta = 0.5,
  prior_s_d_g2_delta = 0.5,
  prior_mu_d_g2_50 = 0,
  prior_s_d_g2_50 = 1,
  prior_mu_d_g2_h = 1,
  prior_s_d_g2_h = 0.5,
  prior_mu_u_g1 = 0.05,
  prior_s_u_g1 = 0.02,
  prior_mu_d_g1 = 0.05,
  prior_s_d_g1 = 0.02
)
initfun <- function() {
  list(
    u_g2_i = runif(1, 0, 2),
    u_g2_delta = runif(1, 0, 2),
    u_g2_50 = runif(1, -2, 2),
    u_g2_h = runif(1, 0, 2),
    d_g2_0 = runif(1, 0, 2),
    d_g2_delta = runif(1, 0, 2),
    d_g2_50 = runif(1, -2, 2),
    d_g2_h = runif(1, 0, 2),
    u_g1 = runif(1, 0, 0.1),
    d_g1 = runif(1, 0, 0.1),
    count_err = matrix(runif(data.list$N * ntypes, -0.5, 0.5), ncol = ntypes)
  )
}


model_posterior <- stanmodel$sample(
  data = c(data.list, model.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)


params_posterior <-
  model_posterior$draws(variables = c("u_g2_i", "u_g2_delta", "u_g2_50", "u_g2_h",
                                      "d_g2_0", "d_g2_delta", "d_g2_50", "d_g2_h",
                                      "u_g1", "d_g1"), format = "df")

################################################################################
# Compare to ground truth
param_truth <- param_truth %>% pivot_longer(everything(), names_to = "stat")

params_posterior_l <- params_posterior %>%
  select(.chain, "u_g2_i", "u_g2_delta", "u_g2_50", "u_g2_h",
         "d_g2_0", "d_g2_delta", "d_g2_50", "d_g2_h",
         "u_g1", "d_g1") %>%
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
u_post <- data.frame(
  crange,
  sapply(sample_post, function(x) {
    b_curve(
      lcrange,
      params_posterior$u_g2_i[x],
      params_posterior$u_g2_delta[x],
      params_posterior$u_g2_h[x],
      params_posterior$u_g2_50[x]
    )
  })
) %>%
  pivot_longer(cols = !crange, names_to = "sample", values_to = "u")
d_post <- data.frame(
  crange,
  sapply(sample_post, function(x) {
    d_curve(
      lcrange,
      params_posterior$d_g2_0[x],
      params_posterior$d_g2_delta[x],
      params_posterior$d_g2_h[x],
      params_posterior$d_g2_50[x]
    )
  })
) %>%
  pivot_longer(cols = !crange, names_to = "sample", values_to = "d")
posterior_samples <- u_post %>%
  inner_join(d_post) %>%
  mutate(crange = ifelse(crange == 0, 0.000001, crange))
posterior_samples %>%
  ggplot(aes(x = crange)) +
  geom_line(aes(y = u, group = sample), color = "darkblue", alpha = 0.1) +
  geom_line(aes(y = d, group = sample), color = "darkred", alpha = 0.1) +
  geom_point(data = curve_truth, aes(x = concentration, y = u_g2), color = "blue", size = 2) +
  geom_point(data = curve_truth, aes(x = concentration, y = d_g2), color = "red", size = 2) +
  scale_x_log10() +
  theme_bw()
