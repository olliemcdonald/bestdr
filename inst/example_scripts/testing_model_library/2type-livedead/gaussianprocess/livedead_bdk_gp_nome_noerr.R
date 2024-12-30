### Gaussian Process Concentration of Birth, Death, Clearance Rates
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

k <- 0.005

param_truth <- data.frame(
  bi = bi, b50 = b50, bh = bh, db = db,
  d0 = d0, d50 = d50, dh = dh, dd = dd,
  k = k
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
  geom_hline(aes(yintercept = param_truth$k), color = "darkgreen") +
  scale_x_log10() +
  labs(x = "Concentration", y = "True Values") +
  theme_bw()


N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(1000, 0) # Number of ancestors


# Simulate process with estipop
set.seed(77030)
dat <- data.frame()
for (i in seq_along(concentrations)) {
  params <- c(b_truth[i], d_truth[i], k)
  # time-homogeneous two-type model
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = c(2, 0)), # birth
    transition(rate = rate(params[2]), parent = 1, offspring = c(0, 1)), # death
    transition(rate = rate(params[3]), parent = 2, offspring = c(0, 0)) # clearance
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
stanfileloc <- system.file("model_library/2type-livedead/gaussianprocess/birthdeathclearance.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
ntypes <- 2
nevents <- 3
dt <- array(unique(dat$dt))

# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
dim(log_conc) <- c(length(log_conc), 1)
num_conc <- length(log_conc)
log_conc[log_conc == -Inf] <- -12

# Concentrations for prediction (in between values allowed)
# concentrations_prediction <- c(10^seq(-4, 4, 0.2))
# concentrations_prediction <- setdiff(concentrations_prediction, concentrations)
# log_conc_prediction <- log(concentrations_prediction)
# nc_predict <- length(log_conc_prediction)
# dim(log_conc_prediction) <- c(length(log_conc_prediction), 1)


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
  nd = 1,
  conc = log_conc
)
model.list <- list(
  e_mat = matrix(c(2, 0, 0, 1, 0, 0), nrow = 3, ncol = 2, byrow = T),
  p_vec = c(1, 1, 2)
)
prior.list <- list(
  alpha_rho = 2,
  beta_rho = 2,
  mu_alpha_b = 0,
  s_alpha_b = 0.2,
  mu_alpha_d = 0,
  s_alpha_d = 0.2,
  mu_k = 0.05,
  s_k = 0.02
  )

initfun <- function() {
  list(
    birth_tilde = runif(num_conc, -1, 1),
    death_tilde = runif(num_conc, -1, 1),
    count_err = matrix(runif(data.list$N * ntypes, -0.5, 0.5), ncol = ntypes)
  )
}


model_posterior <- stanmodel$sample(
  data = c(data.list, model.list, prior.list),
  chains = 4,
  parallel_chains = 4,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)

################################################################################
# Compare to ground truth
params_posterior <- model_posterior$draws(variables = c("b", "d", "k"), format = "df")


params_posterior_l <- params_posterior %>%
  select(.chain, .draw, starts_with("b"), starts_with("d")) %>%
  pivot_longer(!c(".chain", ".draw"), names_to = "stat") %>%
  separate(stat, sep = "\\[", c("stat", "conc_idx")) %>%
  mutate(conc_idx = (sub("\\]", "", conc_idx))) %>%
  left_join(data.frame(
    conc_idx = as.character(c(1:(num_conc + nc_predict))),
    concentration = c(concentrations, concentrations_prediction),
    log_conc = c(log_conc, log_conc_prediction)
  ))


curve_truth_l <- curve_truth %>%
  pivot_longer(!c("concentration", "log_conc"), names_to = "stat") %>%
  mutate(concentration = ifelse(concentration == 0, exp(-12), concentration))


params_posterior_l %>%
  mutate(concentration = ifelse(concentration == 0, exp(-12), concentration)) %>%
  ggplot(aes(x = concentration, y = value)) +
  geom_violin(aes(group = interaction(stat, conc_idx), fill = stat), position = "identity", alpha = 0.05) +
  geom_point(data = curve_truth_l, aes(color = stat), size = 2) +
  scale_x_log10()

params_posterior_l %>%
  filter(.draw %in% sample(nrow(params_posterior), 300)) %>%
  mutate(concentration = ifelse(concentration == 0, exp(-12), concentration)) %>%
  arrange(concentration) %>%
  ggplot(aes(x = concentration, y = value)) +
  geom_line(aes(group = interaction(.draw, stat), color = stat), alpha = 0.1) +
  geom_point(data = curve_truth_l, size = 2) +
  scale_x_log10()
