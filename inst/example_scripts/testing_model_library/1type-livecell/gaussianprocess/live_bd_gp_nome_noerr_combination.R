### Gaussian Process of Birth and Death Rates from Combination Model
### No Error
### No Mixed Effects
################################################################################

library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters:
# Generate curve from 2 separate sigmoidal functions with an interaction effect

# Start with additive model for ease
b0_1 <- 1 # left asymptote
bi_1 <- 0.2 # right asymptote
b50_1 <- -3 # log-dose at midpoint of birthrate (should switch to log scale)
bh_1 <- 0.5 # Hill Coefficient
db_1 <- b0_1 - bi_1 # Difference in asymptotes (for parameterization)

d0_1 <- 0.01
di_1 <- 0.5
d50_1 <- 0
dh_1 <- 1
dd_1 <- di_1 - d0_1

b0_2 <- 1 # left asymptote
bi_2 <- 0.5 # right asymptote
b50_2 <- -1 # log-dose at midpoint of birthrate (should switch to log scale)
bh_2 <- 1.0 # Hill Coefficient
db_2 <- b0_2 - bi_2 # Difference in asymptotes (for parameterization)

d0_2 <- 0.01
di_2 <- 0.4
d50_2 <- 2
dh_2 <- 1
dd_2 <- di_2 - d0_2

param_truth <- data.frame(
  bi_1 = bi_1, b50_1 = b50_1, bh_1 = bh_1, db_1 = db_1,
  d0_1 = d0_1, d50_1 = d50_1, dh_1 = dh_1, dd_1 = dd_1,
  bi_2 = bi_2, b50_2 = b50_2, bh_2 = bh_2, db_2 = db_2,
  d0_2 = d0_2, d50_2 = d50_2, dh_2 = dh_2, dd_2 = dd_2
)

concentrations_1 <- c(0, 10^seq(-4, 4, 1))
concentrations_2 <- c(0, 10^seq(-4, 4, 1))
log_conc_1 <- log(concentrations_1)
log_conc_2 <- log(concentrations_2)
b_curve <- function(lc, bi, db, bh, b50) bi + db / (1 + exp(bh * (lc - b50)))
d_curve <- function(lc, d0, dd, dh, d50) d0 + dd - dd / (1 + exp(dh * (lc - d50)))

b_truth_1 <- b_curve(log_conc_1, bi_1, db_1, bh_1, b50_1)
d_truth_1 <- d_curve(log_conc_1, d0_1, dd_1, dh_1, d50_1)
b_truth_2 <- b_curve(log_conc_2, bi_2, db_2, bh_2, b50_2)
d_truth_2 <- d_curve(log_conc_2, d0_2, dd_2, dh_2, d50_2)

surface_truth <- expand.grid(concentrations_1, concentrations_2)
names(surface_truth) <- c("concentration_1", "concentration_2")
surface_truth <- surface_truth %>%
  mutate(
    log_conc_1 = log(concentration_1),
    log_conc_2 = log(concentration_2),
    b = b_curve(log_conc_1, bi_1, db_1, bh_1, b50_1) + b_curve(log_conc_2, bi_2, db_2, bh_2, b50_2),
    d = d_curve(log_conc_1, d0_1, dd_1, dh_1, d50_1) + d_curve(log_conc_2, d0_2, dd_2, dh_2, d50_2)
  )
surface_truth %>% ggplot(aes(log_conc_1, log_conc_2, z = b)) +
  geom_contour_filled()
surface_truth %>% ggplot(aes(log_conc_1, log_conc_2, z = d)) +
  geom_contour_filled()


N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors


# Simulate process with estipop
set.seed(1717171717)
dat <- data.frame()
for (i in 1:nrow(surface_truth)) {
  params <- c(surface_truth$b[i], surface_truth$d[i])
  # time-homogeneous two-type model
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
    transition(rate = rate(params[2]), parent = 1, offspring = 0) # death
  )

  dat_temp <- branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration_1 <- surface_truth$concentration_1[i]
  dat_temp$concentration_2 <- surface_truth$concentration_2[i]
  dat_temp$log_conc_1 <- log(surface_truth$concentration_1[i])
  dat_temp$log_conc_2 <- log(surface_truth$concentration_2[i])
  dat <- bind_rows(dat, dat_temp)
}

dat <- dat %>%
  rename(sample = rep, live = type1) %>%
  group_by(sample, concentration_1, concentration_2) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    live_prev = lag(live)
  )
dat <- dat %>% filter(!is.na(dt))
### Skip plotting for now since this will never show up

# Gaussian Process Birth-Death No ME no error ------------------------------------
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/gaussianprocess/birthdeath.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

unique_concs <- surface_truth %>%
  select(starts_with("concentration")) %>%
  mutate(conc_idx = row_number())

dat <- dat %>%
  left_join(unique_concs)

# Set up data for inputs
num_conc <- nrow(unique_concs)
log_conc <- cbind(
  log(unique_concs$concentration_1),
  log(unique_concs$concentration_2)
)
log_conc[log_conc == -Inf] <- -12


data.list <- list(
  N = nrow(dat),
  count = dat$live,
  count_prev = dat$live_prev,
  dt = dat$dt,
  c_idx = dat$conc_idx,
  nc = num_conc,
  nd = 2,
  conc = log_conc
)

prior.list <- list(
  alpha_rho = 2,
  beta_rho = 2,
  mu_alpha_b = 0,
  s_alpha_b = 0.1,
  mu_alpha_d = 0,
  s_alpha_d = 0.1
)
initfun <- function() {
  list(
    birth_tilde = runif(num_conc, -1, 1),
    death_tilde = runif(num_conc, -1, 1)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Extract and visualize posteriors ---------------------------------------------
params_posterior <- model_posterior$draws(variables = c("b", "d"), format = "df")


params_posterior_l <- params_posterior %>%
  select(.chain, .draw, starts_with("b"), starts_with("d")) %>%
  pivot_longer(!c(".chain", ".draw"), names_to = "stat") %>%
  separate(stat, sep = "\\[", c("stat", "conc_idx")) %>%
  mutate(conc_idx = (sub("\\]", "", conc_idx))) %>%
  left_join(unique_concs %>%
    mutate(
      conc_idx = as.character(conc_idx),
      log_conc_1 = ifelse(concentration_1 == 0, -12, log(concentration_1)),
      log_conc_2 = ifelse(concentration_2 == 0, -12, log(concentration_2)))
    )


surface_truth_l <- surface_truth %>%
  select(log_conc_1, log_conc_2, b, d) %>%
  pivot_longer(!c(log_conc_1, log_conc_2), names_to = "stat") %>%
  mutate(log_conc_1 = ifelse(log_conc_1 == -Inf, -12, log_conc_1),
         log_conc_2 = ifelse(log_conc_2 == -Inf, -12, log_conc_2))

params_posterior_l %>%
  filter(stat == "b") %>%
  ggplot(aes(x = log_conc_1, y = value, group = log_conc_1)) +
  geom_violin() +
  geom_point(data = surface_truth_l %>% filter(stat == "b")) +
  facet_wrap(log_conc_2 ~ ., scale = "free")

params_posterior_l %>%
  filter(stat == "d") %>%
  ggplot(aes(x = log_conc_1, y = value, group = log_conc_1)) +
  geom_violin() +
  geom_point(data = surface_truth_l %>% filter(stat == "d")) +
  facet_wrap(log_conc_2 ~ ., scale = "free")


params_posterior_l_summary <- params_posterior_l %>%
  group_by(stat, conc_idx, log_conc_1, log_conc_2) %>%
  summarize(
    mean = mean(value),
    error = 2 * sd(value),
    int95 = quantile(value, 0.95),
    int05 = quantile(value, 0.05)
  )

library(plotly)

plot_ly() %>%
  add_trace(
    data = (params_posterior_l_summary %>% filter(stat == "b")),
    x = ~log_conc_1,
    y = ~log_conc_2,
    z = ~mean,
    error_z = list(value = ~error, color = "black"),
    marker = list(color = ~mean),
    mode = "markers"
  ) %>%
  add_trace(
    data = (surface_truth_l %>% filter(stat == "b")),
    x = ~log_conc_1,
    y = ~log_conc_2,
    z = ~value,
    marker = list(size = 2, color="black")
  )


plot_ly() %>%
  add_trace(
    data = (params_posterior_l_summary %>% filter(stat == "d")),
    x = ~log_conc_1,
    y = ~log_conc_2,
    z = ~mean,
    error_z = list(value = ~error, color = "black"),
    marker = list(color = ~mean),
    mode = "markers"
  ) %>%
  add_trace(
    data = (surface_truth_l %>% filter(stat == "d")),
    x = ~log_conc_1,
    y = ~log_conc_2,
    z = ~value,
    marker = list(size = 2, color = "black")
  )

surface_posterior_b <- params_posterior_l_summary %>%
  ungroup() %>%
  filter(stat == "b") %>%
  arrange(log_conc_1, log_conc_2) %>%
  select(log_conc_1, log_conc_2, mean) %>%
  pivot_wider(names_from = log_conc_1, values_from = mean) %>% 
  select(-log_conc_2) %>%
  as.matrix()
surface_posterior_b <- list(x = c(-12, log_conc_1[-1]), y = c(-12, log_conc_1[-1]), z = surface_posterior$z)

surface_posterior_d <- params_posterior_l_summary %>%
  ungroup() %>%
  filter(stat == "d") %>%
  arrange(log_conc_1, log_conc_2) %>%
  select(log_conc_1, log_conc_2, mean) %>%
  pivot_wider(names_from = log_conc_1, values_from = mean) %>%
  select(-log_conc_2) %>%
  as.matrix()
surface_posterior_d <- list(x = c(-12, log_conc_1[-1]), y = c(-12, log_conc_1[-1]), z = surface_posterior$z)

plot_ly(x = surface_posterior_b$x, y = surface_posterior_b$y, z = surface_posterior_b$z) %>% add_surface()
plot_ly(x = surface_posterior_d$x, y = surface_posterior_d$y, z = surface_posterior_d$z) %>% add_surface()

