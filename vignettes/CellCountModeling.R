## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width=7,
  fig.asp = 0.8,
  cache = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bestdr)
library(tidyverse)
library(cmdstanr)
library(estipop)

## ----data-generation----------------------------------------------------------
# Generate curve from 4-parameter logistic function for birth and death
b0 <- 1 # left asymptote
bi <- 0.2 # right asymptote
b50 <- -5 # log-dose at midpoint of birthrate (should switch to log scale)
bh <- 2 # Hill Coefficient
db <- b0 - bi # Difference in asymptotes (for parameterization)

d0 <- 0.01
di <- 0.5
d50 <- 1
dh <- 2
dd <- di - d0

b_curve <- function(lc, bi, db, bh, b50) bi + db / (1 + exp(bh * (lc - b50)))
d_curve <- function(lc, d0, dd, dh, d50) d0 + dd - dd / (1 + exp(dh * (lc - d50)))

param_truth <- data.frame(
  bi = bi, b50 = b50, bh = bh, db = db,
  d0 = d0, d50 = d50, dh = dh, dd = dd
)

crange <- c(0, 10^seq(-4, 4, 0.1))
lc_range <- log(crange)
concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)

b_truth <- b_curve(lc_range, bi, db, bh, b50)
d_truth <- d_curve(lc_range, d0, dd, dh, d50)
g_truth <- b_truth - d_truth

curve_truth <- data.frame(concentration = crange,
                          log_conc = lc_range,
                          b = b_truth, d = d_truth, g = g_truth) %>%
  pivot_longer(-c(concentration, log_conc), names_to = "Rate", values_to = "Value")


curve_truth %>%
  ggplot(aes(x = concentration, y = Value, color = Rate)) +
  geom_line() +
  geom_point(data = curve_truth %>% filter(concentration %in% concentrations)) +
  scale_x_log10() +
  labs(x = "Concentration", y = "Rate") +
  scale_color_manual(values = c("darkblue", "darkred", "purple")) +
  theme_bw()

## ----simulate-----------------------------------------------------------------
N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors

rate_truth <- data.frame(concentration = crange,
                         log_conc = lc_range,
                         b = b_truth, d = d_truth, g = g_truth) %>%
  filter(concentration %in% concentrations)
# Simulate process with estipop
set.seed(9283)
dat <- data.frame()
for (i in seq_along(rate_truth$concentration)) {
  params <- c(rate_truth$b[i], rate_truth$d[i])
  # time-homogeneous two-type model
  model <- estipop::process_model(
    estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = 2), # birth
    estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = 0) # death
  ) # clearance
  
  dat_temp <- estipop::branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- concentrations[i]
  dat_temp$log_conc <- log_conc[i]
  dat <- bind_rows(dat, dat_temp)
}

head(dat)

# Plot to look at data
dat %>% ggplot(aes(x = time, group = as.factor(rep))) +
  geom_line(aes(y = type1), col = "black") +
  facet_wrap(concentration ~ .) +
  theme_bw() +
  scale_y_log10()

## ----preprocess---------------------------------------------------------------
dat_dt <- dat %>%
  group_by(rep, concentration) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    type1_prev = lag(type1)
) %>% filter(!is.na(dt))

## ----estimation-preprocess----------------------------------------------------
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/logistic/birthdeath.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
dat_dt$log_conc_idx <- as.numeric(factor(dat_dt$log_conc, ordered = T))
log_conc <- unique(dat_dt$log_conc)
num_conc <- length(log_conc)

data.list <- list(
  N = nrow(dat_dt),
  count = dat_dt$type1,
  count_prev = dat_dt$type1_prev,
  dt = dat_dt$dt,
  c_idx = dat_dt$log_conc_idx,
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

## ----sampler, warning=FALSE---------------------------------------------------
model_posterior <- stanmodel$sample(
  seed = 117,
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 0,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 500
)

## ----posterior-summary, warning=FALSE-----------------------------------------
model_posterior$summary()
params_posterior <- model_posterior$draws(variables = c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh"), format = "df")

## ----comparison, warning=FALSE------------------------------------------------
# Compare to ground truth
param_truth_l <- param_truth %>% pivot_longer(everything(), names_to = "stat")

params_posterior_l <- params_posterior %>%
  select(.chain, bi, db, b50, bh, d0, dd, d50, dh) %>%
  pivot_longer(!c(".chain"), names_to = "stat")

# Marginal Posteriors
params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = param_truth_l, aes(xintercept = value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()

## ----setup-posterior-curves---------------------------------------------------
n_posteriors <- 500
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
  mutate(crange = ifelse(crange == 0, 0.000001, crange),
         g = b - d)
head(posterior_samples)

## ----plotting-curves----------------------------------------------------------
posterior_samples %>% pivot_longer(-c(crange, sample),
                                   names_to = "param", values_to = "value") %>%
  ggplot(aes(x = crange, y = value, color = param, group = interaction(param, sample))) +
  geom_line(alpha = 0.1, linewidth = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c("darkblue", "darkred", "purple")) +
  theme_bw()

## ----hier-params, warning=FALSE-----------------------------------------------
# parameters
set.seed(1238)
b0 <- rnorm(5, 1, 0.05)
bi <- rnorm(5, 0.2, 0.02)
db <- b0 - bi
b50 <- rnorm(5, -5, 0.5)
bh <- rnorm(5, 2, 0.3)
d0 <- rnorm(5, 0.05, 0.02)
di <- rnorm(5, 0.5, 0.1)
dd <- di - d0
d50 <- rnorm(5, 0, 0.2)
dh <- rnorm(5, 2, 0.4)


param_truth <- data.frame(
  group = 1:length(b0),
  bi = bi, db = db, b50 = b50, bh = bh,
  d0 = d0, dd = dd, d50 = d50, dh = dh
)

concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)

b_curve <- function(lc, bi, db, bh, lb50) bi + db / (1 + exp(bh * (lc - lb50)))
d_curve <- function(lc, d0, dd, dh, ld50) d0 + dd - dd / (1 + exp(dh * (lc - ld50)))

b_truth <- data.frame(sapply(1:length(bi), function(i) b_curve(log_conc, bi[i], db[i], bh[i], b50[i])))
d_truth <- data.frame(sapply(1:length(d0), function(i) d_curve(log_conc, d0[i], dd[i], dh[i], d50[i])))
b_truth$concentration <- concentrations
d_truth$concentration <- concentrations
b_truth <- b_truth %>% pivot_longer(!concentration, names_to="group", values_to="b")
d_truth <- d_truth %>% pivot_longer(!concentration, names_to = "group", values_to="d")
curve_truth <- b_truth %>%
  inner_join(d_truth) %>%
  mutate(group = sub("X", "", group),
         g = b - d)

curve_truth %>%
  pivot_longer(-c(concentration, group), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = concentration, y = value, color = param, group = interaction(param, group))) +
    geom_line() +
    geom_point(size = 2) +
  scale_color_manual(values = c("darkblue", "darkred", "purple")) +
  scale_x_log10() +
  theme_bw()

## ----hier-simulate------------------------------------------------------------
N <- 5 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors


# Simulate process with estipop
set.seed(342)
dat <- data.frame()
for (i in 1:nrow(curve_truth)) {
  params <- c(curve_truth$b[i], curve_truth$d[i])
  # time-homogeneous two-type model
  model <- estipop::process_model(
    estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = 2), # birth
    estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = 0) # death
  ) # clearance
  
  dat_temp <- estipop::branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- curve_truth$concentration[i]
  dat_temp$log_conc <- log(dat_temp$concentration)
  dat_temp$group <- curve_truth$group[i]
  dat <- bind_rows(dat, dat_temp)
}
head(dat)

dat %>% ggplot(aes(x = time, group = interaction(group, as.factor(rep)))) +
  geom_line(aes(y = type1, color = group), alpha = 0.5) +
  facet_wrap(concentration ~ ., scales = "free_y") +
  scale_y_log10() +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#a6761d")) +
  theme_bw()

## ----hier-preproc-------------------------------------------------------------
# Preprocessing
dat_dr <- dat %>%
  group_by(rep, concentration, group) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    type1_prev = lag(type1)
  ) %>%
  filter(!is.na(dt))

# Set up the data for input into Stan
dat_dr$log_conc_idx <- as.numeric(factor(dat_dr$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)
num_groups <- length(unique(dat$group))

## ----hier-estimate, warning=FALSE, message=FALSE------------------------------
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/logistic/birthdeath_mixedeffects.stan",
                           package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

data.list <- list(
  N = nrow(dat_dr),
  count = dat_dr$type1,
  count_prev = dat_dr$type1_prev,
  dt = dat_dr$dt,
  c_idx = dat_dr$log_conc_idx,
  nc = num_conc,
  conc = log_conc,
  g_idx = as.numeric(as.factor(dat_dr$group)),
  ng = num_groups
)
prior.list <- list(theta_1_mu_bi = 1.0,
                   theta_2_mu_bi = 0.4,
                   theta_1_s_bi = 0.5,
                   theta_2_s_bi = 0.2,
                   theta_1_mu_db = 1.0,
                   theta_2_mu_db = 0.4,
                   theta_1_s_db = 0.5,
                   theta_2_s_db = 0.2,
                   theta_1_mu_b50 = 0.0,
                   theta_2_mu_b50 = 0.5,
                   theta_1_s_b50 = 0.5,
                   theta_2_s_b50 = 0.2,
                   theta_1_mu_bh = 1.0,
                   theta_2_mu_bh = 0.5,
                   theta_1_s_bh = 0.5,
                   theta_2_s_bh = 0.2,
                   theta_1_mu_d0 = 0.1,
                   theta_2_mu_d0 = 0.4,
                   theta_1_s_d0 = 0.5,
                   theta_2_s_d0 = 0.2,
                   theta_1_mu_dd = 1.0,
                   theta_2_mu_dd = 0.4,
                   theta_1_s_dd = 0.5,
                   theta_2_s_dd = 0.2,
                   theta_1_mu_d50 = 0.0,
                   theta_2_mu_d50 = 0.5,
                   theta_1_s_d50 = 0.5,
                   theta_2_s_d50 = 0.2,
                   theta_1_mu_dh = 1.0,
                   theta_2_mu_dh = 0.5,
                   theta_1_s_dh = 0.5,
                   theta_2_s_dh = 0.2)
initfun = function() {
  list(
    bi = runif(num_groups, 0.1, 2),
    db = runif(num_groups, 0.1, 2),
    b50 = runif(num_groups, -2, 2),
    bh = runif(num_groups, 0.1, 2),
    d0 = runif(num_groups, 0.1, 2),
    dd = runif(num_groups, 0.1, 2),
    d50 = runif(num_groups, -2, 2),
    dh = runif(num_groups, 0.1, 2),
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

# Sampling
model_posterior <- stanmodel$sample(
  seed=711,
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 0,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 500
)

## ----hier-results-------------------------------------------------------------
params_posterior <- model_posterior$draws(
  variables =
    c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh"), format = "df"
)
hyperparams_posterior <- model_posterior$draws(
  variables =
    c("mu_bi", "mu_db", "mu_b50", "mu_bh", "mu_d0", "mu_dd", "mu_d50", "mu_dh",
      "s_bi", "s_db", "s_b50", "s_bh", "s_d0", "s_dd", "s_d50", "s_dh"), format = "df"
)

params_posterior_l <- params_posterior %>%
  select(.chain, .draw, starts_with("b"), starts_with("d")) %>%
  pivot_longer(!c(".chain", ".draw"), names_to = "stat") %>%
  separate(stat, sep = "\\[", c("stat", "group")) %>%
  mutate(group = (sub("\\]", "", group)))

params_posterior_l %>% group_by(stat, group) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(names_from = "stat", values_from = "value")

ground_truth_l <- param_truth %>%
  pivot_longer(!group, names_to = "stat") %>%
  mutate(group = as.character(group))
params_posterior_l %>%
  ggplot(aes(x = value, color = group, group = interaction(group, .chain))) +
  geom_density() +
  geom_vline(data = ground_truth_l, aes(xintercept = value, color = group), linewidth = 2) +
  facet_wrap(stat ~ ., scales = "free") +
  theme_bw()

## ----hier-posteriorcurves-----------------------------------------------------
# Samples of posterior curves
n_posteriors <- 300
crange <- c(10^seq(-6, 4, 0.1))
lcrange <- log(crange)
posterior_samples <- data.frame()
for (i in 1:num_groups) {
  sample_post <- sample(nrow(params_posterior), n_posteriors)
  posterior_subset <- params_posterior_l %>%
    filter(group == i, .draw %in% sample_post) %>%
    select(-c(.chain, group)) %>%
    pivot_wider(id_cols = .draw, names_from = "stat")
  b_subset <- data.frame(
    group = i,
    crange,
    sapply(1:n_posteriors, function(x) {
      b_curve(
        lcrange,
        posterior_subset$bi[x],
        posterior_subset$db[x],
        posterior_subset$bh[x],
        posterior_subset$b50[x]
      )
    })
  ) %>%
    pivot_longer(cols = !c(group, crange), names_to = "sample", values_to = "b")
  d_subset <- data.frame(
    group = i,
    crange,
    sapply(1:n_posteriors, function(x) {
      d_curve(
        lcrange,
        posterior_subset$d0[x],
        posterior_subset$dd[x],
        posterior_subset$dh[x],
        posterior_subset$d50[x]
      )
    })
  ) %>%
    pivot_longer(cols = !c(group, crange), names_to = "sample", values_to = "d")
  posterior_samples <- posterior_samples %>%
    bind_rows(b_subset %>%
                inner_join(d_subset) %>%
                mutate(crange = ifelse(crange == 0, 0.000001, crange)))
}


posterior_samples %>%
  pivot_longer(-c(group, crange, sample), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = crange, y = value, color = param, group = interaction(param, sample))) +
  geom_line() +
  facet_grid(.~group) +
  scale_color_manual(values = c("darkblue", "darkred")) +
  scale_x_log10() +
  theme_bw()

