### Logistic Concentration of Birth, Death, Clearance Rates
### No Error
### Mixed Effects

################################################################################

library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters: make b50 = d50 and bh = dh create proper gr curve
# Generate data -------------------------------------------------
# Generate Parameters: make b50 = d50 and bh = dh create proper gr curve
set.seed(358811787)
b0 <- rnorm(7, 1, 0.05)
bi <- rnorm(7, 0.2, 0.07)
db <- b0 - bi
b50 <- rnorm(7, 0.02, 0.005)
logb50 <- log(b50) # convert to logscale
bh <- rnorm(7, 2, 0.1)
d0 <- rnorm(7, 0.3, 0.05)
di <- rnorm(7, 0.5, 0.05)
dd <- di - d0
d50 <- rnorm(7, 1, 0.05)
logd50 <- log(d50) # convert to log scale
dh <- bh

kr <- rnorm(7, 0.05, 0.01)

param_truth <- data.frame(
  group = 1:length(b0),
  bi = bi, db = db, b50 = logb50, bh = bh,
  d0 = d0, dd = dd, d50 = logd50, dh = dh,
  k=kr
)

concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)

b_curve <- function(lc, bi, db, bh, lb50) bi + db / (1 + exp(bh * (lc - lb50)))
d_curve <- function(lc, d0, dd, dh, ld50) d0 + dd - dd / (1 + exp(dh * (lc - ld50)))

b_truth <- data.frame(sapply(1:length(bi), function(i) b_curve(log_conc, bi[i], db[i], bh[i], logb50[i])))
d_truth <- data.frame(sapply(1:length(d0), function(i) d_curve(log_conc, d0[i], dd[i], dh[i], logd50[i])))
b_truth$concentration <- concentrations
d_truth$concentration <- concentrations
b_truth <- b_truth %>% pivot_longer(!concentration, names_to = "group", values_to = "b")
d_truth <- d_truth %>% pivot_longer(!concentration, names_to = "group", values_to = "d")
curve_truth <- b_truth %>%
  inner_join(d_truth) %>%
  mutate(group = sub("X", "", group))

curve_truth %>% ggplot(aes(x = concentration, group = group)) +
  geom_point(aes(y = b), color = "blue", size = 2) +
  geom_point(aes(y = d), color = "red", size = 2) +
  geom_line(aes(y = b), color = "blue") +
  geom_line(aes(y = d), color = "red") +
  scale_x_log10()


N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(1000, 0) # Number of ancestors


# Simulate process with estipop
set.seed(3768237)
dat <- data.frame()
for (i in 1:nrow(curve_truth)) {
  params <- c(curve_truth$b[i], curve_truth$d[i], kr[i])
  # time-homogeneous two-type model
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = c(2,0)), # birth
    transition(rate = rate(params[2]), parent = 1, offspring = c(0,1)), # death
    transition(rate = rate(params[2]), parent = 2, offspring = c(0,0)) # clearance
  )

  dat_temp <- branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- curve_truth$concentration[i]
  dat_temp$log_conc <- log(dat_temp$concentration)
  dat_temp$group <- curve_truth$group[i]
  dat <- bind_rows(dat, dat_temp)
}

# Plot to look at data
dat %>%
  pivot_longer(c(type1, type2), names_to = "type") %>%
  ggplot(aes(x = time, group = interaction(type, group, as.factor(rep)))) +
  geom_line(aes(y = value, color = group, linetype = type)) +
  facet_wrap(concentration ~ .) +
  theme_bw() +
  scale_y_log10()


dat <- dat %>%
  rename(sample = rep) %>%
  group_by(sample, concentration, group) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    type1_prev = lag(type1),
    type2_prev = lag(type2)
  )

dat <- dat %>% filter(!is.na(dt))

################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/2type-livedead/logistic/birthdeathclearance_mixedeffects.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
ntypes <- 2
nevents <- 3
dt <- array(unique(dat$dt))

# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)
num_groups <- length(unique(dat$group))

# 4-parameter logistic Birth-Death - growth rate ME no error ------------------------------------
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
  conc = log_conc,
  g_idx = as.numeric(as.factor(dat$group)),
  ng = num_groups
)
model.list <- list(
  e_mat = matrix(c(2, 0, 0, 1, 0, 0), nrow = 3, ncol = 2, byrow = T),
  p_vec = c(1, 1, 2)
)
prior.list <- list(
  theta_1_mu_bi = 1.0,
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
  theta_2_s_dh = 0.2,
  theta_1_mu_k = 0.05,
  theta_2_mu_k = 0.03,
  theta_1_s_k = 0.02,
  theta_2_s_k = 0.01
)
initfun <- function() {
  list(
    bi = runif(num_groups, 0.01, 1),
    db = runif(num_groups, 0.1, 1),
    b50 = runif(num_groups, -2, 1),
    bh = runif(num_groups, 0.1, 1),
    d0 = runif(num_groups, 0.01, 1),
    dd = runif(num_groups, 0.1, 1),
    d50 = runif(num_groups, -2, 1),
    dh = runif(num_groups, 0.1, 1),
    k = runif(num_groups, 0.01, 0.1),
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

model_posterior_pf <- stanmodel$pathfinder(
  data = c(data.list, model.list, prior.list),
  refresh = 1,
  init = initfun,
  draws = 1000
)


params_posterior <- model_posterior$draws(
  variables =
    c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh", "k"), format = "df"
)
hyperparams_posterior <- model_posterior$draws(
  variables =
    c(
      "mu_bi", "mu_db", "mu_b50", "mu_bh", "mu_d0", "mu_dd", "mu_d50", "mu_dh",
      "s_bi", "s_db", "s_b50", "s_bh", "s_d0", "s_dd", "s_d50", "s_dh"
    ), format = "df"
)

################################################################################
# Compare to ground truth


params_posterior_l <- params_posterior %>%
  select(.chain, .draw, starts_with("b"), starts_with("d")) %>%
  pivot_longer(!c(".chain", ".draw"), names_to = "stat") %>%
  separate(stat, sep = "\\[", c("stat", "group")) %>%
  mutate(group = (sub("\\]", "", group)))
head(params_posterior_l)

ground_truth_l <- param_truth %>%
  pivot_longer(!group, names_to = "stat") %>%
  mutate(group = as.character(group))
params_posterior_l %>%
  ggplot(aes(x = value, color = group, group = interaction(group, .chain))) +
  geom_density() +
  geom_vline(data = ground_truth_l, aes(xintercept = value), linewidth = 2) +
  facet_wrap(stat ~ group, scales = "free", ncol = 7) +
  theme_bw()

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
  ggplot(aes(x = crange)) +
  geom_line(aes(y = b, group = sample), color = "darkblue", alpha = 0.1) +
  geom_line(aes(y = d, group = sample), color = "darkred", alpha = 0.1) +
  geom_point(data = curve_truth, aes(x = concentration, y = b), color = "blue", size = 2) +
  geom_point(data = curve_truth, aes(x = concentration, y = d), color = "red", size = 2) +
  facet_wrap(group ~ .) +
  scale_x_log10() +
  theme_bw()
