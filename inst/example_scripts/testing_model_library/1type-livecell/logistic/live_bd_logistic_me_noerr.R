### Logistic Concentration of Birth and Death Rates
### No Error
### Mixed Effects

################################################################################
library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters: make b50 = d50 and bh = dh create proper gr curve
set.seed(12749873)

# parameters
mu_b0_1 <- 1
s_b0_1 <- 0.04
mu_bi_1 <- 0.2
s_bi_1 <- 0.03
mu_b50_1 <- 0.2
s_b50_1 <- 0.04
mu_bh_1 <- 2
s_bh_1 <- 0.3
mu_d0_1 <- 0.3
s_d0_1 <- 0.05
mu_di_1 <- 0.5
s_di_1 <- 0.03
mu_d50_1 <- 01
s_d50_1 <- 0.1
mu_dh_1 <- 2
s_dh_1 <- 0.5


b0 <- rnorm(7, mu_b0_1, s_b0_1)
bi <- rnorm(7, mu_bi_1, s_bi_1)
db <- b0 - bi
b50 <- rnorm(7, mu_b50_1, s_b50_1)
logb50 <- log(b50) # convert to logscale
bh <- rnorm(7, mu_bh_1, s_bh_1)
d0 <- rnorm(7, mu_d0_1, s_d0_1)
di <- rnorm(7, mu_di_1, s_di_1)
dd <- di - d0
d50 <- rnorm(7, mu_d50_1, s_d50_1)
logd50 <- log(d50) # convert to log scale
dh <- rnorm(7, mu_dh_1, s_dh_1)


param_truth <- data.frame(
  group = 1:length(b0),
  bi = bi, db = db, b50 = logb50, bh = bh,
  d0 = d0, dd = dd, d50 = logd50, dh = dh
)

concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)

b_curve <- function(lc, bi, db, bh, lb50) bi + db / (1 + exp(bh * (lc - lb50)))
d_curve <- function(lc, d0, dd, dh, ld50) d0 + dd - dd / (1 + exp(dh * (lc - ld50)))

b_truth <- data.frame(sapply(1:length(bi), function(i) b_curve(log_conc, bi[i], db[i], bh[i], logb50[i])))
d_truth <- data.frame(sapply(1:length(d0), function(i) d_curve(log_conc, d0[i], dd[i], dh[i], logd50[i])))
b_truth$concentration <- concentrations
d_truth$concentration <- concentrations
b_truth <- b_truth %>% pivot_longer(!concentration, names_to="group", values_to="b")
d_truth <- d_truth %>% pivot_longer(!concentration, names_to = "group", values_to="d")
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
init_pop <- 200 # Number of ancestors


# Simulate process with estipop
set.seed(453768423)
dat <- data.frame()
for (i in 1:nrow(curve_truth)) {
    params <- c(curve_truth$b[i], curve_truth$d[i])
    # time-homogeneous two-type model
    model <- process_model(
      transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
      transition(rate = rate(params[2]), parent = 1, offspring = 0) # death
    ) # clearance

    dat_temp <- branch(model, params, init_pop, times, N, silent = T)
    dat_temp$concentration <- curve_truth$concentration[i]
    dat_temp$log_conc <- log(dat_temp$concentration)
    dat_temp$group <- curve_truth$group[i]
    dat <- bind_rows(dat, dat_temp)
  }


dat <- dat %>%
  rename(sample = rep, live = type1) %>%
  group_by(sample, concentration, group) %>%
  mutate(
    dt = signif(time - lag(time), 6),
    live_prev = lag(live)
  )

# Plot to look at data
dat %>% ggplot(aes(x = time, group = interaction(group, as.factor(sample)))) +
  geom_line(aes(y = live, color = group)) +
  facet_wrap(concentration ~ .) +
  theme_bw()

dat <- dat %>% filter(!is.na(dt))

################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/logistic/birthdeath_mixedeffects.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)

# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
num_conc <- length(log_conc)
num_groups <- length(unique(dat$group))

# 4-parameter logistic Birth-Death - growth rate ME no error ------------------------------------
data.list <- list(
  N = nrow(dat),
  count = dat$live,
  count_prev = dat$live_prev,
  dt = dat$dt,
  c_idx = dat$log_conc_idx,
  nc = num_conc,
  conc = log_conc,
  g_idx = as.numeric(as.factor(dat$group)),
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



model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 4,
  parallel_chains = 4,
  refresh = 10,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)


params_posterior <- model_posterior$draws(
  variables =
    c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh"), format = "df"
)
hyperparams_posterior <- model_posterior$draws(
  variables =
    c("mu_bi", "mu_db", "mu_b50", "mu_bh", "mu_d0", "mu_dd", "mu_d50", "mu_dh",
      "s_bi", "s_db", "s_b50", "s_bh", "s_d0", "s_dd", "s_d50", "s_dh"), format = "df"
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
  facet_wrap(stat ~ group, scales = "free", ncol=7) +
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
  facet_wrap(group~.) +
  scale_x_log10() +
  theme_bw()

posterior_summ <- posterior_samples %>%
  ungroup %>%
  group_by(group, crange) %>%
  summarize(
    b_mean = mean(b),
    d_mean = mean(d),
    b_sd = sd(b),
    d_sd = sd(d)
  )

posterior_summ %>%
  ggplot(aes(x = crange)) +
  geom_line(aes(y = b_mean), color = "darkblue", alpha = 0.5) +
  geom_ribbon(aes(ymin = b_mean - 2 * b_sd, ymax = b_mean + 2 * b_sd), fill = "blue", alpha = 0.05) +
  geom_line(aes(y = d_mean), color = "darkred", alpha = 0.1) +
  geom_ribbon(aes(ymin = d_mean - 2 * d_sd, ymax = d_mean + 2 * d_sd), fill = "darkred", alpha = 0.05) +
  geom_point(data = curve_truth, aes(x = concentration, y = b), color = "blue", size = 2) +
  geom_point(data = curve_truth, aes(x = concentration, y = d), color = "red", size = 2) +
  facet_wrap(group ~ .) +
  scale_x_log10() +
  theme_bw()
  