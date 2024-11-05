### Single Dose Estimation of Birth and Death Rates
### No Error
### Mixed Effects

################################################################################
library(tidyverse)
library(estipop)
library(cmdstanr)

################################################################################
# Generate data -------------------------------------------------
# Generate Parameters:
set.seed(3762)
ngroups <- 7
br <- rnorm(ngroups, 1, 0.1) # birth rate
dr <- rnorm(ngroups, 0.5, 0.1)
lambda <- br - dr
N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 200 # Number of ancestors

# Simulate process with estipop
ground_truth <- data.frame(group = 1:length(br), b = br, d = dr)
dat <- data.frame()
for (i in 1:nrow(ground_truth)) {
  params <- c(ground_truth$b[i], ground_truth$d[i])
  model <- process_model(
    transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
    transition(rate = rate(params[2]), parent = 1, offspring = 0) # death
  ) # death
  dat_temp <- branch(model, params, init_pop, times, N, silent = T)

  dat_temp <- dat_temp %>%
    rename(sample = rep, live = type1) %>%
    group_by(sample) %>%
    mutate(
      group = i,
      dt = signif(time - lag(time), 6),
      live_prev = lag(live)
    )
  dat <- dat %>% bind_rows(dat_temp)
}
rm(dat_temp)

dat <- dat %>% filter(!is.na(dt))


dat %>% ggplot(aes(x = time, y = live, color = as.factor(group), group = interaction(group, sample))) +
  geom_line()


################################################################################
# Estimation with cmdstanr
stanfileloc <- system.file("model_library/1type-livecell/singledose/birthdeath_mixedeffects.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)


data.list <- list(
  N = nrow(dat),
  count = dat$live,
  count_prev = dat$live_prev,
  dt = dat$dt,
  g = dat$group,
  ng = ngroups
)
prior.list <- list(
  theta_1_mu_b = 1,
  theta_2_mu_b = 0.5,
  theta_1_mu_d = 1,
  theta_2_mu_d = 0.5,
  theta_1_s_b = 0,
  theta_2_s_b = 0.3,
  theta_1_s_d = 0,
  theta_2_s_d = 0.3
)

initfun <- function() {
  list(
    b = runif(ngroups, 0, 1),
    d = runif(ngroups, 0, 0.5),
    count_err = runif(data.list$N, -0.5, 0.5),
    mu_b = runif(1, 0, 1),
    mu_d = runif(1, 0, 0.5),
    s_b = runif(1, 0, 1),
    s_d = runif(1, 0, 1)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 4,
  parallel_chains=4,
  init = initfun,
  refresh=1,
  iter_warmup = 1000,
  iter_sampling = 1000)

hyperparams_posterior <- model_posterior$draws(variables = c("mu_b", "mu_d", "s_b", "s_d"), format = "df")
params_posterior <- model_posterior$draws(variables = c("b", "d"), format = "df")

################################################################################
# Compare to ground truth


params_posterior_l <- params_posterior %>%
  select(.chain, starts_with("b"), starts_with("d")) %>%
  pivot_longer(!c(".chain"), names_to = "stat") %>%
  separate(stat, sep = "\\[", c("stat", "group")) %>%
  mutate(group = (sub("\\]", "", group)))

ground_truth_l <- ground_truth %>%
  pivot_longer(!group, names_to = "stat") %>%
  mutate(group = as.character(group))
params_posterior_l %>%
  ggplot(aes(x = value, color=group, group = interaction(group, .chain))) +
  geom_density() +
  geom_vline(data = ground_truth_l, aes(xintercept = value), linewidth = 2) +
  facet_grid(group ~ stat, scales = "free") +
  theme_bw()


params_posterior_l %>%
  left_join(ground_truth) %>%
  mutate(diff = value - true_value) %>%
  group_by(.chain, stat) %>%
  summarize(mean_diff = mean(diff)) %>%
  pivot_wider(names_from = stat, values_from = mean_diff)