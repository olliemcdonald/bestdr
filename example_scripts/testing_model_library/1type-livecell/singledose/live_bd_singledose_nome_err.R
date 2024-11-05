### Single Dose Estimation of Birth and Death Rates
### w/Error
### No Mixed Effects

################################################################################

library(tidyverse)
library(estipop)
library(cmdstanr)

################################################################################
# Generate data -------------------------------------------------
# Generate Parameters:
br <- 0.2 # birth rate
dr <- 0.02 # death rate
obs_err <- 0.2
lambda <- br - dr
N <- 10 # Number of samples
tot_time <- 10
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- 500 # Number of ancestors

# Simulate process with estipop
set.seed(4733498)
params <- c(br, dr)
model <- process_model(
  transition(rate = rate(params[1]), parent = 1, offspring = 2), # birth
  transition(rate = rate(params[2]), parent = 1, offspring = 0)  # death
) 
dat <- branch(model, params, init_pop, times, N, silent = T)

dat <- dat %>%
  rename(sample = rep) %>%
  group_by(sample) %>%
  mutate(
    type1_prev = lag(type1),
    dt = signif(time - lag(time), 6),
    type1_mean = type1_prev * exp((br - dr) * dt),
    type1_sd = sqrt(type1_prev * (br+dr) / (br-dr) * (exp(2*(br-dr)*dt) - exp((br-dr)*dt))),
    error = rnorm(n(), mean = 0, sd = obs_err * type1_mean),
    count = ifelse(is.na(type1_prev), type1, round(type1 + error)),
    count_prev = lag(count),
    approx = rnorm(n(), mean = type1_mean, sd = sqrt(type1_sd^2 + (obs_err * type1_mean)^2))
  )

#
dat %>% ggplot(aes(x = time, group = sample)) +
  geom_line(aes(y = count)) +
  geom_line(aes(y = type1), color = "green") +
  #geom_line(aes(y = approx), color = "red") +
  scale_y_log10()

dat <- dat %>% filter(!is.na(dt))

################################################################################
# Estimation with cmdstanr
#stanmodel_noerr <- cmdstan_model("./stan/live_singledose/live_bd_singledose_nome_noerr.stan")
stanfileloc <- system.file("model_library/1type-livecell/singledose/birthdeath_error.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)
  
data.list <- list(
  N = nrow(dat),
  count = dat$count,
  count_prev = dat$count_prev,
  dt = dat$dt
)
prior.list <- list(
  prior_mu_b = 0.2,
  prior_mu_d = 0.1,
  prior_s_b = 0.1,
  prior_s_d = 0.1,
  prior_mu_obs_err = 0.05,
  prior_s_obs_err = 0.05
)
initfun <- function() {
  list(
    b = runif(1, 0, 0.2),
    d = runif(1, 0, 0.2),
    obs_err = runif(1, 0, 0.1)#,
    #count_err = runif(data.list$N, -0.5, 0.5)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 4,
  parallel_chains = 4,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 10
)

model_posterior$summary(variables = c("b", "d", "obs_err"))
params_posterior <- model_posterior$draws(variables = c("b", "d", "obs_err"), format = "df")
params_posterior %>% group_by(.chain) %>% summarize(b_mu = mean(b),
                                                    b_sd = sd(b),
                                                    d_mu = mean(d),
                                                    d_sd = sd(d),
                                                    obs_err_mu = mean(obs_err),
                                                    obs_err_sd = sd(obs_err))

pred <- model_posterior$draws(variables = c("pred"), format = "df")

pred_run <- as.numeric(pred %>% select(-c(.chain, .iteration, .draw)) %>% summarize_all(last))
dat$pred_run <- pred_run

dat %>% ggplot(aes(x = time, group = sample)) +
  geom_line(aes(y = count)) +
  geom_line(aes(y = type1_mean), color = "green") +
  geom_line(aes(y = approx), color = "red") +
  geom_line(aes(y = pred_run), color = "blue") +
  scale_y_log10()

dat %>% filter(time == 5) %>% ggplot() +
  geom_density(aes(x = count)) +
  geom_density(aes(x = type1_mean), color = "green") +
  geom_density(aes(x = approx), color = "red") +
  geom_density(aes(x = pred_run), color = "blue")
################################################################################
# Compare to ground truth
ground_truth <- data.frame(stat = c("b", "d", "obs_err"), true_value = c(br, dr, obs_err))

params_posterior_l <- params_posterior %>% select(.chain, b, d, obs_err) %>%
  pivot_longer(!c(".chain"), names_to="stat")

params_posterior_l %>%
  ggplot(aes(x = value, group = .chain)) +
  geom_density() +
  geom_vline(data = ground_truth, aes(xintercept = true_value), linewidth = 2) +
  facet_wrap(. ~ stat, scales = "free") +
  theme_bw()


params_posterior_l %>%
  left_join(ground_truth) %>%
  mutate(diff = value - true_value) %>%
  group_by(.chain, stat) %>%
  summarize(mean_diff = mean(diff)) %>%
  pivot_wider(names_from = stat, values_from = mean_diff)

params_posterior_l %>%
  group_by(.chain, stat) %>%
    summarize(mean_stat = mean(value)) %>%
    pivot_wider(names_from = stat, values_from = mean_stat)
