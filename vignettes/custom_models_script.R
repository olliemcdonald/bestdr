# Model 1 - Simple Birth Death Process -----------------------------------------

b <- 0.5
d <- 0.1

N <- 10 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(100) # Number of ancestors

set.seed(8437)
params <- c(b, d)
model <- estipop::process_model(
  estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = c(2)), # birth
  estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = c(1)) # death
)
dat <- estipop::branch(model, params, init_pop, times, N, silent = T)

dat %>% ggplot(aes(x = time, y = type1, group = rep)) +
  geom_line() +
  scale_y_log10()


model <- process_model(
  bp_model(transition(name = "birth", parent = 1, offspring = c(2)), # birth
           transition(name = "death", parent = 1, offspring = c(0))), # death
  priors = c(
    "birth ~ normal(mu_b, s_b[0,])",
    "death ~ normal(mu_d, s_d[0,])"
  ),
  parameter_constraints = c("birth[0,]", "death[0,]"),
  observation_error = FALSE
)

mystancode <- create_stan_code(model)
f2 <- write_stan_file(mystancode)
stanmodel <- cmdstanr::cmdstan_model(f2, force_recompile = T)

dat_dt <- dat %>% group_by(rep) %>%
  mutate(dt = time - lag(time),
         type1_prev = lag(type1)) %>%
  filter(!is.na(dt))
dat_dt <- as.data.frame(dat_dt)

data.list <- attach_data(model, dat_dt,
                         time_var = "dt",
                         count_vars = "type1",
                         prevcount_vars = "type1_prev")

prior.list <- list(
  mu_b = 1,
  s_b = 0.3,
  mu_d = 1,
  s_d = 0.3
)
initfun <- function() {
  list(
    birth = runif(1, 0, 2),
    death = runif(1, 0, 2),
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 0,
  init = initfun,
  iter_warmup = 500,
  iter_sampling = 500
)

model_posterior$summary()
params_posterior <- model_posterior$draws(variables = c("birth", "death"), format = "df")

true_data <- data.frame(stat = c("birth", "death"), value = c(b, d))
params_posterior %>%
  pivot_longer(c("birth", "death"), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = value, fill = stat)) +
  geom_density() +
  geom_vline(data = true_data, aes(xintercept = value))

# Model 2 - 2-type Birth-Death-Mutation ----------------------------------------
# Type 1 parameters
##  Birth - Constant
##  Death - 4-parameter logistic
##  Mutation - Constant
# Type 2 parameters
##  Birth - Constant
##  Death - 4-parameter logistic


## Generate data -------------------------------------------------
# Generate Parameters:
# Generate curve from 4-parameter logistic function for birth and death
b <- 0.10

d0_1 <- 0.005
di_1 <- 0.12
d50_1 <- -2
dh_1 <- 1
dd_1 <- di_1 - d0_1

d0_2 <- 0.002
di_2 <- 0.10
d50_2 <- 0
dh_2 <- 0.5
dd_2 <- di_2 - d0_2

u_12 <- 0.005

param_truth <- data.frame(
  birth = b,
  d0_1 = d0_1, d50_1 = d50_1, dh_1 = dh_1, dd_1 = dd_1,
  d0_2 = d0_2, d50_2 = d50_2, dh_2 = dh_2, dd_2 = dd_2,
  u_12 = u_12
)
print(param_truth)

concentrations <- c(0, 10^seq(-4, 4, 0.1))
log_conc <- log(concentrations)
d_curve <- function(lc, d0, dd, dh, d50) d0 + dd - dd / (1 + exp(dh * (lc - d50)))

d1_truth <- d_curve(log_conc, d0_1, dd_1, dh_1, d50_1)
d2_truth <- d_curve(log_conc, d0_2, dd_2, dh_2, d50_2)

data.frame(concentrations, b = b, d1 = d1_truth, d2 = d2_truth, u_12 = u_12) %>%
  pivot_longer(-concentrations, names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = concentrations, y = value, color = stat)) +
  geom_line() +
  scale_x_log10()

concentrations <- c(0, 10^seq(-4, 4, 1))
log_conc <- log(concentrations)
d1_truth <- d_curve(log_conc, d0_1, dd_1, dh_1, d50_1)
d2_truth <- d_curve(log_conc, d0_2, dd_2, dh_2, d50_2)


N <- 3 # Number of samples
tot_time <- 3
dt <- 0.1 # Time between measurements
times <- seq(0, tot_time, dt)
init_pop <- c(500, 100) # Number of ancestors

# Simulate process with estipop
set.seed(29837)
dat <- data.frame()
for (i in seq_along(concentrations)) {
  params <- c(b, d1_truth[i], u_12, b, d2_truth[i])
  # time-homogeneous two-type model
  model <- estipop::process_model(
    estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = c(2, 0)), # 1 birth
    estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = c(0, 0)), # 1 death
    estipop::transition(rate = estipop::rate(params[3]), parent = 1, offspring = c(0, 1)), # 1->2 mutation
    estipop::transition(rate = estipop::rate(params[4]), parent = 2, offspring = c(0, 2)), # 1 birth
    estipop::transition(rate = estipop::rate(params[5]), parent = 2, offspring = c(0, 0))  # 1 death
  )
  
  dat_temp <- estipop::branch(model, params, init_pop, times, N, silent = T)
  dat_temp$concentration <- concentrations[i]
  dat <- dplyr::bind_rows(dat, dat_temp)
}

dat %>% pivot_longer(c(type1, type2), names_to = "type", values_to = "count") %>%
  ggplot(aes(x = time, y = count, color = factor(concentration), group = interaction(rep, factor(concentration)))) +
  geom_line() +
  scale_color_manual(values = colorRampPalette(c("gold", "darkblue"))( 10 )) +
  facet_wrap(type ~ ., scales = "free_y")

dat_dr <- dat %>%
  dplyr::group_by(rep, concentration) %>%
  dplyr:: mutate(
    dt = signif(time - dplyr::lag(time), 6),
    type1_prev = lag(type1),
    type2_prev = lag(type2)
  ) %>% dplyr::filter(!is.na(dt))


## Define model -------------------------------------------------
model <- process_model(
  bp_model(transition(name = "birth_1", parent = 1, offspring = c(2, 0)),
           transition(name = "death_1", parent = 1, offspring = c(0, 0),
                      model = "death_1 ~ d0_1 + dd_1  - dd_1 / (1 + exp(dh_1 * (x - d50_1)))"),
           transition(name = "mut_12", parent = 1, offspring = c(0, 1)),
           transition(name = "birth_2", parent = 1, offspring = c(0, 2)),
           transition(name = "death_2", parent = 1, offspring = c(0, 0),
                      model = "death_2 ~ d0_2 + dd_2  - dd_2 / (1 + exp(dh_2 * (x - d50_2)))")),
  parameter_constraints = c("birth_1[0,]", "birth_2[0,]", "mut_12[0,]",
                            "d0_1[0,]", "dd_1[0,]", "dh_1[0,]",
                            "d0_2[0,]", "dd_2[0,]", "dh_2[0,]"),
  priors = c(
    "birth_1 ~ normal(mu_b1, s_b1[0,])",
    "birth_2 ~ normal(mu_b2, s_b2[0,])",
    "mut_12 ~ normal(mu_u, s_u[0,])",
    "d0_1 ~ normal(mu_d0_1, s_d0_1[0,])",
    "dd_1 ~ normal(mu_dd_1, s_dd_1[0,])",
    "d50_1 ~ normal(mu_d50_1, s_d50_1[0,])",
    "dh_1 ~ normal(mu_dh_1, s_dh_1[0,])",
    "d0_2 ~ normal(mu_d0_2, s_d0_2[0,])",
    "dd_2 ~ normal(mu_dd_2, s_dd_2[0,])",
    "d50_2 ~ normal(mu_d50_2, s_d50_2[0,])",
    "dh_2 ~ normal(mu_dh_2, s_dh_2[0,])"
  ),
  predictor_names = "x",
  observation_error = FALSE
)
mystancode <- create_stan_code(model)
f2 <- write_stan_file(mystancode)
stanmodel <- cmdstanr::cmdstan_model(f2, force_recompile = T)

# Adding the data
data.list <- attach_data(model, dat_dr,
                         time_var = "dt",
                         count_vars = c("type1", "type2"),
                         prevcount_vars = c("type1_prev", "type2_prev"))

prior.list <- list(
  mu_b1 = 0.1,
  s_b1 = 0.05,
  mu_b2 = 0.1,
  s_b2 = 0.05,
  mu_u = 0.01,
  s_u = 0.02,
  mu_d0_1 = 0.01,
  s_d0_1 = 0.1,
  mu_dd_1 = 0.1,
  s_dd_1 = 0.1,
  mu_d50_1 = 0,
  s_d50_1 = 1,
  mu_dh_1 = 1,
  s_dh_1 = 1,
  mu_d0_2 = 0.01,
  s_d0_2 = 0.1,
  mu_dd_2 = 0.1,
  s_dd_2 = 0.1,
  mu_d50_2 = 0,
  s_d50_2 = 1,
  mu_dh_2 = 1,
  s_dh_2 = 1
)

model_posterior <- stanmodel$optimize(
  data = c(data.list, prior.list)
  #init = initfun,
)

model_posterior$summary()
param_truth

# Model 3 - 2-type Birth-Death-Mutation Hierarchical ---------------------------
# Type 1 parameters
##  Birth - Constant
##  Death - 4-parameter logistic
##  Mutation - Constant

model <- process_model(
  bp_model(transition(name = "birth_1", parent = 1, offspring = c(2, 0)),
           transition(name = "death_1", parent = 1, offspring = c(0, 0),
                      model = "death_1 ~ d0_1 + dd_1  - dd_1 / (1 + exp(dh_1 * (x - d50_1)))"),
           transition(name = "mut_12", parent = 1, offspring = c(0, 1)),
           transition(name = "birth_2", parent = 1, offspring = c(0, 2)),
           transition(name = "death_2", parent = 1, offspring = c(0, 0),
                      model = "death_2 ~ d0_2 + dd_2  - dd_2 / (1 + exp(dh_2 * (x - d50_2)))")),
  parameter_constraints = c("birth_1[0,]", "birth_2[0,]", "mut_12[0,]",
                            "d0_1[0,]", "dd_1[0,]", "dh_1[0,]",
                            "d0_2[0,]", "dd_2[0,]", "dh_2[0,]"),
  priors = c(
    "birth_1 ~ normal(mu_b1, s_b1[0,])",
    "birth_2 ~ normal(mu_b2, s_b2[0,])",
    "mut_12 ~ normal(mu_u, s_u[0,])",
    "d0_1 ~ normal(mu_d0_1, s_d0_1[0,])",
    "dd_1 ~ normal(mu_dd_1, s_dd_1[0,])",
    "d50_1 ~ normal(mu_d50_1, s_d50_1[0,])",
    "dh_1 ~ normal(mu_dh_1, s_dh_1[0,])",
    "d0_2 ~ normal(mu_d0_2, s_d0_2[0,])",
    "dd_2 ~ normal(mu_dd_2, s_dd_2[0,])",
    "d50_2 ~ normal(mu_d50_2, s_d50_2[0,])",
    "dh_2 ~ normal(mu_dh_2, s_dh_2[0,])",
    
    "mu_b1 ~ normal(theta_1_mu_birth1, theta_2_mu_birth1)",
    "s_b1 ~ normal(theta_1_s_birth1, theta_2_s_birth1)",
    "mu_d0_1 ~ normal(theta_1_mu_d0_1, theta_2_mu_d0_1)",
    "s_d0_1 ~ normal(theta_1_s_d0_1, theta_2_s_d0_1)",
    "mu_dd_1 ~ normal(theta_1_mu_dd_1, theta_2_mu_dd_1)",
    "s_dd_1 ~ normal(theta_1_s_dd_1, theta_2_s_dd_1)",
    "mu_d50_1 ~ normal(theta_1_mu_d50_1, theta_2_mu_d50_1)",
    "s_d50_1 ~ normal(theta_1_s_d50_1, theta_2_s_d50_1)",
    "mu_dh_1 ~ normal(theta_1_mu_dh_1, theta_2_mu_dh_1)",
    "s_dh_1 ~ normal(theta_1_s_dh_1, theta_2_s_dh_1)"
  ),
  hierarchical = c("birth_1", "d0_1", "dd_1", "d50_1", "dh_1"),
  predictor_names = "x",
  observation_error = FALSE
)
mystancode <- create_stan_code(model)
f2 <- write_stan_file(mystancode)
stanmodel <- cmdstanr::cmdstan_model(f2, force_recompile = T)
