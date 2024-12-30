# Define Model with formula and priors on parameters
model <- process_model(bp_model(transition(name = "birth", parent = 1, offspring = c(2, 0), model = "birth ~ bi + db / (1 + exp(bh * (x - b50)))"), # birth
           transition(name = "death", parent = 1, offspring = c(0, 1), model = "death ~ d0 + dd  - dd / (1 + exp(dh * (x - d50)))"), # death
           transition(name = "k", parent = 2, offspring = c(0, 0))),
  parameter_constraints = c("bi[0,]", "db[0,]", "bh[0,]", "d0[0,]", "dd[0,]", "dh[0,]"),
  priors = c(
    "bi ~ normal(prior_mu_bi, prior_s_bi[0,])",
    "db ~ normal(prior_mu_db, prior_s_db[0,])",
    "b50 ~ normal(prior_mu_b50, prior_s_b50[0,])",
    "bh ~ normal(prior_mu_bh, prior_s_bh[0,])",
    "d0 ~ normal(prior_mu_d0, prior_s_d0[0,])",
    "dd ~ normal(prior_mu_dd, prior_s_dd[0,])",
    "d50 ~ normal(prior_mu_d50, prior_s_d50[0,])",
    "dh ~ normal(prior_mu_dh, prior_s_dh[0,])",
    "k ~ normal(prior_mu_k, prior_s_k[0,])"
  ),
  predictor_names = "x",
  observation_error = FALSE
)

model <- process_model(
  bp_model(transition(name = "birth", parent = 1, offspring = c(2)), # birth
           transition(name = "death", parent = 1, offspring = c(0))), # death
  priors = c(
    "birth ~ normal(prior_mu_b, prior_s_b[0,])",
    "death ~ normal(prior_mu_d, prior_s_d[0,])"
  ),
  observation_error = FALSE
)

model <- process_model(
  bp_model(transition(name = "birth", parent = 1, offspring = c(2)), # birth
           transition(name = "death", parent = 1, offspring = c(0))), # death
  priors = c(
    "birth ~ normal(mu_b, s_b[0,])",
    "death ~ normal(mu_d, s_d[0,])",
    "mu_b ~ normal(theta_1_mu_b, theta_2_mu_b)",
    "s_b ~ normal(theta_1_s_b, theta_2_s_b)",
    "mu_d ~ normal(theta_1_mu_d, theta_2_mu_d)",
    "s_d ~ normal(theta_1_s_d, theta_2_s_d)"
  ),
  observation_error = FALSE,
  hierarchical = c("birth", "death")
)

model <- process_model(
  bp_model(transition(name = "birth", parent = 1, offspring = c(2), model = "birth ~ bi + db / (1 + exp(bh * (x - b50)))"), # birth
           transition(name = "death", parent = 1, offspring = c(0), model = "death ~ ddeath")), # death
  priors = c(
    "bi ~ normal(prior_mu_bi, prior_s_bi[0,])",
    "db ~ normal(prior_mu_db, prior_s_db[0,])",
    "b50 ~ normal(prior_mu_b50, prior_s_b50[0,])",
    "bh ~ normal(prior_mu_bh, prior_s_bh[0,])",
    "ddeath ~ normal(d1, d2)",
    "d1 ~ normal(prior_mu_d1, prior_s_d1)",
    "d2 ~ normal(prior_mu_d2, prior_s_d2)"
  ),
  observation_error = FALSE,
  predictor_names = "x",
  hierarchical = "ddeath"
)


# add additional option to change x to user defined (but also specified)

# TEST
mystancode <- create_stan_code(model)
# Next set up data - note: may need to change x in the model formula based on data predictor

f2 <- write_stan_file(mystancode)
testmodel <- cmdstanr::cmdstan_model(f2, force_recompile = T)
dat$log_conc[dat$log_conc == -Inf] <- -12
data.list <- attach_data(model, dat,
                         count_vars = "type1",
                         prevcount_vars = "type1_prev",
                         predictor_vars = "log_conc")

data.list <- attach_data(model, dat,
  count_vars = c("type1", "type2"),
  prevcount_vars = c("type1_prev", "type2_prev"),
  predictor_vars = "log_conc"
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
  prior_s_dh = 0.5,
  prior_mu_k = 0.05,
  prior_s_k = 0.03
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
    k = runif(1, 0, 0.05),
    count_err = matrix(runif(data.list$N * data.list$ntypes, -0.5, 0.5), ncol = data.list$ntypes)
  )
}

model_posterior <- testmodel$sample(
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)

params_posterior <- model_posterior$draws(variables = c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh", "k"), format = "df")
head(params_posterior)

############################################################

### Mixed effects specification (Work from here)

model <- process_model(
  bp_model(
    transition(name = "birth", parent = 1, offspring = c(2, 0), model = "birth ~ bi + db / (1 + exp(bh * (x - b50)))"),
    transition(name = "death", parent = 1, offspring = c(0, 1), model = "death ~ d0 + dd  - dd / (1 + exp(dh * (x - d50)))"),
    transition(name = "clearance", parent = 2, offspring = c(0, 0), model = "clearance ~ k")
  ),
  parameter_constraints = c("bi[0,]", "db[0,]", "bh[0,]", "d0[0,]", "dd[0,]", "dh[0,]", "k[0,]"),
  priors = c(
    "bi ~ normal(mu_bi, s_bi[0,])",
    "db ~ normal(mu_db, s_db[0,])",
    "b50 ~ normal(mu_b50, s_b50[0,])",
    "bh ~ normal(mu_bh, s_bh[0,])",
    "d0 ~ normal(mu_d0, s_d0[0,])",
    "dd ~ normal(mu_dd, s_dd[0,])",
    "d50 ~ normal(mu_d50, s_d50[0,])",
    "dh ~ normal(mu_dh, s_dh[0,])",
    "k ~ normal(mu_k, s_k[0,])",

    "mu_bi ~ normal(theta_1_mu_bi, theta_2_mu_bi)",
    "s_bi[0,] ~ normal(theta_1_s_bi, theta_2_s_bi)",
    "mu_db ~ normal(theta_1_mu_db, theta_2_mu_db)",
    "s_db[0,] ~ normal(theta_1_s_db, theta_2_s_db)",
    "mu_b50 ~ normal(theta_1_mu_b50, theta_2_mu_b50)",
    "s_b50[0,] ~ normal(theta_1_s_b50, theta_2_s_b50)",
    "mu_bh ~ normal(theta_1_mu_bh, theta_2_mu_bh)",
    "s_bh[0,] ~ normal(theta_1_s_bh, theta_2_s_bh)",
    "mu_d0 ~ normal(theta_1_mu_d0, theta_2_mu_d0)",
    "s_d0[0,] ~ normal(theta_1_s_d0, theta_2_s_d0)",
    "mu_dd ~ normal(theta_1_mu_dd, theta_2_mu_dd)",
    "s_dd[0,] ~ normal(theta_1_s_dd, theta_2_s_dd)",
    "mu_d50 ~ normal(theta_1_mu_d50, theta_2_mu_d50)",
    "s_d50[0,] ~ normal(theta_1_s_d50, theta_2_s_d50)",
    "mu_dh ~ normal(theta_1_mu_dh, theta_2_mu_dh)",
    "s_dh[0,] ~ normal(theta_1_s_dh, theta_2_s_dh)",
    "mu_k ~ normal(theta_1_mu_k, theta_2_mu_k)",
    "s_k[0,] ~ normal(theta_1_s_k, theta_2_s_k)"
  ),
  hierarchical = c("bi", "db", "b50", "bh", "d0", "dd", "d50", "dh", "k"),
  predictor_names = "x",
  observation_error = FALSE
)


mystancode <- create_stan_code(model)
f2 <- write_stan_file(mystancode)
testmodel <- cmdstanr::cmdstan_model(f2)
data.list <- attach_data(model, dat,
  count_vars = c("type1", "type2"),
  prevcount_vars = c("type1_prev", "type2_prev"),
  predictor_vars = "log_conc",
  hierarchical_group = "group"
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
    bi = runif(data.list$ng, 0.01, 1),
    db = runif(data.list$ng, 0.1, 1),
    b50 = runif(data.list$ng, -2, 1),
    bh = runif(data.list$ng, 0.1, 1),
    d0 = runif(data.list$ng, 0.01, 1),
    dd = runif(data.list$ng, 0.1, 1),
    d50 = runif(data.list$ng, -2, 1),
    dh = runif(data.list$ng, 0.1, 1),
    k = runif(data.list$ng, 0.01, 0.1),
    count_err = matrix(runif(data.list$N * data.list$ntypes, -0.5, 0.5), ncol = data.list$ntypes)
  )
}

model_posterior <- testmodel$sample(
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 1,
  init = initfun,
  iter_warmup = 1000,
  iter_sampling = 1000
)
