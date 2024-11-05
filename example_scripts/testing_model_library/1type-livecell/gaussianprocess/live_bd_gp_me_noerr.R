### Gaussian Process Concentration of Birth and Death Rates
### No Error
### Mixed Effects

########################################################################
# Set Project Directory Location---------------------------------
libloc <- Sys.getenv("R_LIBS_USER")
if (grepl("michorlab", libloc)) {
  setwd("/michorlab/mcdonald/bp-clt")
} else {
  setwd("~/Dropbox/michor/projects/bp-clt")
}

library(tidyverse)
library(estipop)
library(cmdstanr)

# Generate data -------------------------------------------------
# Generate Parameters: make b50 = d50 and bh = dh create proper gr curve
set.seed(75486)
b0 <- rnorm(7, 1, 0.02)
bi <- rnorm(7, 0.2, 0.03)
db <- b0 - bi
b50 <- rnorm(7, 0.2, 0.04)
logb50 <- log(b50) # convert to logscale
bh <- rnorm(7, 2, 0.1)
d0 <- rnorm(7, 0.3, 0.05)
di <- rnorm(7, 0.5, 0.02)
dd <- di - d0
d50 <- rnorm(7, 1, 0.05)
logd50 <- log(d50) # convert to log scale
dh <- bh


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
stanfileloc <- system.file("model_library/1type-livecell/gaussianprocess/birthdeath_mixedeffects.stan", package = "bestdr")
stanmodel <- cmdstan_model(stanfileloc)


# Set up data for inputs
dat$log_conc_idx <- as.numeric(factor(dat$log_conc, ordered = T))
log_conc <- unique(dat$log_conc)
dim(log_conc) <- c(length(log_conc), 1)
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
  nd = 1,
  g_idx = as.numeric(as.factor(dat$group)),
  ng = num_groups
)
prior.list <- list(
  theta_1_alpha_rho = 2,
  theta_2_alpha_rho = 1,
  theta_1_beta_rho = 2,
  theta_2_beta_rho = 1,
  
  theta_1_mu_alpha_b = 0,
  theta_2_mu_alpha_b = 0.2,
  theta_1_s_alpha_b = 0.2,
  theta_2_s_alpha_b = 0.1,
  theta_1_mu_alpha_d = 0,
  theta_2_mu_alpha_d = 0.2,
  theta_1_s_alpha_d = 0.2,
  theta_2_s_alpha_d = 0.1
)
initfun <- function() {
  list(
    count_err = runif(data.list$N, -0.5, 0.5)
  )
}

model_posterior <- stanmodel$sample(
  data = c(data.list, prior.list),
  chains = 1,
  parallel_chains = 1,
  refresh = 1,
  init = initfun,
  iter_warmup = 10,
  iter_sampling = 10
)
