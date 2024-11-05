### Logistic Concentration of Birth, Death, Clearance Rates
### No Error
### No Mixed Effects

library(magrittr)
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


  N <- 10 # Number of samples
  tot_time <- 3
  dt <- 0.1 # Time between measurements
  times <- seq(0, tot_time, dt)
  init_pop <- c(1000, 0) # Number of ancestors


  # Simulate process with estipop
  set.seed(90210)
  dat <- data.frame()
  for (i in seq_along(concentrations)) {
    params <- c(b_truth[i], d_truth[i], k)
    # time-homogeneous two-type model
    model <- estipop::process_model(
      estipop::transition(rate = estipop::rate(params[1]), parent = 1, offspring = c(2, 0)), # birth
      estipop::transition(rate = estipop::rate(params[2]), parent = 1, offspring = c(0, 1)), # death
      estipop::transition(rate = estipop::rate(params[3]), parent = 2, offspring = c(0, 0)) # clearance
    )

    dat_temp <- estipop::branch(model, params, init_pop, times, N, silent = T)
    dat_temp$concentration <- concentrations[i]
    dat_temp$log_conc <- log_conc[i]
    dat <- dplyr::bind_rows(dat, dat_temp)
  }

  dat <- dat %>%
    dplyr::rename(sample = rep) %>%
    dplyr::group_by(sample, concentration) %>%
    dplyr:: mutate(
      dt = signif(time - dplyr::lag(time), 6),
      type1_prev = lag(type1),
      type2_prev = lag(type2)
    )


  dat <- dat %>% dplyr::filter(!is.na(dt))
