//  Logistic Birth-Death Model with birth and death rates as parameters
//  Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;               // time between observations
  array[N] int<lower=1> c_idx;      // index (int) of each concentration
  int<lower=1> nc;            // number of unique doses
  array[nc] real conc;              // value of each concentration
  array[N] int<lower=1> g_idx;      // group id (random effect)
  int<lower=1> ng;            // number of groups
  
  // Parameters for priors of hyperparameters
  real theta_1_mu_bi;              // theta_1 for prior distribution of mean of bi
  real<lower=0> theta_2_mu_bi;     // theta_2 for prior distribution of mean of bi
  real theta_1_s_bi;               // theta_1 for prior distribution of stdev of bi
  real<lower=0> theta_2_s_bi;      // theta_2 for prior distribution of stdev of bi
  
  real theta_1_mu_db;              // theta_1 for prior distribution of mean of db
  real<lower=0> theta_2_mu_db;     // theta_2 for prior distribution of mean of db
  real theta_1_s_db;               // theta_1 for prior distribution of stdev of db
  real<lower=0> theta_2_s_db;      // theta_2 for prior distribution of stdev of db
  
  real theta_1_mu_b50;              // theta_1 for prior distribution of mean of b50
  real<lower=0> theta_2_mu_b50;     // theta_2 for prior distribution of mean of b50
  real theta_1_s_b50;               // theta_1 for prior distribution of stdev of b50
  real<lower=0> theta_2_s_b50;      // theta_2 for prior distribution of stdev of b50
  
  real theta_1_mu_bh;              // theta_1 for prior distribution of mean of bh
  real<lower=0> theta_2_mu_bh;     // theta_2 for prior distribution of mean of bh
  real theta_1_s_bh;               // theta_1 for prior distribution of stdev of bh
  real<lower=0> theta_2_s_bh;      // theta_2 for prior distribution of stdev of bh
  
  real theta_1_mu_d0;              // theta_1 for prior distribution of mean of d0
  real<lower=0> theta_2_mu_d0;     // theta_2 for prior distribution of mean of d0
  real theta_1_s_d0;               // theta_1 for prior distribution of stdev of d0
  real<lower=0> theta_2_s_d0;      // theta_2 for prior distribution of stdev of d0
  
  real theta_1_mu_dd;              // theta_1 for prior distribution of mean of dd
  real<lower=0> theta_2_mu_dd;     // theta_2 for prior distribution of mean of dd
  real theta_1_s_dd;               // theta_1 for prior distribution of stdev of dd
  real<lower=0> theta_2_s_dd;      // theta_2 for prior distribution of stdev of dd
  
  real theta_1_mu_d50;              // theta_1 for prior distribution of mean of d50
  real<lower=0> theta_2_mu_d50;     // theta_2 for prior distribution of mean of d50
  real theta_1_s_d50;               // theta_1 for prior distribution of stdev of d50
  real<lower=0> theta_2_s_d50;      // theta_2 for prior distribution of stdev of d50
  
  real theta_1_mu_dh;              // theta_1 for prior distribution of mean of dh
  real<lower=0> theta_2_mu_dh;     // theta_2 for prior distribution of mean of dh
  real theta_1_s_dh;               // theta_1 for prior distribution of stdev of dh
  real<lower=0> theta_2_s_dh;      // theta_2 for prior distribution of stdev of dh
  
  real<lower=0> prior_mu_obs_err;  // mean parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error
}
parameters{
  // parameters
  vector<lower=0>[ng] bi;
  vector<lower=0>[ng] db;
  vector[ng] b50;
  vector<lower=0>[ng] bh;
  vector<lower=0>[ng] d0;
  vector<lower=0>[ng] dd;
  vector[ng] d50;
  vector<lower=0>[ng] dh;
  
  real<lower=0> obs_err;
  
  // hyperparameters
  real<lower=0> mu_bi;
  real<lower=0> s_bi;
  real<lower=0> mu_db;
  real<lower=0> s_db;
  real mu_b50;
  real<lower=0> s_b50;
  real<lower=0> mu_bh;
  real<lower=0> s_bh;
  real<lower=0> mu_d0;
  real<lower=0> s_d0;
  real<lower=0> mu_dd;
  real<lower=0> s_dd;
  real mu_d50;
  real<lower=0> s_d50;
  real<lower=0> mu_dh;
  real<lower=0> s_dh;
  // To account for rounding to nearest cell
  //vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  array[ng,nc] real<lower=0> b;
  array[ng,nc] real<lower=0> d;
  //vector[N] z;
  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;
  
  
  for (g in 1:ng)
  {
    for (c in 1:nc)
    {
      if (conc[c] == negative_infinity())
      {
        b[g,c] = bi[g] + db[g];
        d[g,c] = d0[g];
      }
      else
      {
        b[g,c] = bi[g] + db[g] / (1 + exp(bh[g] * (conc[c] - b50[g])));
        d[g,c] = d0[g] + dd[g] -  dd[g] / (1 + exp(dh[g] * (conc[c] - d50[g])));
      }

    }
  }

  
}
model{
  vector[N] sigma;
  vector[N] mu;

  // specified prior distributions
  mu_bi ~ normal(theta_1_mu_bi, theta_2_mu_bi);
  s_bi ~ normal(theta_1_s_bi, theta_2_s_bi);
  mu_db ~ normal(theta_1_mu_db, theta_2_mu_db);
  s_db ~ normal(theta_1_s_db, theta_2_s_db);
  mu_b50 ~ normal(theta_1_mu_b50, theta_2_mu_b50);
  s_b50 ~ normal(theta_1_s_b50, theta_2_s_b50);
  mu_bh ~ normal(theta_1_mu_bh, theta_2_mu_bh);
  s_bh ~ normal(theta_1_s_bh, theta_2_s_bh);
  mu_d0 ~ normal(theta_1_mu_d0, theta_2_mu_d0);
  s_d0 ~ normal(theta_1_s_d0, theta_2_s_d0);
  mu_dd ~ normal(theta_1_mu_dd, theta_2_mu_dd);
  s_dd ~ normal(theta_1_s_dd, theta_2_s_dd);
  mu_d50 ~ normal(theta_1_mu_d50, theta_2_mu_d50);
  s_d50 ~ normal(theta_1_s_d50, theta_2_s_d50);
  mu_dh ~ normal(theta_1_mu_dh, theta_2_mu_dh);
  s_dh ~ normal(theta_1_s_dh, theta_2_s_dh);
  
  // b and d are vectors but shouldn't matter
  bi ~ normal(mu_bi, s_bi);
  db ~ normal(mu_db, s_db);
  b50 ~ normal(mu_b50, s_b50);
  bh ~ normal(mu_bh, s_bh);
  d0 ~ normal(mu_d0, s_d0);
  dd ~ normal(mu_dd, s_dd);
  d50 ~ normal(mu_d50, s_d50);
  dh ~ normal(mu_dh, s_dh);
  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);

  // analytical values for mean and variance of a single type birth death process
  for ( n in 1:N ) {
    sigma[n] = sqrt(count_prev[n] * (b[g_idx[n],c_idx[n]]+d[g_idx[n],c_idx[n]]) / (b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]]) *
      (exp(2*(b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]])*dt[n]) - exp((b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]])*dt[n])));
    mu[n] = count_prev[n] * exp((b[g_idx[n],c_idx[n]] - d[g_idx[n],c_idx[n]]) * dt[n]);
  }
  // likelihood function
  count ~ normal( mu , sqrt(sigma.^2.0 + (obs_err * mu).^2.0) );
}
/*
generated quantities{
  vector[N] log_lik;
  
  {
    real sigma;
    real mu;
    for (n in 1:N)
    {
      sigma = sqrt(count_prev[n] * (b[g_idx[n],c_idx[n]]+d[g_idx[n],c_idx[n]]) / (b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]]) *
        (exp(2*(b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]])*dt[n]) - exp((b[g_idx[n],c_idx[n]]-d[g_idx[n],c_idx[n]])*dt[n])));
      mu = count_prev[n] * exp((b[g_idx[n],c_idx[n]] - d[g_idx[n],c_idx[n]]) * dt[n]);
      log_lik[n] = normal_lpdf(count[n]| mu, sigma);
    }
  }
  
}
*/
