//  Logistic Birth-Death Model with birth and death rates as parameters
//  Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;               // time between observations
  array[N] int<lower=1> c_idx_1;      // index (int) of each concentration
  int<lower=1> nc1;            // number of unique doses
  vector[nc1] conc1;              // value of each concentration
  array[N] int<lower=1> c_idx_2;      // index (int) of each concentration
  int<lower=1> nc2;            // number of unique doses
  vector[nc2] conc2;              // value of each concentration
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
  
  real theta_1_mu_b50_1;              // theta_1 for prior distribution of mean of b50_1
  real<lower=0> theta_2_mu_b50_1;     // theta_2 for prior distribution of mean of b50_1
  real theta_1_s_b50_1;               // theta_1 for prior distribution of stdev of b50_1
  real<lower=0> theta_2_s_b50_1;      // theta_2 for prior distribution of stdev of b50_1
  
  real theta_1_mu_bh_1;              // theta_1 for prior distribution of mean of bh_1
  real<lower=0> theta_2_mu_bh_1;     // theta_2 for prior distribution of mean of bh_1
  real theta_1_s_bh_1;               // theta_1 for prior distribution of stdev of bh_1
  real<lower=0> theta_2_s_bh_1;      // theta_2 for prior distribution of stdev of bh_1

  real theta_1_mu_b50_2;              // theta_1 for prior distribution of mean of b50_2
  real<lower=0> theta_2_mu_b50_2;     // theta_2 for prior distribution of mean of b50_2
  real theta_1_s_b50_2;               // theta_1 for prior distribution of stdev of b50_2
  real<lower=0> theta_2_s_b50_2;      // theta_2 for prior distribution of stdev of b50_2
  
  real theta_1_mu_bh_2;              // theta_1 for prior distribution of mean of bh_2
  real<lower=0> theta_2_mu_bh_2;     // theta_2 for prior distribution of mean of bh_2
  real theta_1_s_bh_2;               // theta_1 for prior distribution of stdev of bh_2
  real<lower=0> theta_2_s_bh_2;      // theta_2 for prior distribution of stdev of bh_2
  
  real theta_1_mu_d0;              // theta_1 for prior distribution of mean of d0
  real<lower=0> theta_2_mu_d0;     // theta_2 for prior distribution of mean of d0
  real theta_1_s_d0;               // theta_1 for prior distribution of stdev of d0
  real<lower=0> theta_2_s_d0;      // theta_2 for prior distribution of stdev of d0
  
  real theta_1_mu_dd;              // theta_1 for prior distribution of mean of dd
  real<lower=0> theta_2_mu_dd;     // theta_2 for prior distribution of mean of dd
  real theta_1_s_dd;               // theta_1 for prior distribution of stdev of dd
  real<lower=0> theta_2_s_dd;      // theta_2 for prior distribution of stdev of dd
  
  real theta_1_mu_d50_1;              // theta_1 for prior distribution of mean of d50_1
  real<lower=0> theta_2_mu_d50_1;     // theta_2 for prior distribution of mean of d50_1
  real theta_1_s_d50_1;               // theta_1 for prior distribution of stdev of d50_1
  real<lower=0> theta_2_s_d50_1;      // theta_2 for prior distribution of stdev of d50_1
  
  real theta_1_mu_dh_1;              // theta_1 for prior distribution of mean of dh_1
  real<lower=0> theta_2_mu_dh_1;     // theta_2 for prior distribution of mean of dh_1
  real theta_1_s_dh_1;               // theta_1 for prior distribution of stdev of dh_1
  real<lower=0> theta_2_s_dh_1;      // theta_2 for prior distribution of stdev of dh_1

    real theta_1_mu_d50_2;              // theta_1 for prior distribution of mean of d50_2
  real<lower=0> theta_2_mu_d50_2;     // theta_2 for prior distribution of mean of d50_2
  real theta_1_s_d50_2;               // theta_1 for prior distribution of stdev of d50_2
  real<lower=0> theta_2_s_d50_2;      // theta_2 for prior distribution of stdev of d50_2
  
  real theta_1_mu_dh_2;              // theta_1 for prior distribution of mean of dh_2
  real<lower=0> theta_2_mu_dh_2;     // theta_2 for prior distribution of mean of dh_2
  real theta_1_s_dh_2;               // theta_1 for prior distribution of stdev of dh_2
  real<lower=0> theta_2_s_dh_2;      // theta_2 for prior distribution of stdev of dh_2
  
  real<lower=0> prior_mu_obs_err;  // mean parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error
}
parameters{
  // parameters
  vector<lower=0>[ng] bi;
  vector<lower=0>[ng] db;
  vector[ng] b50_1;
  vector<lower=0>[ng] bh_1;
  vector[ng] b50_2;
  vector<lower=0>[ng] bh_2;
  vector<lower=0>[ng] d0;
  vector<lower=0>[ng] dd;
  vector[ng] d50_1;
  vector<lower=0>[ng] dh_1;
  vector[ng] d50_2;
  vector<lower=0>[ng] dh_2;
  
  real<lower=0> obs_err;
  
  // hyperparameters
  real<lower=0> mu_bi;
  real<lower=0> s_bi;
  real<lower=0> mu_db;
  real<lower=0> s_db;
  real mu_b50_1;
  real<lower=0> s_b50_1;
  real mu_b50_2;
  real<lower=0> s_b50_2;
  real<lower=0> mu_bh_1;
  real<lower=0> s_bh_1;
  real<lower=0> mu_bh_2;
  real<lower=0> s_bh_2;
  real<lower=0> mu_d0;
  real<lower=0> s_d0;
  real<lower=0> mu_dd;
  real<lower=0> s_dd;
  real mu_d50_1;
  real<lower=0> s_d50_1;
  real<lower=0> mu_dh_1;
  real<lower=0> s_dh_1;
  real mu_d50_2;
  real<lower=0> s_d50_2;
  real<lower=0> mu_dh_2;
  real<lower=0> s_dh_2;
  // To account for rounding to nearest cell
  //vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  array[ng] matrix<lower=0>[nc1,nc2] b;
  array[ng] matrix<lower=0>[nc1,nc2] d;
  //vector[N] z;
  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;
  
  
  for (g in 1:ng)
  {
    for (c1 in 1:nc1)
    {
      for (c2 in 1:nc2)
      {
        if ((conc1[c1] == negative_infinity()) || (conc2[c2] == negative_infinity()))
        {
          b[g][c1,c2] = bi[g] + db[g];
          d[g][c1,c2] = d0[g];
        }
        else
        {
          b[g][c1,c2] = bi[g] + db[g] / ((1 + exp(bh_1[g] * (conc1[c1] - b50_1[g]))) * (1 + exp(bh_2[g] * (conc2[c2] - b50_2[g]))));
          d[g][c1,c2] = d0[g] + dd[g] - dd[g] / ((1 + exp(dh_1[g] * (conc1[c1] - d50_1[g]))) * (1 + exp(dh_2[g] * (conc2[c2] - d50_2[g]))));
        }
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
  mu_b50_1 ~ normal(theta_1_mu_b50_1, theta_2_mu_b50_1);
  s_b50_1 ~ normal(theta_1_s_b50_1, theta_2_s_b50_1);
  mu_bh_1 ~ normal(theta_1_mu_bh_1, theta_2_mu_bh_1);
  s_bh_1 ~ normal(theta_1_s_bh_1, theta_2_s_bh_1);
  mu_b50_2 ~ normal(theta_1_mu_b50_2, theta_2_mu_b50_2);
  s_b50_2 ~ normal(theta_1_s_b50_2, theta_2_s_b50_2);
  mu_bh_2 ~ normal(theta_1_mu_bh_2, theta_2_mu_bh_2);
  s_bh_2 ~ normal(theta_1_s_bh_2, theta_2_s_bh_2);
  mu_d0 ~ normal(theta_1_mu_d0, theta_2_mu_d0);
  s_d0 ~ normal(theta_1_s_d0, theta_2_s_d0);
  mu_dd ~ normal(theta_1_mu_dd, theta_2_mu_dd);
  s_dd ~ normal(theta_1_s_dd, theta_2_s_dd);
  mu_d50_1 ~ normal(theta_1_mu_d50_1, theta_2_mu_d50_1);
  s_d50_1 ~ normal(theta_1_s_d50_1, theta_2_s_d50_1);
  mu_dh_1 ~ normal(theta_1_mu_dh_1, theta_2_mu_dh_1);
  s_dh_1 ~ normal(theta_1_s_dh_1, theta_2_s_dh_1);
  mu_d50_2 ~ normal(theta_1_mu_d50_2, theta_2_mu_d50_2);
  s_d50_2 ~ normal(theta_1_s_d50_2, theta_2_s_d50_2);
  mu_dh_2 ~ normal(theta_1_mu_dh_2, theta_2_mu_dh_2);
  s_dh_2 ~ normal(theta_1_s_dh_2, theta_2_s_dh_2);
  // b and d are vectors but shouldn't matter
  bi ~ normal(mu_bi, s_bi);
  db ~ normal(mu_db, s_db);
  b50_1 ~ normal(mu_b50_1, s_b50_1);
  bh_1 ~ normal(mu_bh_1, s_bh_1);
  b50_2 ~ normal(mu_b50_2, s_b50_2);
  bh_2 ~ normal(mu_bh_2, s_bh_2);
  d0 ~ normal(mu_d0, s_d0);
  dd ~ normal(mu_dd, s_dd);
  d50_1 ~ normal(mu_d50_1, s_d50_1);
  dh_1 ~ normal(mu_dh_1, s_dh_1);
  d50_2 ~ normal(mu_d50_2, s_d50_2);
  dh_2 ~ normal(mu_dh_2, s_dh_2);
  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);

  // analytical values for mean and variance of a single type birth death process
  for ( n in 1:N ) {
    sigma[n] = sqrt(count_prev[n] * (b[g_idx[n]][c_idx_1[n],c_idx_2[n]]+d[g_idx[n]][c_idx_1[n],c_idx_2[n]]) / (b[g_idx[n]][c_idx_1[n],c_idx_2[n]]-d[g_idx[n]][c_idx_1[n],c_idx_2[n]]) *
      (exp(2*(b[g_idx[n]][c_idx_1[n],c_idx_2[n]]-d[g_idx[n]][c_idx_1[n],c_idx_2[n]])*dt[n]) - exp((b[g_idx[n]][c_idx_1[n],c_idx_2[n]]-d[g_idx[n]][c_idx_1[n],c_idx_2[n]])*dt[n])));
    mu[n] = count_prev[n] * exp((b[g_idx[n]][c_idx_1[n],c_idx_2[n]] - d[g_idx[n]][c_idx_1[n],c_idx_2[n]]) * dt[n]);
  }
  // likelihood function
  count ~ normal( mu , sqrt(sigma.^2.0 + (obs_err * mu).^2.0) );
}
