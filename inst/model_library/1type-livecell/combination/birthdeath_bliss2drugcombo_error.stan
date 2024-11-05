// Logistic combination Birth-Death Model with birth and death rates as parameters
// No Mixed Effects
// Attempt with errors added in that is constant across concentrations
data{
  int<lower=1> N;             // number of samples
  vector<lower=0>[N] count;            // total count at current time
  vector<lower=0>[N] count_prev;       // total count at previous time
  vector<lower=0>[N] dt;               // time between observations
  array[N] int<lower=1> c_idx_1;      // index (int) of each concentration
  int<lower=1> nc1;            // number of unique doses
  vector[nc1] conc1;              // value of each concentration
  array[N] int<lower=1> c_idx_2;      // index (int) of each concentration
  int<lower=1> nc2;            // number of unique doses
  vector[nc2] conc2;              // value of each concentration
  
  real prior_mu_bi;            // mean for prior distribution of bi (right asymptote)
  real<lower=0> prior_s_bi;    // standard deviation for prior distribution of bi
  real prior_mu_db;            // diff of asymptotes
  real<lower=0> prior_s_db;    
  real prior_mu_b50;           // concentration at midpoint / inflection point
  real<lower=0> prior_s_b50;   
  real prior_mu_bh;            // slope at inflection point
  real<lower=0> prior_s_bh;    
  
  real prior_mu_d0;            // mean for prior distribution of d0 (left asymptote)
  real<lower=0> prior_s_d0;    // standard deviation for prior distribution of bi
  real prior_mu_dd;            // diff of asymptotes
  real<lower=0> prior_s_dd;    
  real prior_mu_d50;           // concentration at midpoint / inflection point
  real<lower=0> prior_s_d50;   
  real prior_mu_dh;            // slope at inflection point
  real<lower=0> prior_s_dh;

  real<lower=0> prior_mu_obs_err;  // mean parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error
   
}
parameters{
  //  using log birth rate and log death rate as priors since has better convergence issues
  real<lower=0> bi;
  real<lower=0> db;
  real b50_1;
  real b50_2;
  real<lower=0> bh_1;
  real<lower=0> bh_2;
  real<lower=0> d0;
  real<lower=0> dd;
  real d50_1;
  real d50_2;
  real<lower=0> dh_1;
  real<lower=0> dh_2;
  real<lower=0> obs_err;
  // To account for rounding to nearest cell
  //vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  matrix<lower=0>[nc1,nc2] b;
  matrix<lower=0>[nc1,nc2] d;
  // vector[N] z;
  // // True value is sum of rounded value (data) and roundoff error
  // //z = count + count_err;
  
  for (c1 in 1:nc1)
  {
    for(c2 in 1:nc2){
      if ((conc1[c1] == negative_infinity()) || (conc2[c2] == negative_infinity()))
      {
        b[c1,c2] = bi + db;
        d[c1,c2] = d0;
      }
      else
      {
        b[c1,c2] = bi + db / ((1 + exp(bh_1 * (conc1[c1] - b50_1))) * (1 + exp(bh_2 * (conc2[c2] - b50_2))));
        d[c1,c2] = d0 + dd - dd / ((1 + exp(dh_1 * (conc1[c1] - d50_1))) * (1 + exp(dh_2 * (conc2[c2] - d50_2))));
      }
    }
  }
}
model{
  vector[N] sigma;
  vector[N] mu;

  // specified prior distributions
  bi ~ normal(prior_mu_bi, prior_s_bi);
  db ~ normal(prior_mu_db, prior_s_db);
  b50_1 ~ normal(prior_mu_b50, prior_s_b50);
  b50_2 ~ normal(prior_mu_b50, prior_s_b50);
  bh_1 ~ normal(prior_mu_bh, prior_s_bh);
  bh_2 ~ normal(prior_mu_bh, prior_s_bh);
  d0 ~ normal(prior_mu_d0, prior_s_d0);
  dd ~ normal(prior_mu_dd, prior_s_dd);
  d50_1 ~ normal(prior_mu_d50, prior_s_d50);
  d50_2 ~ normal(prior_mu_d50, prior_s_d50);
  dh_1 ~ normal(prior_mu_dh, prior_s_dh);
  dh_2 ~ normal(prior_mu_dh, prior_s_dh);
  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);

  // analytical values for mean and variance of a single type birth death process
  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b[c_idx_1[i], c_idx_2[i]]+d[c_idx_1[i], c_idx_2[i]]) / (b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]]) * (exp(2*(b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]])*dt[i]) - exp((b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]])*dt[i])));
    mu[i] = count_prev[i] * exp((b[c_idx_1[i], c_idx_2[i]] - d[c_idx_1[i], c_idx_2[i]]) * dt[i]);
  }
  // likelihood function
  count ~ normal( mu , sqrt(sigma.^2.0 + (obs_err * mu).^2.0) );
}
// generated quantities{
//   vector[N] log_lik;
//   
//   {
//     real sigma;
//     real mu;
//     for (i in 1:N)
//     {
//       sigma = sqrt(count_prev[i] * (b[c_idx_1[i], c_idx_2[i]]+d[c_idx_1[i], c_idx_2[i]]) / (b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]]) * (exp(2*(b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]])*dt[i]) - exp((b[c_idx_1[i], c_idx_2[i]]-d[c_idx_1[i], c_idx_2[i]])*dt[i])));
//       mu = count_prev[i] * exp((b[c_idx_1[i], c_idx_2[i]] - d[c_idx_1[i], c_idx_2[i]]) * dt[i]);
//       log_lik[i] = normal_lpdf(z[i]| mu, sigma + mu*obs_err);
//     }
//   }
// }