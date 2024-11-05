// Logistic Birth-Death Model with birth and death rates as parameters
// No Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector<lower=0>[N] count;            // total count at current time
  vector<lower=0>[N] count_prev;       // total count at previous time
  vector<lower=0>[N] dt;               // time between observations
  array[N] int<lower=1> c_idx;      // index (int) of each concentration
  int<lower=1> nc;            // number of unique doses
  array[nc] real conc;              // value of each concentration
  
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
}
parameters{
  //  using log birth rate and log death rate as priors since has better convergence issues
  real<lower=0> bi;
  real<lower=0> db;
  real b50;
  real<lower=0> bh;
  real<lower=0> d0;
  real<lower=0> dd;
  real d50;
  real<lower=0> dh;
  // To account for rounding to nearest cell
  vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  vector<lower=0>[nc] b;
  vector<lower=0>[nc] d;
  vector<lower=0>[N] z;
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;
  

  for (i in 1:nc)
  {
    if (conc[i] == negative_infinity())
    {
      b[i] = bi + db;
      d[i] = d0;
    }
    else
    {
      b[i] = bi + db / (1 + exp(bh * (conc[i] - b50)));
      d[i] = d0 + dd - dd / (1 + exp(dh * (conc[i] - d50)));
    }
  }
}
model{
  vector[N] sigma;
  vector[N] mu;

  // specified prior distributions
  bi ~ normal(prior_mu_bi, prior_s_bi) T[0, ];
  db ~ normal(prior_mu_db, prior_s_db) T[0, ];
  b50 ~ normal(prior_mu_b50, prior_s_b50);
  bh ~ normal(prior_mu_bh, prior_s_bh) T[0, ];
  d0 ~ normal(prior_mu_d0, prior_s_d0) T[0, ];
  dd ~ normal(prior_mu_dd, prior_s_dd) T[0, ];
  d50 ~ normal(prior_mu_d50, prior_s_d50);
  dh ~ normal(prior_mu_dh, prior_s_dh) T[0, ];

  // analytical values for mean and variance of a single type birth death process
  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b[c_idx[i]]+d[c_idx[i]]) / (b[c_idx[i]]-d[c_idx[i]]) * (exp(2*(b[c_idx[i]]-d[c_idx[i]])*dt[i]) - exp((b[c_idx[i]]-d[c_idx[i]])*dt[i])));
    mu[i] = count_prev[i] * exp((b[c_idx[i]] - d[c_idx[i]]) * dt[i]);
  }
  // likelihood function
  z ~ normal( mu , sigma );
}
// generated quantities{
//   vector[N] log_lik;
  
//   {
//     real sigma;
//     real mu;
//     for (i in 1:N)
//     {
//       sigma = sqrt(count_prev[i] * (b[c_idx[i]]+d[c_idx[i]]) / (b[c_idx[i]]-d[c_idx[i]]) * (exp(2*(b[c_idx[i]]-d[c_idx[i]])*dt[i]) - exp((b[c_idx[i]]-d[c_idx[i]])*dt[i])));
//       mu = count_prev[i] * exp((b[c_idx[i]] - d[c_idx[i]]) * dt[i]);
//       log_lik[i] = normal_lpdf(z[i]| mu, sigma);
//     }
//   }
// }