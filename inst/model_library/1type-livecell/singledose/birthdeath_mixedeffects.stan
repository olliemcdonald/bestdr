//  Single Dose Birth-Death Model with birth and death rates as parameters
//  Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;               // time between observations
  array[N] int<lower=1> g;          // group id (random effect)
  int<lower=1> ng;            // number of groups
  
  // Parameters for priors of hyperparameters
  real theta_1_mu_b;   // theta_1 for prior distribution of mean of b
  real theta_2_mu_b;   // theta_2 for prior distribution of mean of b
  real theta_1_mu_d;   // theta_1 for prior distribution of mean of d
  real theta_2_mu_d;   // theta_2 for prior distribution of mean of d
  real theta_1_s_b;    // theta_1 for prior distribution of stdev of b
  real theta_2_s_b;    // theta_2 for prior distribution of stdev of b
  real theta_1_s_d;    // theta_1 for prior distribution of stdev of d
  real theta_2_s_d;    // theta_2 for prior distribution of stdev of d
}
parameters{
  vector<lower=0>[ng] b;
  vector<lower=0>[ng] d;
  // To account for rounding to nearest cell
  vector<lower=-0.5, upper=0.5>[N] count_err;
  // hyperparameters
  real<lower=0> mu_b;
  real<lower=0> mu_d;
  real<lower=0> s_b;
  real<lower=0> s_d;
}
transformed parameters{
  vector[N] z;
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;
}
model{
  vector[N] sigma;
  vector[N] mu;

  // specified hyperprior distributions
  mu_b ~ normal(theta_1_mu_b, theta_2_mu_b);
  s_b ~ normal(theta_1_s_b, theta_2_s_b);
  mu_d ~ normal(theta_1_mu_d, theta_2_mu_d);
  s_d ~ normal(theta_1_s_d, theta_2_s_d);
  
  // b and d priors
  b ~ normal(mu_b, s_b);
  d ~ normal(mu_d, s_d);

  // analytical values for mean and variance of a single type birth death process
  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b[g[i]]+d[g[i]]) / (b[g[i]]-d[g[i]]) *
      (exp(2*(b[g[i]]-d[g[i]])*dt[i]) - exp((b[g[i]]-d[g[i]])*dt[i])));
    mu[i] = count_prev[i] * exp((b[g[i]] - d[g[i]]) * dt[i]);
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
//       sigma = sqrt(count_prev[i] * (b[g[i]]+d[g[i]]) / (b[g[i]]-d[g[i]]) *
//         (exp(2*(b[g[i]]-d[g[i]])*dt[i]) - exp((b[g[i]]-d[g[i]])*dt[i])));
//       mu = count_prev[i] * exp((b[g[i]] - d[g[i]]) * dt[i]);
//       log_lik[i] = normal_lpdf(z[i]| mu, sigma);
//     }
//   }
// }