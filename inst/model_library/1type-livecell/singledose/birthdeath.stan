// Single Dose Birth-Death Model with birth and death rates as parameters
// No Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;                // time between observations
  real prior_mu_b;            // mean for prior distribution of lb
  real prior_mu_d;            // mean for prior distribution of ld
  real<lower=0> prior_s_b;    // standard deviation for prior distribution of lb
  real<lower=0> prior_s_d;    // standard deviation for prior distribution of ld
}
parameters{
  //  using log birth rate and log death rate as priors since has better convergence issues
  real<lower=0> b;
  real<lower=0> d;
  // To account for rounding to nearest cell
  vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  vector[N] z;
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;
}
model{
  vector[N] sigma;
  vector[N] mu;

  // specified prior distributions
  b ~ normal( prior_mu_b, prior_s_b );
  d ~ normal( prior_mu_d, prior_s_d );

  // analytical values for mean and variance of a single type birth death process
  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b+d) / (b-d) * (exp(2*(b-d)*dt[i]) - exp((b-d)*dt[i])));
    mu[i] = count_prev[i] * exp((b - d) * dt[i]);
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
//       sigma = sqrt(count_prev[i] * (b+d) / (b-d) * (exp(2*(b-d)*dt[i]) - exp((b-d)*dt[i])));
//       mu = count_prev[i] * exp((b - d) * dt[i]);
//       log_lik[i] = normal_lpdf(z[i]| mu, sigma);
//     }
//   }
  
// }