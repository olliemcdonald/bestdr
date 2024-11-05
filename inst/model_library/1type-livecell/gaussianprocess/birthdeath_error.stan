// Gaussian Process Birth-Death Model for N unique drugs in combination
// with birth and death rates as parameters
// No Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;               // time between observations
  array[N] int<lower=1> c_idx;      // index (int) of each concentration
  int<lower=1> nc;            // number (of unique combos) of concentrations
  int<lower=1> nd;            // number of unique drugs / treatment groups
  array[nc] vector[nd] conc;  // value of each drug's concentration - each conc
                              // is a vector containing concentration of each drug
  // Priors
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  real mu_alpha_b;
  real<lower=0> s_alpha_b;
  real mu_alpha_d;
  real<lower=0> s_alpha_d;

  real<lower=0> prior_mu_obs_err;  // mean parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error

}
parameters{
  vector[nc] birth_tilde;
  vector[nc] death_tilde;
  real<lower=0> rho_b;
  real<lower=0> alpha_b;
  real<lower=0> rho_d;
  real<lower=0> alpha_d;
  real<lower=0> obs_err;
  // To account for rounding to nearest cell
  //vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  vector<lower=0>[nc] b;
  vector<lower=0>[nc] d;
  //vector[N] z;
  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;
  
  { // Keep matrices local so they don't get stored
    matrix[nc, nc] K_b = gp_exp_quad_cov(conc, alpha_b, rho_b) + diag_matrix(rep_vector(1e-8, nc));
    matrix[nc, nc] L_K_b = cholesky_decompose(K_b);
    matrix[nc, nc] K_d = gp_exp_quad_cov(conc, alpha_d, rho_d) + diag_matrix(rep_vector(1e-8, nc));
    matrix[nc, nc] L_K_d = cholesky_decompose(K_d);
    
    b = exp(L_K_b * birth_tilde);
    d = exp(L_K_d * death_tilde);
  }
  
}
model{
  vector[N] sigma;
  vector[N] mu;
  
  rho_b ~ inv_gamma(alpha_rho, beta_rho);
  alpha_b ~ normal(mu_alpha_b, s_alpha_b);
  rho_d ~ inv_gamma(alpha_rho, beta_rho);
  alpha_d ~ normal(mu_alpha_d, s_alpha_d);

  birth_tilde ~ std_normal();
  death_tilde ~ std_normal();

  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);


  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b[c_idx[i]]+d[c_idx[i]]) / (b[c_idx[i]]-d[c_idx[i]]) * (exp(2*(b[c_idx[i]]-d[c_idx[i]])*dt[i]) - exp((b[c_idx[i]]-d[c_idx[i]])*dt[i])));
    mu[i] = count_prev[i] * exp((b[c_idx[i]] - d[c_idx[i]]) * dt[i]);
  }
  // likelihood function
  count ~ normal( mu , sqrt(sigma.^2.0 + (obs_err * mu).^2.0) );
}
/* generated quantities{
  vector[N] log_lik;
  vector<lower=0>[nc_predict] b_predict;
  vector<lower=0>[nc_predict] d_predict;
  
  {
    real sigma;
    real mu;
    
    for ( i in 1:N ) {
      sigma = sqrt(count_prev[i] * (b[c_idx[i]]+d[c_idx[i]]) / (b[c_idx[i]]-d[c_idx[i]]) * (exp(2*(b[c_idx[i]]-d[c_idx[i]])*dt[i]) - exp((b[c_idx[i]]-d[c_idx[i]])*dt[i])));
      mu = count_prev[i] * exp((b[c_idx[i]] - d[c_idx[i]]) * dt[i]);
      log_lik[i] = normal_lpdf(z[i] | mu, sigma);
    }
  }

}*/
