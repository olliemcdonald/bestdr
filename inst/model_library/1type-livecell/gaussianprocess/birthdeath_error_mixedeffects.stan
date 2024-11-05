// Gaussian Process Birth-Death Model for N unique drugs in combination
// with birth and death rates as parameters
// Mixed Effects
data{
  int<lower=1> N;             // number of samples
  vector[N] count;            // total count at current time
  vector[N] count_prev;       // total count at previous time
  vector[N] dt;               // time between observations
  array[N] int<lower=1> c_idx;      // index (int) of each concentration
  int<lower=1> nc;            // number (of unique combos) of concentrations
  int<lower=1> nd;            // number of unique drugs / treatment groups
  array[nc] vector[nd] conc;        // value of each drug's concentration - each conc
                              // is a vector containing concentration of each drug
  array[N] int<lower=1> g_idx;      // group id (random effect)
  int<lower=1> ng;            // number of groups
  
  real theta_1_alpha_rho;
  real theta_2_alpha_rho;
  real theta_1_beta_rho;
  real theta_2_beta_rho;
  
  real theta_1_mu_alpha_b;
  real theta_2_mu_alpha_b;
  real theta_1_s_alpha_b;
  real theta_2_s_alpha_b;
  real theta_1_mu_alpha_d;
  real theta_2_mu_alpha_d;
  real theta_1_s_alpha_d;
  real theta_2_s_alpha_d;

  real<lower=0> prior_mu_obs_err;  // mean parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error
}
parameters{
  array[ng] vector[nc] birth_tilde;
  array[ng] vector[nc] death_tilde;
  array[ng] real<lower=0> rho_b;
  array[ng] real<lower=0> alpha_b;
  array[ng] real<lower=0> rho_d;
  array[ng] real<lower=0> alpha_d;
  real<lower=0> obs_err;
  
  // hyperparameters
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  real mu_alpha_b;
  real<lower=0> s_alpha_b;
  real mu_alpha_d;
  real<lower=0> s_alpha_d;

}
transformed parameters{
  array[ng] vector<lower=0>[nc] b;
  array[ng] vector<lower=0>[nc] d;

  
  {  // Keep matrices local so they don't get stored
    matrix[nc, nc] K_b;
    matrix[nc, nc] L_K_b;
    matrix[nc, nc] K_d;
    matrix[nc, nc] L_K_d;
    
    for(g in 1:ng){
      K_b = gp_exp_quad_cov(conc, alpha_b[g], rho_b[g]) + diag_matrix(rep_vector(1e-8, nc));
      K_d = gp_exp_quad_cov(conc, alpha_d[g], rho_d[g]) + diag_matrix(rep_vector(1e-8, nc));
      L_K_b = cholesky_decompose(K_b);
      L_K_d = cholesky_decompose(K_d);
      
      b[g] = exp(L_K_b * birth_tilde[g]);
      d[g] = exp(L_K_d * death_tilde[g]);
    }
  }
  
}
model{
  vector[N] sigma;
  vector[N] mu;
  
  alpha_rho ~ normal(theta_1_alpha_rho, theta_2_alpha_rho);
  beta_rho ~ normal(theta_1_beta_rho, theta_2_beta_rho);
  mu_alpha_b ~ normal(theta_1_mu_alpha_b, theta_2_mu_alpha_b);
  s_alpha_b ~ normal(theta_1_s_alpha_b, theta_2_s_alpha_b);
  mu_alpha_d ~ normal(theta_1_mu_alpha_d, theta_2_mu_alpha_d);
  s_alpha_d ~ normal(theta_1_s_alpha_d, theta_2_s_alpha_d);
  
  rho_b ~ inv_gamma(alpha_rho, beta_rho);
  alpha_b ~ normal(mu_alpha_b, s_alpha_b);
  rho_d ~ inv_gamma(alpha_rho, beta_rho);
  alpha_d ~ normal(mu_alpha_d, s_alpha_d);

  for(g in 1:ng){
    birth_tilde[g] ~ std_normal();
    death_tilde[g] ~ std_normal();
  }

  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);
  
  for ( i in 1:N ) {
    sigma[i] = sqrt(count_prev[i] * (b[g_idx[i], c_idx[i]]+d[g_idx[i], c_idx[i]]) / (b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]]) * (exp(2*(b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]])*dt[i]) - exp((b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]])*dt[i])));
    mu[i] = count_prev[i] * exp((b[g_idx[i], c_idx[i]] - d[g_idx[i], c_idx[i]]) * dt[i]);
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
    for ( i in 1:N ) {
      sigma = sqrt(count_prev[i] * (b[g_idx[i], c_idx[i]]+d[g_idx[i], c_idx[i]]) / (b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]]) * (exp(2*(b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]])*dt[i]) - exp((b[g_idx[i], c_idx[i]]-d[g_idx[i], c_idx[i]])*dt[i])));
      mu = count_prev[i] * exp((b[g_idx[i], c_idx[i]] - d[g_idx[i], c_idx[i]]) * dt[i]);
      log_lik[i] = normal_lpdf(z[i] | mu, sigma);
    }
  }
}
*/
