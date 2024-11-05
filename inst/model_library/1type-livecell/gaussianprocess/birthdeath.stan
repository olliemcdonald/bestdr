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
  array[nc] vector[nd] conc;        // value of each drug's concentration
  
  real mu_gr;
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  real mu_alpha_gr;
  real<lower=0> s_alpha_gr;
  
  real<lower=0> s_sigma_gr;
}
transformed data{
  vector[nc] m_gr = rep_vector(mu_gr, nc);
}
parameters{
  vector[nc] growth_tilde;
  real<lower=0> rho_gr;
  real<lower=0> alpha_gr;
  real<lower=0> sigma_gr;
  // To account for rounding to nearest cell
  vector<lower=-0.5, upper=0.5>[N] count_err;
}
transformed parameters{
  vector[N] z;
  vector[nc] gr;
  {
    matrix[nc, nc] K_gr = gp_exp_quad_cov(conc, alpha_gr, rho_gr) + diag_matrix(rep_vector(1e-10, nc));
    matrix[nc, nc] L_K_gr = cholesky_decompose(K_gr);
    gr = L_K_gr * growth_tilde + m_gr;
  }
  
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;
  
}
model{
  vector[N] mu;
  
  rho_gr ~ inv_gamma(alpha_rho, beta_rho);
  alpha_gr ~ normal(mu_alpha_gr, s_alpha_gr);
  sigma_gr ~ cauchy(0, s_sigma_gr);

  growth_tilde ~ std_normal();

  for ( i in 1:N ) {
    mu[i] = count_prev[i] * exp((gr[c_idx[i]]) * dt[i]);
  }
  // likelihood function
  z ~ normal( mu , sigma_gr );
}
// generated quantities{
//   vector[N] log_lik;
  
//   {
//     vector[N] mu;
    
//     for ( i in 1:N ) {
//       mu[i] = count_prev[i] * exp((gr[c_idx[i]]) * dt[i]);
//       log_lik[i] = normal_lpdf(z[i] | mu[i], sigma_gr);
//     }
//   }
// }
