functions{
  //one giant ODE system for solving all the moements at once
  //theta = unrolled birth matrix + death matrix
  //state - unrolled vector of first moment matrix + second moment matrix
  vector moment_ode(real t, vector state, int ntypes, int nevents, matrix e_mat, matrix r_mat) {// rdata param not used
    matrix[ntypes,ntypes] a_mat;             //1st moment ode coefficients
    matrix[ntypes,ntypes] b_mat;             //individual offspring mean matrix
    array[ntypes] matrix[ntypes,ntypes] c_mat;          //array of matrices of second derivates of offspring PGFs
    matrix[ntypes,ntypes*ntypes] beta_mat;        //2nd moment ode coefficients
    matrix[ntypes,ntypes] m_t;            //mean matrix at time t
    matrix[ntypes,ntypes] m_t_prime;      //first moment derivatives
    matrix[ntypes,ntypes*ntypes] d_t;          //second moments at time t. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[ntypes,ntypes*ntypes] d_t_prime;    //second moment derivatives. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[nevents,ntypes] e_mat_tmp;
    matrix[ntypes,ntypes] b2;
    vector[ntypes] lambda;

    for(i in 1:ntypes){
      lambda[i] = sum(col(r_mat,i));
      b_mat[i] = r_mat[,i]'*e_mat;
      a_mat[i] = b_mat[i];
      a_mat[i,i] = b_mat[i,i] - lambda[i];
    }

    for(a in 1:ntypes){//the ancestor type we are multiplying by
      for(i in 1:nevents){
        e_mat_tmp[i] = e_mat[i]*e_mat[i,a];
      }

      for(i in 1:ntypes){
        b2[i] = (r_mat[,i]'*e_mat_tmp); //computing the TRANSPOSE of previous convention
        c_mat[i][a] = b2[i];
        c_mat[i][a,a] = b2[i,a] - b_mat[i,a];
      }
    }


    //unpack the moment matrices from the state vector
    m_t = to_matrix(head(state, ntypes*ntypes), ntypes, ntypes);
    d_t = to_matrix(segment(state,ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);//read this in row-major order

    for(i in 1:ntypes){
      beta_mat[i] = to_row_vector((m_t')*c_mat[i]*m_t);
    }

    //first moment ODE
    m_t_prime = a_mat*m_t;

    //second moment ODE
    d_t_prime = a_mat*d_t + beta_mat;

    //collapse derivatives into 1d array and return
    return append_row(to_vector(m_t_prime),to_vector(d_t_prime));
  }
}
data{
  int<lower=0> ntypes; //number of types
  int<lower=0> nevents; //number of events
  int<lower=0> N; //number of datapoints
  
  matrix<lower=0>[N,ntypes] count; //endpoint populations
  matrix<lower=0>[N,ntypes] count_prev; //inital population vectors
  matrix<lower=0>[nevents,ntypes] e_mat; //birth events matrix
  array[nevents] int<lower=0> p_vec; //parents for each birth event
  
  int<lower=0> n_dt_unique; //number of distinct timepoints.
  array[N] int<lower=0> times_idx; //index of the duration of each run
  array[n_dt_unique] real<lower=0> dt; //the actual inter_event times
  
 array[N] int<lower=0> c_idx;      // index (int) of each concentration
  int<lower=0> nc;            // number of unique doses
  array[nc] real conc;              // value of each concentration
  
  array[N] int<lower=1> g_idx;          // index for group
  int<lower=1> ng;                // number of groups
  
  real mu_b;
  real mu_d;

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

  real theta_1_mu_k;              // theta_1 for prior distribution of mean of k
  real<lower=0> theta_2_mu_k;     // theta_2 for prior distribution of mean of k
  real theta_1_s_k;               // theta_1 for prior distribution of stdev of k
  real<lower=0> theta_2_s_k;      // theta_2 for prior distribution of stdev of k
  
  // Parameter for error - assume this is the same for all cell lines
  real<lower=0> s_obs_err;  // scale parameter for prior distribution of error
}
transformed data{
  array[ng] vector[nc] m_b;
  array[ng] vector[nc] m_d;

  vector[ntypes*ntypes + ntypes*ntypes*ntypes] init_state; //state is unrolled vector of first + second moments
  vector[ntypes * ntypes] m_init = rep_vector(0.0,ntypes*ntypes);
  vector[ntypes*ntypes*ntypes] d_init = rep_vector(0.0,ntypes*ntypes*ntypes);

  //evaluations of the moments at t=0. Mostly zeros, since type a -> b_mat is impossible in 0 time
  for(i in 1:ntypes){
    m_init[i+(i-1)*ntypes] = 1;
    d_init[i+(i-1)*ntypes*(ntypes+1)] = 1;
  }
  init_state = append_row(m_init,d_init);
}
parameters{
  // parameters
  array[ng] vector[nc] birth_tilde;
  array[ng] vector[nc] death_tilde;
  vector[ng] k_raw;
  array[ng] real<lower=0> rho_b;
  array[ng] real<lower=0> alpha_b;
  array[ng] real<lower=0> rho_d;
  array[ng] real<lower=0> alpha_d;
  real mu_k;
  real<lower=0> s_k;
  real<lower=0> obs_err;

  // hyperparameters
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  real mu_alpha_b;
  real<lower=0> s_alpha_b;
  real mu_alpha_d;
  real<lower=0> s_alpha_d;
  real mu_alpha_k;
  real<lower=0> s_alpha_k;
  //matrix<lower=-0.5, upper=0.5>[N,ntypes] count_err;  // To account for rounding to nearest cell

}
transformed parameters{
  array[ng]vector<lower=0>[nc] b;
  array[ng] vector<lower=0>[nc] d;
  vector<lower=0>[ng] k;
  //matrix[N,ntypes] z;

  array[ng, nc, n_dt_unique] vector[ntypes*ntypes + ntypes*ntypes*ntypes] moments; //raw single-ancestor moments vector evolving over time

  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;
  
  k = mu_k + s_k * k_raw;

  {  // block to keep variables local
    matrix[nevents,ntypes] r_mat;
    array[nevents*2*ntypes] real theta;
    matrix[nc, nc] K_b;
    matrix[nc, nc] L_K_b;
    matrix[nc, nc] K_d;
    matrix[nc, nc] L_K_d;
    
    for(g in 1:ng){
      K_b = gp_exp_quad_cov(conc, alpha_b[g], rho_b[g]) + diag_matrix(rep_vector(1e-10, nc));
      L_K_b = cholesky_decompose(K_b);
      K_d = gp_exp_quad_cov(conc, alpha_d[g], rho_d[g]) + diag_matrix(rep_vector(1e-10, nc));
      L_K_d = cholesky_decompose(K_d);
      
      b[g] = L_K_b * birth_tilde[g] + m_b[g];
      d[g] = L_K_d * death_tilde[g] + m_d[g];
      
      for(c in 1:nc){
        r_mat = rep_matrix(0, nevents,ntypes);
        r_mat[1, p_vec[1]] = b[g,c];
        r_mat[2, p_vec[2]] = d[g,c];
        r_mat[3, p_vec[3]] = k[g];
        moments[g,c] = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-4, 1e-3, 1000, ntypes, nevents, e_mat, r_mat);
      }
    }
  }
}
model{
  // specified prior distributions
  alpha_rho ~ normal(theta_1_alpha_rho, theta_2_alpha_rho);
  beta_rho ~ normal(theta_1_beta_rho, theta_2_beta_rho);
  mu_alpha_b ~ normal(theta_1_mu_alpha_b, theta_2_mu_alpha_b);
  s_alpha_b ~ normal(theta_1_s_alpha_b, theta_2_s_alpha_b);
  mu_alpha_d ~ normal(theta_1_mu_alpha_d, theta_2_mu_alpha_d);
  s_alpha_d ~ normal(theta_1_s_alpha_d, theta_2_s_alpha_d);
  mu_k ~ normal(theta_1_mu_k, theta_2_mu_k);
  s_k ~ normal(theta_1_s_k, theta_2_s_k);
  
  rho_b ~ inv_gamma(alpha_rho, beta_rho);
  alpha_b ~ normal(mu_alpha_b, s_alpha_b);
  rho_d ~ inv_gamma(alpha_rho, beta_rho);
  alpha_d ~ normal(mu_alpha_d, s_alpha_d);
  k_raw ~ std_normal();
  
  obs_err ~ normal(0, s_obs_err);

  {  // block to keep these variables local
    array[ng, nc, n_dt_unique] matrix[ntypes,ntypes] m_t; //first moment matrices
    array[ng, nc, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff

    for(g in 1:ng){
      for(c in 1:nc){
        for(ti in 1:n_dt_unique){
          m_t[g,c,ti] = to_matrix(head(moments[g,c,ti],ntypes*ntypes), ntypes,ntypes);
          d_t[g,c,ti] = to_matrix(segment(moments[g,c,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
        }
      }
    }

    for(n in 1:N){
      int ti = times_idx[n];
      int c = c_idx[n];
      int g = g_idx[n];

      //plug in the inital conditions
      mu_t = (m_t[g,c,ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[g,c,ti],ntypes,ntypes);
      for(i in 1:ntypes){
        sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],i)); //subtract to get covariance from 2nd moment
        for(j in 1:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],j)); //subtract to get covariance from 2nd moment
          sigma_t[j,i] = sigma_t[i,j];
        }
      }
  
      count[n] ~ multi_normal(mu_t, sigma_t + (obs_err*diag_matrix(mu_t)).^2);
    }
  }
}
/*
generated quantities{
  array[N] real log_lik;

  {  // block to keep these variables local
    array[ng, nc, n_dt_unique] matrix[ntypes,ntypes] m_t; //first moment matrices
    array[ng, nc, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff

    for(g in 1:ng){
      for(c in 1:nc){
        for(ti in 1:n_dt_unique){
          m_t[g,c,ti] = to_matrix(head(moments[g,c,ti],ntypes*ntypes), ntypes,ntypes);
          d_t[g,c,ti] = to_matrix(segment(moments[g,c,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
        }
      }
    }

    for(n in 1:N){
      int ti = times_idx[n];
      int c = c_idx[n];
      int g = g_idx[n];

      //plug in the inital conditions
      mu_t = (m_t[g,c,ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[g,c,ti],ntypes,ntypes);
      for(i in 1:ntypes){
        sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],i)); //subtract to get covariance from 2nd moment
        for(j in 1:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],j)); //subtract to get covariance from 2nd moment
          sigma_t[j,i] = sigma_t[i,j];
        }
      }
  
      log_lik[n] = multi_normal_lpdf(z[n] | mu_t, sigma_t + (obs_err*diag_matrix(mu_t)).^2);
    }
  }
}
*/