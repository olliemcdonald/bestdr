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
      b_mat[i] = (r_mat[,i]'*e_mat);
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
  int<lower=0> n_dt_unique; //number of distinct timepoints.
  
  matrix<lower=0>[N,ntypes] count; //endpoint populations
  matrix<lower=0>[N,ntypes] count_prev; //inital population vectors
  matrix<lower=0>[nevents,ntypes] e_mat; //birth events matrix
  array[nevents] int<lower=0> p_vec; //parents for each birth event
  array[N] int<lower=0> times_idx; //index of the duration of each run
  array[n_dt_unique] real<lower=0> dt; //the actual inter_event times
  array[N] int<lower=0> c_idx;      // index (int) of each concentration
  int<lower=0> nc;            // number of unique doses
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
  
  real prior_mu_k;            // mean for prior distribution of k
  real<lower=0> prior_s_k;    // standard deviation for prior distribution of k
  real<lower=0> prior_mu_obs_err; // location parameter for prior distribution of error
  real<lower=0> prior_s_obs_err;  // scale parameter for prior distribution of error

}
transformed data{
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
  real<lower=0> bi;
  real<lower=0> db;
  real b50;
  real<lower=0> bh;
  real<lower=0> d0;
  real<lower=0> dd;
  real d50;
  real<lower=0> dh;
	real<lower=0> k;
	real<lower=0> obs_err;
  ////matrix<lower=-0.5, upper=0.5>[N,ntypes] count_err;  // To account for rounding to nearest cell
}
transformed parameters{
  // To account for rounding to nearest cell
  //matrix[N,ntypes] z;
  array[nc, n_dt_unique] vector[ntypes*ntypes + ntypes*ntypes*ntypes] moments; //raw single-ancestor moments vector evolving over time

  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;

  { // block to make these variables local
    matrix[nevents,ntypes] r_mat;
    for(c in 1:nc){
      r_mat = rep_matrix(0, nevents,ntypes);
  
      if(conc[c] == negative_infinity())
      {
        r_mat[1, p_vec[1]] = bi + db;
        r_mat[2, p_vec[2]] = d0;
      }
      else
      {
        r_mat[1, p_vec[1]] = bi + db / (1 + exp(bh * (conc[c] - b50)));
    		r_mat[2, p_vec[2]] = d0 + dd - dd / (1 + exp(dh * (conc[c] - d50)));
      }
      r_mat[3, p_vec[3]] = k;
  
      moments[c] = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-4, 1e-3, 1000, ntypes, nevents, e_mat, r_mat);
    }
  }
}
model{
  bi ~ normal(prior_mu_bi, prior_s_bi);
  db ~ normal(prior_mu_db, prior_s_db);
  b50 ~ normal(prior_mu_b50, prior_s_b50);
  bh ~ normal(prior_mu_bh, prior_s_bh);
  d0 ~ normal(prior_mu_d0, prior_s_d0);
  dd ~ normal(prior_mu_dd, prior_s_dd);
  d50 ~ normal(prior_mu_d50, prior_s_d50);
  dh ~ normal(prior_mu_dh, prior_s_dh);
  k ~ normal(prior_mu_k, prior_s_k);
  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);

  { // block to make these variables local
    array[nc, n_dt_unique] matrix[ntypes,ntypes] m_t; //first moment matrices
    array[nc, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff
    
    for(c in 1:nc){
      for(ti in 1:n_dt_unique){
        m_t[c,ti] = to_matrix(head(moments[c,ti],ntypes*ntypes), ntypes,ntypes);
        d_t[c,ti] = to_matrix(segment(moments[c,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
      }
    }


    for(n in 1:N){
      int ti = times_idx[n];
      int c = c_idx[n];

      //plug in the inital conditions
      mu_t = (m_t[c,ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[c,ti],ntypes,ntypes);
      for(i in 1:ntypes){
        sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[c,ti],i).*col(m_t[c,ti],i)); //subtract to get covariance from 2nd moment
        for(j in 1:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[c,ti],i).*col(m_t[c,ti],j)); //subtract to get covariance from 2nd moment
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

  { // block to make these variables local
    array[nc, n_dt_unique] matrix[ntypes,ntypes] m_t; //first moment matrices
    array[nc, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff
    
    for(c in 1:nc){
      for(ti in 1:n_dt_unique){
        m_t[c,ti] = to_matrix(head(moments[c,ti],ntypes*ntypes), ntypes,ntypes);
        d_t[c,ti] = to_matrix(segment(moments[c,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
      }
    }


    for(n in 1:N){
      int ti = times_idx[n];
      int c = c_idx[n];

      //plug in the inital conditions
      mu_t = (m_t[c,ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[c,ti],ntypes,ntypes);
      for(i in 1:ntypes){
        sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[c,ti],i).*col(m_t[c,ti],i)); //subtract to get covariance from 2nd moment
        for(j in 1:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[c,ti],i).*col(m_t[c,ti],j)); //subtract to get covariance from 2nd moment
          sigma_t[j,i] = sigma_t[i,j];
        }
      }
  
      log_lik[n] = multi_normal_lpdf(z[n] | mu_t, sigma_t + (obs_err*diag_matrix(mu_t)).^2);
      
    }
  }
}
*/