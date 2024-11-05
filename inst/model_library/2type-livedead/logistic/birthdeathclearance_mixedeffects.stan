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
  
  // Parameters for priors of hyperparameters
  real theta_1_mu_bi;              // theta_1 for prior distribution of mean of bi
  real<lower=0> theta_2_mu_bi;     // theta_2 for prior distribution of mean of bi
  real theta_1_s_bi;               // theta_1 for prior distribution of stdev of bi
  real<lower=0> theta_2_s_bi;      // theta_2 for prior distribution of stdev of bi
  
  real theta_1_mu_db;              // theta_1 for prior distribution of mean of db
  real<lower=0> theta_2_mu_db;     // theta_2 for prior distribution of mean of db
  real theta_1_s_db;               // theta_1 for prior distribution of stdev of db
  real<lower=0> theta_2_s_db;      // theta_2 for prior distribution of stdev of db
  
  real theta_1_mu_b50;              // theta_1 for prior distribution of mean of b50
  real<lower=0> theta_2_mu_b50;     // theta_2 for prior distribution of mean of b50
  real theta_1_s_b50;               // theta_1 for prior distribution of stdev of b50
  real<lower=0> theta_2_s_b50;      // theta_2 for prior distribution of stdev of b50
  
  real theta_1_mu_bh;              // theta_1 for prior distribution of mean of bh
  real<lower=0> theta_2_mu_bh;     // theta_2 for prior distribution of mean of bh
  real theta_1_s_bh;               // theta_1 for prior distribution of stdev of bh
  real<lower=0> theta_2_s_bh;      // theta_2 for prior distribution of stdev of bh
  
  real theta_1_mu_d0;              // theta_1 for prior distribution of mean of d0
  real<lower=0> theta_2_mu_d0;     // theta_2 for prior distribution of mean of d0
  real theta_1_s_d0;               // theta_1 for prior distribution of stdev of d0
  real<lower=0> theta_2_s_d0;      // theta_2 for prior distribution of stdev of d0
  
  real theta_1_mu_dd;              // theta_1 for prior distribution of mean of dd
  real<lower=0> theta_2_mu_dd;     // theta_2 for prior distribution of mean of dd
  real theta_1_s_dd;               // theta_1 for prior distribution of stdev of dd
  real<lower=0> theta_2_s_dd;      // theta_2 for prior distribution of stdev of dd
  
  real theta_1_mu_d50;              // theta_1 for prior distribution of mean of d50
  real<lower=0> theta_2_mu_d50;     // theta_2 for prior distribution of mean of d50
  real theta_1_s_d50;               // theta_1 for prior distribution of stdev of d50
  real<lower=0> theta_2_s_d50;      // theta_2 for prior distribution of stdev of d50
  
  real theta_1_mu_dh;              // theta_1 for prior distribution of mean of dh
  real<lower=0> theta_2_mu_dh;     // theta_2 for prior distribution of mean of dh
  real theta_1_s_dh;               // theta_1 for prior distribution of stdev of dh
  real<lower=0> theta_2_s_dh;      // theta_2 for prior distribution of stdev of dh

  real theta_1_mu_k;              // theta_1 for prior distribution of mean of k
  real<lower=0> theta_2_mu_k;     // theta_2 for prior distribution of mean of k
  real theta_1_s_k;               // theta_1 for prior distribution of stdev of k
  real<lower=0> theta_2_s_k;      // theta_2 for prior distribution of stdev of k
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
    // parameters
  vector<lower=0>[ng] bi;
  vector<lower=0>[ng] db;
  vector[ng] b50;
  vector<lower=0>[ng] bh;
  vector<lower=0>[ng] d0;
  vector<lower=0>[ng] dd;
  vector[ng] d50;
  vector<lower=0>[ng] dh;
  vector<lower=0>[ng] k;
  
  // hyperparameters
  real mu_bi;
  real<lower=0> s_bi;
  real mu_db;
  real<lower=0> s_db;
  real mu_b50;
  real<lower=0> s_b50;
  real mu_bh;
  real<lower=0> s_bh;
  real mu_d0;
  real<lower=0> s_d0;
  real mu_dd;
  real<lower=0> s_dd;
  real mu_d50;
  real<lower=0> s_d50;
  real mu_dh;
  real<lower=0> s_dh;
  real mu_k;
  real<lower=0> s_k;
  
  matrix<lower=-0.5, upper=0.5>[N,ntypes] count_err;  // To account for rounding to nearest cell
}
transformed parameters{  
  // To account for rounding to nearest cell
  matrix[N,ntypes] z;

  array[ng, nc, n_dt_unique] vector[ntypes*ntypes + ntypes*ntypes*ntypes] moments; //raw single-ancestor moments vector evolving over time
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;

  {  // block to keep variables local
    matrix[nevents,ntypes] r_mat;
    for(g in 1:ng){
      for(c in 1:nc){
        r_mat = rep_matrix(0, nevents,ntypes);
  
        if(conc[c] == negative_infinity())
        {
          r_mat[1, p_vec[1]] = bi[g] + db[g];
          r_mat[2, p_vec[2]] = d0[g];
        }
        else
        {
          r_mat[1, p_vec[1]] = bi[g] + db[g] / (1 + exp(bh[g] * (conc[c]- b50[g])));
      		r_mat[2, p_vec[2]] = d0[g] + dd[g] - dd[g] / (1 + exp(dh[g] * (conc[c] - d50[g])));
        }
        r_mat[3, p_vec[3]] = k[g];
    
        moments[g,c] = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-3, 1e-4, 1000, ntypes, nevents, e_mat, r_mat);
      }
    }
  }
}
model{
  // specified prior distributions
  mu_bi ~ normal(theta_1_mu_bi, theta_2_mu_bi);
  s_bi ~ normal(theta_1_s_bi, theta_2_s_bi);
  mu_db ~ normal(theta_1_mu_db, theta_2_mu_db);
  s_db ~ normal(theta_1_s_db, theta_2_s_db);
  mu_b50 ~ normal(theta_1_mu_b50, theta_2_mu_b50);
  s_b50 ~ normal(theta_1_s_b50, theta_2_s_b50);
  mu_bh ~ normal(theta_1_mu_bh, theta_2_mu_bh);
  s_bh ~ normal(theta_1_s_bh, theta_2_s_bh);
  mu_d0 ~ normal(theta_1_mu_d0, theta_2_mu_d0);
  s_d0 ~ normal(theta_1_s_d0, theta_2_s_d0);
  mu_dd ~ normal(theta_1_mu_dd, theta_2_mu_dd);
  s_dd ~ normal(theta_1_s_dd, theta_2_s_dd);
  mu_d50 ~ normal(theta_1_mu_d50, theta_2_mu_d50);
  s_d50 ~ normal(theta_1_s_d50, theta_2_s_d50);
  mu_dh ~ normal(theta_1_mu_dh, theta_2_mu_dh);
  s_dh ~ normal(theta_1_s_dh, theta_2_s_dh);
  mu_k ~ normal(theta_1_mu_k, theta_2_mu_k);
  s_k ~ normal(theta_1_s_k, theta_2_s_k);
  
  bi ~ normal(mu_bi, s_bi);
  db ~ normal(mu_db, s_db);
  b50 ~ normal(mu_b50, s_b50);
  bh ~ normal(mu_bh, s_bh);
  d0 ~ normal(mu_d0, s_d0);
  dd ~ normal(mu_dd, s_dd);
  d50 ~ normal(mu_d50, s_d50);
  dh ~ normal(mu_dh, s_dh);
  k ~ normal(mu_k, s_k);


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
  
      z[n] ~ multi_normal(mu_t, sigma_t);
    }
  }
}
// generated quantities{
//   array[N] real log_lik;

//   {  // block to keep these variables local
//     array[ng, nc, n_dt_unique] matrix[ntypes,ntypes] m_t; //first moment matrices
//     array[ng, nc, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
//     vector[ntypes] mu_t; //mean vectors for each datapoint
//     matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
//     matrix[ntypes,ntypes] temp; //for copying stuff

//     for(g in 1:ng){
//       for(c in 1:nc){
//         for(ti in 1:n_dt_unique){
//           m_t[g,c,ti] = to_matrix(head(moments[g,c,ti],ntypes*ntypes), ntypes,ntypes);
//           d_t[g,c,ti] = to_matrix(segment(moments[g,c,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
//         }
//       }
//     }

//     for(n in 1:N){
//       int ti = times_idx[n];
//       int c = c_idx[n];
//       int g = g_idx[n];

//       //plug in the inital conditions
//       mu_t = (m_t[g,c,ti]')*to_vector(count_prev[n]);
//       temp = to_matrix(count_prev[n]*d_t[g,c,ti],ntypes,ntypes);
//       for(i in 1:ntypes){
//         sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],i)); //subtract to get covariance from 2nd moment
//         for(j in 1:ntypes){
//           sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[g,c,ti],i).*col(m_t[g,c,ti],j)); //subtract to get covariance from 2nd moment
//           sigma_t[j,i] = sigma_t[i,j];
//         }
//       }
  
//       log_lik[n] = multi_normal_lpdf(z[n] | mu_t, sigma_t);
//     }
//   }
// }