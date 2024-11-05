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

    for(anc in 1:ntypes){//the ancestor type we are multiplying by
      for(i in 1:nevents){
        e_mat_tmp[i] = e_mat[i]*e_mat[i,anc];
      }

      for(i in 1:ntypes){
        b2[i] = (r_mat[,i]'*e_mat_tmp); //computing the TRANSPOSE of previous convention
        c_mat[i][anc] = b2[i];
        c_mat[i][anc,anc] = b2[i,anc] - b_mat[i,anc];
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
  
  matrix[N,ntypes] count; //endpoint populations
  matrix[N,ntypes] count_prev; //inital population vectors
  matrix<lower=0>[nevents,ntypes] e_mat; //birth events matrix
  array[nevents] int<lower=0> p_vec; //parents for each birth event
  array[N] int<lower=0> times_idx; //index of the duration of each run
  array[n_dt_unique] real dt; //the actual inter_event times

  vector[nevents] prior_mu;
  vector<lower=0>[nevents] prior_s;

  real<lower=0> prior_mu_obs_err;  // location parameter for prior distribution of error
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
  vector<lower=0>[nevents] theta;
	real<lower=0> obs_err;
  //matrix<lower=-0.5, upper=0.5>[N,ntypes] count_err;  // To account for rounding to nearest cell
}
transformed parameters{
  //matrix[N,ntypes] z;
  array[n_dt_unique] vector[ntypes*ntypes + ntypes*ntypes*ntypes] moments; //raw single-ancestor moments vector evolving over time

  // True value is sum of rounded value (data) and roundoff error
  //z = count + count_err;

  { // make variables local so don't save them
    matrix[nevents,ntypes] r_mat = rep_matrix(0, nevents,ntypes);
    for(i in 1:nevents){
      r_mat[i, p_vec[i]] = theta[i];
    }
    moments = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-6, 1e-6, 1000, ntypes, nevents, e_mat, r_mat);
  }
}
model{
  for(i in 1:nevents){
    theta[i] ~ normal(prior_mu[i], prior_s[i]);
  }
  obs_err ~ normal(prior_mu_obs_err, prior_s_obs_err);

  { // local to this block
    array[n_dt_unique] matrix[ntypes,ntypes] m_t; //fisrt moment matrices
    array[n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff

    for(ti in 1:n_dt_unique){
      m_t[ti] = to_matrix(head(moments[ti],ntypes*ntypes), ntypes,ntypes);
      d_t[ti] = to_matrix(segment(moments[ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
    }

    for(n in 1:N){
      int ti = times_idx[n];
  
      //plug in the inital conditions
      mu_t = (m_t[ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[ti],ntypes,ntypes);
      for(i in 1:ntypes){
          sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[ti],i).*col(m_t[ti],i)); //subtract to get covariance from 2nd moment
        for(j in i:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[ti],i).*col(m_t[ti],j)); //subtract to get covariance from 2nd moment
          sigma_t[j,i] = sigma_t[i,j];
        }
      }
  
      count[n] ~ multi_normal(mu_t, (sigma_t) + (obs_err*diag_matrix(mu_t)).^2);
    }
  }
}