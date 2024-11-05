/*
Multitype Branching Process CLT Stan code template w/ mixed effects
Creator: Thomas O McDonald (mcdonald@jimmy.harvard.edu)

DO NOT ALTER OR WILL BREAK PACKAGE
*/
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

  matrix[N,ntypes] count; //endpoint populations
  matrix[N,ntypes] count_prev; //inital population vectors
  matrix<lower=0>[nevents,ntypes] e_mat; //birth events matrix
  array[nevents] int<lower=0> p_vec; //parents for each birth event
  array[N] int<lower=0> times_idx; //index of the duration of each run
  array[n_dt_unique] real dt; //the actual inter_event times

  array[N] int<lower=0> x_idx;      // index (int) of each concentration
  int<lower=0> nx;            // number of unique doses
  array[nx] real xval;              // value of each concentration


  array[N] int<lower=1> g_idx;          // index for group
  int<lower=1> ng;                // number of groups

  /*******************/
  // HYPERPRIOR BLOCK
  // Includes the hyperparameters for the hyperpriors
%s
  /*******************/

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
  /*******************/
  // PARAMETER BLOCK
  // Includes parameters and hyperparameters
%s
  /*******************/

  matrix<lower=-0.5, upper=0.5>[N,ntypes] count_err;  // To account for rounding to nearest cell
}
transformed parameters{
  // To account for rounding to nearest cell
  matrix[N,ntypes] z;
  array[ng, nx, n_dt_unique] vector[ntypes*ntypes + ntypes*ntypes*ntypes] moments; //raw single-ancestor moments vector evolving over time

  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;

  { // block to make these variables local
    matrix[nevents,ntypes] r_mat;
    for(g in 1:ng){
      for(xi in 1:nx){
      r_mat = rep_matrix(0, nevents,ntypes);

      /*******************/
      // RATE BLOCK
      // Note that the parameters that are random effects need a [g] after
%s
      /*******************/

        //moments[xi] = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-4, 1e-3, 1000, ntypes, nevents, e_mat, r_mat);
        moments[g,xi] = ode_bdf_tol(moment_ode, init_state, 0, dt, 1e-3, 1e-4, 1000, ntypes, nevents, e_mat, r_mat);
      }
    }
  }
}
model{
  /*******************/
  // PRIOR_BLOCK
  // 1.) priors on the hyperparameters
  // 2.) hyperparameter priors
%s
  /*******************/

  {  // block to keep these variables local
    array[ng, nx, n_dt_unique] matrix[ntypes,ntypes] m_t; //fisrt moment matrices
    array[ng, nx, n_dt_unique] matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i)
    vector[ntypes] mu_t; //mean vectors for each datapoint
    matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
    matrix[ntypes,ntypes] temp; //for copying stuff

    for(g in 1:ng){
      for(xi in 1:nx){
        for(ti in 1:n_dt_unique){
          m_t[g,xi,ti] = to_matrix(head(moments[g,xi,ti],ntypes*ntypes), ntypes,ntypes);
          d_t[g,xi,ti] = to_matrix(segment(moments[g,xi,ti],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
        }
      }
    }

    for(n in 1:N){
      int ti = times_idx[n];
      int xi = x_idx[n];
      int g = g_idx[n];

      //plug in the inital conditions
      mu_t = (m_t[g,xi,ti]')*to_vector(count_prev[n]);
      temp = to_matrix(count_prev[n]*d_t[g,xi,ti],ntypes,ntypes);
      for(i in 1:ntypes){
        sigma_t[i,i] = temp[i,i] - count_prev[n]*(col(m_t[g,xi,ti],i).*col(m_t[g,xi,ti],i)); //subtract to get covariance from 2nd moment
        for(j in 1:ntypes){
          sigma_t[i,j] = temp[i,j] - count_prev[n]*(col(m_t[g,xi,ti],i).*col(m_t[g,xi,ti],j)); //subtract to get covariance from 2nd moment
          sigma_t[j,i] = sigma_t[i,j];
        }
      }

      /*******************/
      // BLOCK THAT CHANGES GIVEN ADDITION OF OBSERVATION ERROR
      %s
      /*******************/
    }
  }
}
