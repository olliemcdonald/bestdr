/*
Single Type Birth-Death CLT Stan code template w/o mixed effects
Creator: Thomas O McDonald (mcdonald@jimmy.harvard.edu)

DO NOT ALTER OR WILL BREAK PACKAGE
*/
data{
  int<lower=0> N; //number of datapoints
  vector[N] count; //endpoint populations
  vector[N] count_prev; //inital population vectors
  vector[N] dt;               // time between observations
  array[N] int<lower=0> x_idx;      // index (int) of each concentration
  int<lower=0> nx;            // number of unique doses
  array[nx] real xval;              // value of each concentration

  /*******************/
  // HYPERPRIOR BLOCK
  // Includes the hyperparameters for the hyperpriors
%s
  /*******************/

}
parameters{
  /*******************/
  // PARAMETER BLOCK
  // Includes parameters and hyperparameters
%s
  /*******************/

  vector<lower=-0.5, upper=0.5>[N] count_err;  // To account for rounding to nearest cell
}
transformed parameters{
  /*******************/
  // TRANSFORMED PARAMETER BLOCK
  // Includes parameters
%s
  /*******************/
  vector<lower=0>[N] z;
  // True value is sum of rounded value (data) and roundoff error
  z = count + count_err;

  for(xi in 1:nx){
  /*******************/
  // RATE BLOCK
  // Note that the parameters that are random effects need a [g] after
%s
  /*******************/
  }
}
model{
  vector[N] sigma;
  vector[N] mu;

  /*******************/
  // PRIOR_BLOCK
  // 1.) priors on the hyperparameters
  // 2.) hyperparameter priors
%s
  /*******************/

  // analytical values for mean and variance of a single type birth death process
  for ( i in 1:N ) {
    /*******************/
    // BIRTH DEATH MODEL BLOCK
    // formula for mean and variance as a function of the rate parameters
%s
    /*******************/
  }

  /*******************/
  // BLOCK THAT CHANGES GIVEN ADDITION OF OBSERVATION ERROR
  %s
  /*******************/
}
