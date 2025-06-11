#include /include/common_functions.stan

data {
  // data
  // are individual observation or summary statistics used
  int is_ss;

  // sample sizes
  int<lower=0> N;
  // individual observations
  vector[is_ss == 0 ? N : 0] x;
  // summary statistics
  vector[is_ss == 1 ? 1 : 0] ss_mean;
  vector[is_ss == 1 ? 1 : 0] ss_sd;

  // model type (1 if the parameter is estimated, 0 if parameter is fixed)
  int estimate_mu;
  int estimate_sigma;

  // range of the parameters
  vector[estimate_mu    == 1 ? 2 : 0] bounds_mu;
  vector[estimate_sigma == 1 ? 2 : 0] bounds_sigma;
  array[estimate_mu    == 1 ? 2 : 0] int bounds_type_mu;
  array[estimate_sigma == 1 ? 2 : 0] int bounds_type_sigma;

  // prior distribution specification of the parameteres
  array[estimate_mu     == 0 ? 1 : 0] real fixed_mu;
  array[estimate_sigma  == 0 ? 1 : 0] real fixed_sigma;
  vector[estimate_mu    == 1 ? 3 : 0] prior_parameters_mu;
  vector[estimate_sigma == 1 ? 3 : 0] prior_parameters_sigma;
  int prior_type_mu;
  int prior_type_sigma;
}
parameters{
  array[estimate_mu]    real<lower = coefs_lb(bounds_type_mu, bounds_mu),       upper = coefs_ub(bounds_type_mu, bounds_mu)>       mu_est;
  array[estimate_sigma] real<lower = coefs_lb(bounds_type_sigma, bounds_sigma), upper = coefs_ub(bounds_type_sigma, bounds_sigma)> sigma_est;
}
transformed parameters {
  real mu;
  real sigma;

  if(estimate_mu == 1){
    mu = mu_est[1];
  }else{
    mu = fixed_mu[1];
  }
  if(estimate_sigma == 1){
    sigma = sigma_est[1];
  }else{
    sigma = fixed_sigma[1];
  }
}
model {
  // priors for mu and sigma
  if(estimate_mu    == 1) target += set_prior(mu_est[1],     prior_type_mu,     prior_parameters_mu,     bounds_type_mu,     bounds_mu);
  if(estimate_sigma == 1) target += set_prior(sigma_est[1],  prior_type_sigma,  prior_parameters_sigma,  bounds_type_sigma,  bounds_sigma);

  // likelihood of the data
  if(is_ss == 0){
    target += normal_lpdf(x | mu, sigma);
  }else{
    target += -N / 2.0 * log(2 * pi() * pow(sigma, 2)) - 1 / (2 * pow(sigma, 2)) * ((N - 1) * pow(ss_sd, 2) + N * (ss_mean - mu)^2);
  }
}
