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

  // range of the parameters
  vector[2] bounds_mu;
  vector[2] bounds_sigma;
  array[2] int bounds_type_mu;
  array[2] int bounds_type_sigma;

  // prior distribution specification of the parameteres
  vector[3] prior_parameters_mu;
  vector[3] prior_parameters_sigma;
  int prior_type_mu;
  int prior_type_sigma;
}
parameters{
  real<lower = coefs_lb(bounds_type_mu, bounds_mu),       upper = coefs_ub(bounds_type_mu, bounds_mu)>       mu;
  real<lower = coefs_lb(bounds_type_sigma, bounds_sigma), upper = coefs_ub(bounds_type_sigma, bounds_sigma)> sigma;
}
model {
  // priors for mu and sigma2
  target += set_prior(mu,     prior_type_mu,     prior_parameters_mu,     bounds_type_mu,     bounds_mu);
  target += set_prior(sigma,  prior_type_sigma,  prior_parameters_sigma,  bounds_type_sigma,  bounds_sigma);

  // likelihood of the data
  if(is_ss == 0){
    target += normal_lpdf(x | mu, sigma);
  }else{
    target += -N / 2.0 * log(2 * pi() * pow(sigma, 2)) - 1 / (2 * pow(sigma, 2)) * ((N - 1) * pow(ss_sd, 2) + N * (ss_mean - mu)^2);
  }
}
