#include /include/common_functions.stan

data {
  // sample sizes
  int<lower=0> N;
  // individual observations
  vector[N] x;

  // range of the parameters
  vector[2] bounds_mu;
  vector[2] bounds_sigma;
  vector[2] bounds_nu;
  array[2] int bounds_type_mu;
  array[2] int bounds_type_sigma;
  array[2] int bounds_type_nu;

  // prior distribution specification of the parameteres
  vector[3] prior_parameters_mu;
  vector[3] prior_parameters_sigma;
  vector[3] prior_parameters_nu;
  int prior_type_mu;
  int prior_type_sigma;
  int prior_type_nu;
}
parameters{
  real<lower = coefs_lb(bounds_type_mu, bounds_mu),       upper = coefs_ub(bounds_type_mu, bounds_mu)>       mu;
  real<lower = coefs_lb(bounds_type_sigma, bounds_sigma), upper = coefs_ub(bounds_type_sigma, bounds_sigma)> sigma;
  real<lower = coefs_lb(bounds_type_nu, bounds_nu),       upper = coefs_ub(bounds_type_nu, bounds_nu)>       nu_p; // nu = nu_p + 2
}
transformed parameters {
  real nu    = nu_p + 2;
  real scale = sigma / sqrt(nu / (nu - 2.0));
}
model {
  // priors for mu and sigma2
  target += set_prior(mu,     prior_type_mu,     prior_parameters_mu,     bounds_type_mu,     bounds_mu);
  target += set_prior(sigma,  prior_type_sigma,  prior_parameters_sigma , bounds_type_sigma,  bounds_sigma);
  target += set_prior(nu_p,   prior_type_nu,     prior_parameters_nu,     bounds_type_nu,     bounds_nu);

  // likelihood of the data
  target += student_t_lpdf(x | nu, mu, scale);
}
