#include /include/common_functions.stan

data {
  // sample sizes
  int<lower=0> N;
  // individual observations
  vector[N] x;

  // model type (1 if the parameter is estimated, 0 if parameter is fixed)
  int estimate_mu;
  int estimate_sigma;
  int estimate_nu;

  // range of the parameters
  vector[estimate_mu    == 1 ? 2 : 0] bounds_mu;
  vector[estimate_sigma == 1 ? 2 : 0] bounds_sigma;
  vector[estimate_nu    == 1 ? 2 : 0] bounds_nu;
  array[estimate_mu    == 1 ? 2 : 0] int bounds_type_mu;
  array[estimate_sigma == 1 ? 2 : 0] int bounds_type_sigma;
  array[estimate_nu    == 1 ? 2 : 0] int bounds_type_nu;

  // fixed values for non-estimated parameters
  array[estimate_mu    == 0 ? 1 : 0] real fixed_mu;
  array[estimate_sigma == 0 ? 1 : 0] real fixed_sigma;
  array[estimate_nu    == 0 ? 1 : 0] real fixed_nu;

  // prior distribution specification of the parameters
  vector[estimate_mu    == 1 ? 3 : 0] prior_parameters_mu;
  vector[estimate_sigma == 1 ? 3 : 0] prior_parameters_sigma;
  vector[estimate_nu    == 1 ? 3 : 0] prior_parameters_nu;
  int prior_type_mu;
  int prior_type_sigma;
  int prior_type_nu;

  int sample_priors;
}
parameters{
  array[estimate_mu]    real<lower = coefs_lb(bounds_type_mu, bounds_mu),       upper = coefs_ub(bounds_type_mu, bounds_mu)>       mu_est;
  array[estimate_sigma] real<lower = coefs_lb(bounds_type_sigma, bounds_sigma), upper = coefs_ub(bounds_type_sigma, bounds_sigma)> sigma_est;
  array[estimate_nu]    real<lower = coefs_lb(bounds_type_nu, bounds_nu),       upper = coefs_ub(bounds_type_nu, bounds_nu)>       nu_est; // df = nu_est + 2
}
transformed parameters {
  real mu;
  real scale;
  real nu;

  if (estimate_nu == 1) {
    nu = nu_est[1] + 2;
  } else {
    nu = fixed_nu[1];
  }
  if (estimate_mu == 1) {
    mu = mu_est[1];
  } else {
    mu = fixed_mu[1];
  }
  if (estimate_sigma == 1) {
    scale = sigma_est[1] / sqrt(nu / (nu - 2.0));
  } else {
    scale = fixed_sigma[1] / sqrt(nu / (nu - 2.0));
  }
}
model {
  // priors for mu and sigma2
  if (estimate_mu    == 1) target += set_prior(mu_est[1],     prior_type_mu,     prior_parameters_mu,     bounds_type_mu,     bounds_mu);
  if (estimate_sigma == 1) target += set_prior(sigma_est[1],  prior_type_sigma,  prior_parameters_sigma,  bounds_type_sigma,  bounds_sigma);
  if (estimate_nu    == 1) target += set_prior(nu_est[1],     prior_type_nu,     prior_parameters_nu,     bounds_type_nu,     bounds_nu);

  // likelihood of the data
  if (sample_priors == 0) target += student_t_lpdf(x | nu, mu, scale);
}
