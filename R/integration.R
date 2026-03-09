# ==============================================================================
# Numerical Integration Backend for Bayesian Capability Analysis
# ==============================================================================

# Helper function for log-space difference of exponentials
log_diff_exp <- function(x, y) {
  ifelse(x <= y, -Inf, x + log1p(-exp(y - x)))
}

# ==============================================================================
# Prior Classes
# ==============================================================================

#' Create a conjugate prior for Normal-InverseGamma model
#' @param mu0 Prior mean location
#' @param k0 Prior precision multiplier (0 = noninformative)
#' @param alpha0 Prior shape for sigma^2
#' @param beta0 Prior rate for sigma^2
#' @return PriorConjugate object
#' @export
create_prior_conjugate <- function(mu0 = 0, k0 = 0, alpha0 = -0.5, beta0 = 0) {
  structure(list(mu0 = mu0, k0 = k0, alpha0 = alpha0, beta0 = beta0),
            class = "PriorConjugate")
}

#' Create a generic prior with custom log-density function
#' @param log_dens_fn Function(mu, sigma) returning log prior density
#' @param bayestools_priors Optional list with original BayesTools prior objects (mu and sigma)
#' @return PriorGeneric object
#' @keywords internal
create_prior_generic <- function(log_dens_fn, bayestools_priors = NULL) {
  structure(list(log_dens = log_dens_fn, bayestools_priors = bayestools_priors),
            class = "PriorGeneric")
}

#' Create a semi-conjugate prior (conjugate mu, non-conjugate sigma)
#' @param mu0 Prior mean location
#' @param k0 Prior precision multiplier (0 = noninformative)
#' @param log_dens_sigma Function(sigma) returning log prior density for sigma
#' @param bayestools_priors Optional list with original BayesTools prior objects
#' @return PriorSemiConjugateMu object
#' @keywords internal
create_prior_semi_mu <- function(mu0, k0, log_dens_sigma, bayestools_priors = NULL) {
  structure(list(mu0 = mu0, k0 = k0, log_dens_sigma = log_dens_sigma,
                 bayestools_priors = bayestools_priors),
            class = "PriorSemiConjugateMu")
}

#' Create a semi-conjugate prior (non-conjugate mu, conjugate sigma)
#' @param alpha0 Prior shape for sigma^2
#' @param beta0 Prior rate for sigma^2
#' @param log_dens_mu Function(mu) returning log prior density for mu
#' @param bayestools_priors Optional list with original BayesTools prior objects
#' @return PriorSemiConjugateSigma object
#' @keywords internal
create_prior_semi_sigma <- function(alpha0, beta0, log_dens_mu, bayestools_priors = NULL) {
  structure(list(alpha0 = alpha0, beta0 = beta0, log_dens_mu = log_dens_mu,
                 bayestools_priors = bayestools_priors),
            class = "PriorSemiConjugateSigma")
}

# ==============================================================================
# Sufficient Statistics & Posterior Update Helpers
# ==============================================================================

.extract_suff_stats <- function(data, cached_state) {
  if (!is.null(cached_state)) {
    list(n = cached_state$n, x_bar = cached_state$x_bar, SS = cached_state$sse)
  } else {
    n <- length(data)
    if (n > 0) {
      x_bar <- mean(data)
      list(n = n, x_bar = x_bar, SS = sum((data - x_bar)^2))
    } else {
      list(n = 0L, x_bar = 0, SS = 0)
    }
  }
}

.nig_posterior <- function(prior, n, x_bar, SS) {
  if (n > 0) {
    k_n     <- prior$k0 + n
    mu_n    <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
    alpha_n <- prior$alpha0 + n / 2
    beta_n  <- prior$beta0 + 0.5 * SS + (prior$k0 * n * (x_bar - prior$mu0)^2) / (2 * k_n)
  } else {
    k_n     <- prior$k0
    mu_n    <- prior$mu0
    alpha_n <- prior$alpha0
    beta_n  <- prior$beta0
  }
  list(k_n = k_n, mu_n = mu_n, alpha_n = alpha_n, beta_n = beta_n)
}

.semi_mu_posterior <- function(prior, n, x_bar, sse) {
  k_n   <- prior$k0 + n
  mu_n  <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
  sse_n <- sse + prior$k0 * n * (x_bar - prior$mu0)^2 / k_n
  list(k_n = k_n, mu_n = mu_n, sse_n = sse_n)
}

# ==============================================================================
# Metric Constraints
# ==============================================================================

#' Get constraint functions for a capability metric
#' @param metric One of "Cp", "Cpk", "Cpm", "Cpc", "Cpu", "Cpl"
#' @param c Threshold value
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param target Target value
#' @return List with s_max_fn and mu_b_fn_vec (vectorized)
#' @keywords internal
get_metric_constraints <- function(metric, c, LSL, USL, target) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  list(
    s_max_fn = function() {
      if (metric %in% c("Cp", "Cpk", "Cpm", "Cpc")) return(tol / (6 * c))
      return(Inf)
    },
    # Vectorized version: takes vector of sigma, returns list(lower, upper)
    mu_b_fn_vec = function(s) {
      n <- length(s)
      if (metric == "Cp") {
        return(list(lower = rep(-Inf, n), upper = rep(Inf, n)))
      }
      if (metric == "Cpk") {
        return(list(lower = LSL + 3 * c * s, upper = USL - 3 * c * s))
      }
      if (metric == "Cpu") {
        return(list(lower = rep(-Inf, n), upper = USL - 3 * c * s))
      }
      if (metric == "Cpl") {
        return(list(lower = LSL + 3 * c * s, upper = rep(Inf, n)))
      }
      if (metric %in% c("Cpm", "Cpc")) {
        T_val <- if (metric == "Cpc") mid else target
        R <- tol / (6 * c)
        w <- sqrt(pmax(0, R^2 - s^2))
        lower <- T_val - w
        upper <- T_val + w
        invalid <- s >= R
        lower[invalid] <- 0
        upper[invalid] <- -1
        return(list(lower = lower, upper = upper))
      }
      stop("Unknown metric: ", metric)
    }
  )
}

# ==============================================================================
# Metric Value Computation
# ==============================================================================

#' Compute capability metric value for given (mu, sigma) pairs
#' @param mu Mean value(s) (vectorized)
#' @param sigma Standard deviation value(s) (vectorized)
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param target Target value (required for Cpm)
#' @param metric One of "Cp", "Cpk", "Cpu", "Cpl", "Cpm", "Cpc"
#' @return Metric value(s) (same length as mu/sigma)
#' @keywords internal
compute_metric_value <- function(mu, sigma, LSL, USL, target, metric) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  switch(metric,
    "Cp"  = tol / (6 * sigma),
    "Cpu" = (USL - mu) / (3 * sigma),
    "Cpl" = (mu - LSL) / (3 * sigma),
    "Cpk" = pmin((USL - mu) / (3 * sigma), (mu - LSL) / (3 * sigma)),
    "Cpm" = tol / (6 * sqrt(sigma^2 + (mu - target)^2)),
    "Cpc" = tol / (6 * sqrt(sigma^2 + (mu - mid)^2)),
    stop("Unknown metric: ", metric)
  )
}

# ==============================================================================
# Posterior Moments of Capability Metrics
# ==============================================================================

#' Compute posterior mean and SD of a capability metric
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param use_analytic Use analytic formulas (TRUE) or numerical fallback (FALSE)
#' @param cached_state Pre-computed state for PriorGeneric
#' @return List with mean and sd of the posterior distribution of the metric
#' @keywords internal
compute_metric_moments <- function(data, LSL, USL, prior, metric = "Cpk",
                                   target = NULL, use_analytic = TRUE,
                                   cached_state = NULL) {
  UseMethod("compute_metric_moments", prior)
}

#' Compute grid bounds for a metric using sigma quantiles (prior-only case)
#'
#' When sampling from priors only (no data), we can compute grid bounds directly
#' from the sigma prior quantiles using BayesTools::quant, avoiding expensive
#' 2D numerical integration for moment computation.
#'
#' @param prior_sigma BayesTools prior object for sigma
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param p_low Lower quantile probability (default 0.001)
#' @param p_high Upper quantile probability (default 0.999)
#' @return List with x_start and x_end for the grid
#' @keywords internal
.compute_metric_grid_bounds_from_quantiles <- function(prior_sigma, LSL, USL,
                                                        metric, target = NULL,
                                                        p_low = 0.001,
                                                        p_high = 0.999) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  # Get sigma quantiles using BayesTools::quant
  sigma_low <- BayesTools::quant(prior_sigma, p_low)
  sigma_high <- BayesTools::quant(prior_sigma, p_high)

  # For metrics inversely related to sigma (Cp, Cpk, Cpu, Cpl, Cpm, Cpc),

# metric_high corresponds to sigma_low and vice versa
  metric_at_sigma_low <- switch(metric,
    "Cp" = tol / (6 * sigma_low),
    "Cpk" = tol / (6 * sigma_low),  # Upper bound (assumes mu = mid)
    "Cpu" = (USL - mid) / (3 * sigma_low),
    "Cpl" = (mid - LSL) / (3 * sigma_low),
    "Cpm" = tol / (6 * sigma_low),  # Upper bound (assumes mu = target)
    "Cpc" = tol / (6 * sigma_low),  # Upper bound approximation
    tol / (6 * sigma_low)  # Default
  )

  metric_at_sigma_high <- switch(metric,
    "Cp" = tol / (6 * sigma_high),
    "Cpk" = tol / (6 * sigma_high),
    "Cpu" = (USL - mid) / (3 * sigma_high),
    "Cpl" = (mid - LSL) / (3 * sigma_high),
    "Cpm" = tol / (6 * sigma_high),
    "Cpc" = tol / (6 * sigma_high),
    tol / (6 * sigma_high)
  )

  x_start <- max(0, metric_at_sigma_high)
  x_end <- metric_at_sigma_low

  # Ensure valid range
  if (x_end <= x_start || !is.finite(x_end)) {
    x_end <- max(x_start + 3, 10)
  }

  list(x_start = x_start, x_end = x_end)
}

#' @export
compute_metric_moments.PriorConjugate <- function(data, LSL, USL, prior,
                                                   metric = "Cpk", target = NULL,
                                                   use_analytic = TRUE,
                                                   cached_state = NULL) {
  ss  <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  if (!use_analytic) {
    # Numerical fallback: 2D integration
    return(.compute_moments_numerical_conjugate(
      mu_n, k_n, alpha_n, beta_n, LSL, USL, target, metric
    ))
  }

  # Analytic computation
  # E[1/sigma] for Inverse-Gamma(alpha, beta): sqrt(beta) * Gamma(alpha-0.5) / Gamma(alpha) / sqrt(2)
  # Using the chi-square parameterization: sigma^2 = 2*beta/y where y ~ chi^2(2*alpha)
  # E[1/sigma] = E[sqrt(y/(2*beta))] = E[sqrt(y)] / sqrt(2*beta)
  # E[sqrt(y)] for y ~ chi^2(df) = sqrt(2) * Gamma((df+1)/2) / Gamma(df/2)
  df_p <- 2 * alpha_n
  E_inv_sigma <- sqrt(2) * exp(lgamma((df_p + 1) / 2) - lgamma(df_p / 2)) / sqrt(2 * beta_n)
  E_inv_sigma2 <- df_p / (2 * beta_n)  # E[y/(2*beta)] = E[y]/(2*beta) = df/(2*beta)

  if (metric == "Cp") {
    # Cp = tol / (6*sigma)
    E1 <- (tol / 6) * E_inv_sigma
    # E[Cp^2] = (tol/6)^2 * E[1/sigma^2]
    E2 <- (tol / 6)^2 * E_inv_sigma2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpu") {
    # Cpu = (USL - mu) / (3*sigma)
    # E[Cpu] = E[(USL - mu)/sigma] / 3
    # mu|sigma ~ N(mu_n, sigma^2/k_n), so E[mu|sigma] = mu_n
    # E[(USL - mu)/sigma] = (USL - mu_n) * E[1/sigma]
    E1 <- (USL - mu_n) / 3 * E_inv_sigma

    # E[Cpu^2] = E[(USL - mu)^2 / sigma^2] / 9
    # (USL - mu)^2 = (USL - mu_n)^2 - 2*(USL - mu_n)*(mu - mu_n) + (mu - mu_n)^2
    # E[(mu - mu_n)^2 | sigma] = sigma^2/k_n
    # E[(USL - mu)^2 / sigma^2] = (USL - mu_n)^2 * E[1/sigma^2] + E[1/k_n] = (USL - mu_n)^2 * E[1/sigma^2] + 1/k_n
    E2 <- ((USL - mu_n)^2 * E_inv_sigma2 + 1 / k_n) / 9
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpl") {
    # Cpl = (mu - LSL) / (3*sigma)
    E1 <- (mu_n - LSL) / 3 * E_inv_sigma
    E2 <- ((mu_n - LSL)^2 * E_inv_sigma2 + 1 / k_n) / 9
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpk") {
    # Cpk = min(Cpu, Cpl) = (tol/2 - |mu - mid|) / (3*sigma)
    # E[Cpk] = (tol/2)/3 * E[1/sigma] - (1/3) * E[|mu - mid|/sigma]

    # E[|mu - mid|/sigma] requires integrating over sigma
    # (mu - mid)|sigma ~ N(mu_n - mid, sigma^2/k_n)
    # |X| where X ~ N(delta, tau^2): E[|X|] = tau*sqrt(2/pi)*exp(-delta^2/(2*tau^2)) + delta*(1 - 2*Phi(-delta/tau))
    delta <- mu_n - mid
    # tau = sigma/sqrt(k_n), so we need E[f(sigma)] over marginal of sigma

    # Integrate E[|mu - mid| | sigma] * p(sigma) over sigma
    # Using y = 2*beta/sigma^2 ~ chi^2(df_p), sigma = sqrt(2*beta/y)
    integrand_abs <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      # Mean of |X| for X ~ N(delta, tau^2)
      abs_mean <- tau * sqrt(2 / pi) * exp(-delta^2 / (2 * tau^2)) +
                  delta * (1 - 2 * stats::pnorm(-delta / tau))
      abs_mean / sigma * stats::dchisq(y, df_p)
    }
    E_abs_div_sigma <- stats::integrate(integrand_abs, 0, Inf, rel.tol = 1e-6)$value

    E1 <- (tol / 2) / 3 * E_inv_sigma - E_abs_div_sigma / 3

    # E[Cpk^2] - more complex, use 1D numerical integration
    integrand_sq <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      # E[(tol/2 - |mu - mid|)^2 | sigma] / (9*sigma^2)
      # = E[(tol/2)^2 - tol*|mu-mid| + |mu-mid|^2 | sigma] / (9*sigma^2)
      abs_mean <- tau * sqrt(2 / pi) * exp(-delta^2 / (2 * tau^2)) +
                  delta * (1 - 2 * stats::pnorm(-delta / tau))
      # E[|X|^2] = E[X^2] = delta^2 + tau^2
      abs2_mean <- delta^2 + tau^2
      cpk2_given_sigma <- ((tol / 2)^2 - tol * abs_mean + abs2_mean) / (9 * sigma^2)
      cpk2_given_sigma * stats::dchisq(y, df_p)
    }
    E2 <- stats::integrate(integrand_sq, 0, Inf, rel.tol = 1e-6)$value

    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric %in% c("Cpm", "Cpc")) {
    # Cpm = tol / (6 * sqrt(sigma^2 + (mu - T)^2))
    # No closed form - use 1D numerical integration over sigma (via chi-square)
    T_val <- if (metric == "Cpc") mid else target
    delta_T <- mu_n - T_val

    # Use bounded integration to avoid divergence
    # Chi-square y has most mass near df_p, so integrate from small epsilon to large upper bound
    y_upper <- max(100, df_p * 10)

    z_pts <- c(-3, -2, -1, 0, 1, 2, 3)
    w_pts <- stats::dnorm(z_pts)
    w_pts <- w_pts / sum(w_pts)

    integrand_cpm <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      x_pts <- tcrossprod(tau, z_pts) + delta_T
      sigma2 <- sigma^2
      inv_sqrt_vals <- 1 / sqrt(x_pts^2 + sigma2)
      inner <- drop(inv_sqrt_vals %*% w_pts)
      inner * stats::dchisq(y, df_p)
    }
    E_inv_sqrt <- stats::integrate(integrand_cpm, 1e-6, y_upper,
                            rel.tol = 1e-4, subdivisions = 200)$value
    E1 <- (tol / 6) * E_inv_sqrt

    integrand_cpm2 <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      x_pts <- tcrossprod(tau, z_pts) + delta_T
      sigma2 <- sigma^2
      inv_vals <- 1 / (x_pts^2 + sigma2)
      inner <- drop(inv_vals %*% w_pts)
      inner * stats::dchisq(y, df_p)
    }
    E_inv <- stats::integrate(integrand_cpm2, 1e-6, y_upper,
                       rel.tol = 1e-4, subdivisions = 200)$value
    E2 <- (tol / 6)^2 * E_inv

    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }


  stop("Unknown metric: ", metric)
}

#' Numerical fallback for conjugate prior moments (2D integration)
#' @keywords internal
.compute_moments_numerical_conjugate <- function(mu_n, k_n, alpha_n, beta_n,
                                                  LSL, USL, target, metric) {
  df_p <- 2 * alpha_n
  h_max <- stats::dchisq(max(df_p - 2, 1e-6), df_p, log = TRUE)

  # Vectorized log posterior computation
  log_post_vec <- function(mu, sigma) {
    log_post <- rep(-Inf, length(mu))
    valid <- sigma > 0
    if (any(valid)) {
      y <- 2 * beta_n / sigma[valid]^2
      # Posterior: N(mu | mu_n, sigma^2/k_n) * chi^2(y | df_p) * |dy/dsigma|
      # |dy/dsigma| = 4*beta_n / sigma^3
      log_post[valid] <- stats::dnorm(mu[valid], mu_n, sigma[valid] / sqrt(k_n), log = TRUE) +
                         stats::dchisq(y, df_p, log = TRUE) + log(4 * beta_n) - 3 * log(sigma[valid])
    }
    log_post
  }

  # Fully vectorized integrand
  integrand <- function(x, power = 1) {
    mu <- x[1, ]
    sigma <- x[2, ]

    log_post <- log_post_vec(mu, sigma)
    metric_val <- compute_metric_value(mu, sigma, LSL, USL, target, metric)

    result <- exp(log_post - h_max) * (metric_val^power)
    result[!is.finite(result)] <- 0
    matrix(result, nrow = 1)
  }

  # Vectorized normalizing constant integrand
  integrand_norm <- function(x) {
    mu <- x[1, ]
    sigma <- x[2, ]

    log_post <- log_post_vec(mu, sigma)
    result <- exp(log_post - h_max)
    result[!is.finite(result)] <- 0
    matrix(result, nrow = 1)
  }

  # Integration bounds
  sigma_mode <- sqrt(2 * beta_n / max(df_p - 2, 1))
  sigma_max <- sigma_mode * 5
  mu_lower <- mu_n - 5 * sigma_max / sqrt(k_n)
  mu_upper <- mu_n + 5 * sigma_max / sqrt(k_n)

  Z <- cubature::pcubature(
    integrand_norm,
    lowerLimit = c(mu_lower, 1e-10),
    upperLimit = c(mu_upper, sigma_max),
    tol = 1e-4, vectorInterface = TRUE
  )$integral

  E1 <- cubature::pcubature(
    function(x) integrand(x, power = 1),
    lowerLimit = c(mu_lower, 1e-10),
    upperLimit = c(mu_upper, sigma_max),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  E2 <- cubature::pcubature(
    function(x) integrand(x, power = 2),
    lowerLimit = c(mu_lower, 1e-10),
    upperLimit = c(mu_upper, sigma_max),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}


#' @export
compute_metric_moments.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                                         metric = "Cpk", target = NULL,
                                                         use_analytic = TRUE,
                                                         cached_state = NULL) {
  ss <- .extract_suff_stats(data, NULL)
  n <- ss$n; x_bar <- ss$x_bar; sse <- ss$SS
  smp <- .semi_mu_posterior(prior, n, x_bar, sse)
  k_n <- smp$k_n; mu_n <- smp$mu_n; sse_n <- smp$sse_n
  sd_data <- sqrt(sse_n / max(n, 1))
  n_eff <- if (prior$k0 == 0) n - 1 else n

  spec_tol <- USL - LSL
  M <- (LSL + USL) / 2
  if (is.null(target)) target <- M

  log_sw <- function(sigma) {
    -n_eff * log(sigma) - sse_n / (2 * sigma^2) + prior$log_dens_sigma(sigma)
  }

  # Integrate over t = log(sigma) for numerical stability (peaked integrand).
  t_mode <- stats::optimize(function(t) -log_sw(exp(t)),
                            c(log(max(sd_data * 0.01, 1e-10)), log(max(sd_data * 100, 10))))$minimum
  h_max  <- log_sw(exp(t_mode))
  t_lo   <- t_mode - 30
  t_hi   <- t_mode + 30

  base_w <- function(t) {
    sigma <- exp(t)
    exp(log_sw(sigma) - h_max + t)
  }

  Z <- stats::integrate(base_w, t_lo, t_hi, rel.tol = 1e-5)$value

  if (metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    # Analytical E[metric|sigma] and E[metric^2|sigma] under mu|sigma ~ N(mu_n, sigma^2/k_n)
    inner_moments <- function(sigma) {
      sd_mu <- sigma / sqrt(k_n)
      a <- mu_n - LSL; cc <- USL - mu_n; b <- sd_mu
      inv3s <- 1 / (3 * sigma)
      switch(metric,
        "Cp"  = list(E1 = spec_tol * inv3s / 2, E2 = (spec_tol * inv3s / 2)^2),
        "Cpu" = list(E1 = cc * inv3s, E2 = (cc^2 + b^2) * inv3s^2),
        "Cpl" = list(E1 = a * inv3s,  E2 = (a^2 + b^2) * inv3s^2),
        "Cpk" = {
          zs  <- (M - mu_n) / sd_mu
          Phi <- stats::pnorm(zs); phi <- stats::dnorm(zs)
          list(
            E1 = (a * Phi + cc * (1 - Phi) - 2 * b * phi) * inv3s,
            E2 = ((a^2 + b^2) * Phi + (cc^2 + b^2) * (1 - Phi) -
                    2 * b * (a + cc) * phi) * inv3s^2
          )
        }
      )
    }
    E1 <- stats::integrate(function(t) {
      m <- inner_moments(exp(t)); base_w(t) * m$E1
    }, t_lo, t_hi, rel.tol = 1e-5)$value / Z
    E2 <- stats::integrate(function(t) {
      m <- inner_moments(exp(t)); base_w(t) * m$E2
    }, t_lo, t_hi, rel.tol = 1e-5)$value / Z
  } else {
    # GH quadrature for smooth metrics (Cpm, Cpc).
    # 9-point probabilist nodes (physicist × sqrt(2)) for N(mu_n, sd_mu^2).
    gh_nodes <- sqrt(2) * c(-3.190994, -2.266581, -1.468554, -0.723551, 0,
                              0.723551,  1.468554,  2.266581,  3.190994)
    gh_wts   <- c(3.961e-05, 4.944e-03, 8.847e-02, 4.326e-01, 7.202e-01,
                  4.326e-01, 8.847e-02, 4.944e-03, 3.961e-05)
    gh_wts   <- gh_wts / sum(gh_wts)
    n_gh     <- length(gh_nodes)

    gh_integrand <- function(t, power) {
      sigma <- exp(t); w <- base_w(t)
      ok <- is.finite(w) & w > 0
      result <- numeric(length(t))
      if (!any(ok)) return(result)
      sd_mu <- sigma[ok] / sqrt(k_n)
      mu_mat <- tcrossprod(sd_mu, gh_nodes) + mu_n
      sigma_rep <- rep(sigma[ok], each = n_gh)
      m_vals <- compute_metric_value(as.vector(t(mu_mat)), sigma_rep, LSL, USL, target, metric)
      m_mat <- matrix(m_vals, ncol = n_gh, byrow = TRUE)
      result[ok] <- w[ok] * drop(m_mat^power %*% gh_wts)
      result
    }
    E1 <- stats::integrate(function(t) gh_integrand(t, 1),
                           t_lo, t_hi, rel.tol = 1e-5, subdivisions = 200)$value / Z
    E2 <- stats::integrate(function(t) gh_integrand(t, 2),
                           t_lo, t_hi, rel.tol = 1e-5, subdivisions = 200)$value / Z
  }

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

#' @export
compute_metric_moments.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                            metric = "Cpk", target = NULL,
                                                            use_analytic = TRUE,
                                                            cached_state = NULL) {
  # Case 3: Non-conjugate mu, conjugate sigma (InvGamma)
  # Use 2D integration: outer over mu, inner over sigma via chi-square quadrature.
  # Previous version used E[sigma|mu] as plug-in (Jensen bias for 1/sigma metrics).

  n <- length(data)
  x_bar <- if (n > 0) mean(data) else 0
  sse <- if (n > 0) sum((data - x_bar)^2) else 0

  alpha0 <- prior$alpha0
  beta0 <- prior$beta0
  alpha_n <- alpha0 + n / 2
  df_p <- 2 * alpha_n

  sd_data <- sqrt(sse / max(n - 1, 1))
  alpha_0_times_logbeta0 <- if (prior$beta0 == 0) 0 else alpha0 * log(beta0)

  # 7-point Gauss-Hermite nodes/weights for inner sigma integration
  gh7_nodes <- c(-2.651961, -1.673552, -0.816288, 0, 0.816288, 1.673552, 2.651961)
  gh7_wts   <- c(0.0009718, 0.054536, 0.42560, 0.81026, 0.42560, 0.054536, 0.0009718)
  gh7_wts   <- gh7_wts / sum(gh7_wts)
  n_gh      <- length(gh7_nodes)

  .inner_sigma_expect <- function(mu_vec, beta_n_vec, power) {
    # For each mu, sigma^2|mu ~ InvGamma(alpha_n, beta_n).
    # Transform: y = 2*beta_n/sigma^2 ~ chi^2(df_p).
    # Use GH-like quadrature on log(y) centered at mode.
    n_mu <- length(mu_vec)
    result <- numeric(n_mu)
    y_mode <- max(df_p - 2, 0.5)
    log_y_mode <- log(y_mode)
    # SD of log(chi^2) ~ sqrt(2/df_p) for large df_p, use broader spread for small df
    log_y_sd <- sqrt(2 / max(df_p, 1))

    log_y_pts <- log_y_mode + log_y_sd * gh7_nodes  # n_gh points
    y_pts <- exp(log_y_pts)                          # length n_gh

    # chi^2 density at these y points (constant across mu)
    log_dchisq <- stats::dchisq(y_pts, df_p, log = TRUE)
    # Jacobian: dy = y * d(log_y), so weight includes y
    log_w <- log_dchisq + log_y_pts + log(log_y_sd)  # absorb constants
    w_base <- exp(log_w - max(log_w))
    w_base <- w_base / sum(w_base)

    for (i in seq_len(n_mu)) {
      sigma_pts <- sqrt(2 * beta_n_vec[i] / y_pts)
      m_vals <- compute_metric_value(
        rep(mu_vec[i], n_gh), sigma_pts, LSL, USL, target, metric
      )
      result[i] <- sum(w_base * m_vals^power)
    }
    result
  }

  integrand_Ek <- function(mu, power) {
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    log_marginal <- lgamma(alpha_n) - lgamma(alpha0) + alpha_0_times_logbeta0 - alpha_n * log(beta_n)
    log_prior <- prior$log_dens_mu(mu)
    weight <- exp(log_marginal + log_prior)
    ok <- is.finite(weight) & weight > 0
    result <- numeric(length(mu))
    if (!any(ok)) return(result)
    inner <- .inner_sigma_expect(mu[ok], beta_n[ok], power)
    result[ok] <- weight[ok] * inner
    result
  }

  mu_low <- x_bar - 10 * sd_data
  mu_high <- x_bar + 10 * sd_data

  Z <- stats::integrate(function(mu) {
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    log_marginal <- lgamma(alpha_n) - lgamma(alpha0) + alpha_0_times_logbeta0 - alpha_n * log(beta_n)
    exp(log_marginal + prior$log_dens_mu(mu))
  }, mu_low, mu_high, rel.tol = 1e-5)$value

  E1 <- stats::integrate(function(mu) integrand_Ek(mu, 1), mu_low, mu_high,
                          rel.tol = 1e-4, subdivisions = 200)$value / Z
  E2 <- stats::integrate(function(mu) integrand_Ek(mu, 2), mu_low, mu_high,
                          rel.tol = 1e-4, subdivisions = 200)$value / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

#' @export
compute_metric_moments.PriorGeneric <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 use_analytic = TRUE,
                                                 cached_state = NULL) {
  # For generic priors, always use numerical 2D integration
  # (use_analytic is ignored - no closed form available)

  if (is.null(cached_state)) {
    cached_state <- precompute_generic_state(data, prior)
  }

  log_post_vec <- cached_state$log_post_vec
  h_max <- cached_state$h_max
  uni_s <- cached_state$uni_s
  map_mu <- cached_state$map_mu
  Z <- cached_state$Z

  # Fully vectorized integrand
  integrand <- function(x, power = 1) {
    mu <- x[1, ]
    sigma <- x[2, ]

    # Vectorized metric and log posterior computation
    metric_val <- compute_metric_value(mu, sigma, LSL, USL, target, metric)
    log_p <- log_post_vec(mu, sigma)

    result <- exp(log_p - h_max) * (metric_val^power)
    result[!is.finite(result)] <- 0
    matrix(result, nrow = 1)
  }

  mu_lower <- map_mu - 5 * uni_s
  mu_upper <- map_mu + 5 * uni_s

  # xgrid <- as.matrix(expand.grid(
  #   mu = seq(mu_lower, mu_upper, length.out = 5),
  #   sigma = seq(1e-10, uni_s, length.out = 5)
  # ))
  # log_post_vals <- integrand(t(xgrid))
  # df <- data.frame(
  #   mu = xgrid[, 1],
  #   sigma = xgrid[, 2],
  #   log_post = as.numeric(log_post_vals)
  # )
  # ggplot2::ggplot(data = df, ggplot2::aes(x = mu, y = sigma, fill = log(log_post))) +
  #   ggplot2::geom_tile() +
  #   ggplot2::scale_fill_viridis_c()

  # Both the integral and Z use the same h_max normalization, so they cancel out
  E1 <- cubature::pcubature(
    function(x) integrand(x, power = 1),
    lowerLimit = c(mu_lower, 1e-10),
    upperLimit = c(mu_upper, uni_s),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  E2 <- cubature::pcubature(
    function(x) integrand(x, power = 2),
    lowerLimit = c(mu_lower, 1e-10),
    upperLimit = c(mu_upper, uni_s),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}


# ==============================================================================
# Solver Factory: Returns function P(Index > c)
# ==============================================================================

#' Create solver function for P(Index > c)
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param cached_state Pre-computed state for PriorGeneric (optional)
#' @return Function that takes threshold c and returns P(Index > c)
#' @keywords internal
make_solver <- function(data, LSL, USL, prior, metric = "Cpk", target = NULL, ...) {
  UseMethod("make_solver", prior)
}

#' @export
make_solver.PriorConjugate <- function(data, LSL, USL, prior, metric = "Cpk",
                                        target = NULL, cached_state = NULL, ...) {
  ss   <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n
  df_p <- 2 * alpha_n

  # Pre-compute global h_max (chi-square mode density for numerical stability)
  y_mode <- max(df_p - 2, 1e-6)
  h_max_global <- stats::dchisq(y_mode, df_p, log = TRUE)

  function(c) {
    if (c <= 0) return(1.0)

    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)

    y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)

    # Vectorized log integrand
    log_int <- function(y) {
      # All operations vectorized
      sigma <- sqrt((2 * beta_n) / y)
      sd_mu <- sigma / sqrt(k_n)

      # Get bounds for all sigma values at once
      mb <- constr$mu_b_fn_vec(sigma)
      mb_L <- mb$lower
      mb_U <- mb$upper

      # Identify valid intervals
      valid <- mb_L < mb_U

      # Initialize result with -Inf

      result <- rep(-Inf, length(y))

      if (!any(valid)) return(result)

      # Compute z-scores (vectorized)
      z_U <- ifelse(is.infinite(mb_U), Inf, (mb_U - mu_n) / sd_mu)
      z_L <- ifelse(is.infinite(mb_L), -Inf, (mb_L - mu_n) / sd_mu)

      # log_diff_exp(stats::pnorm(z_U, log.p=TRUE), stats::pnorm(z_L, log.p=TRUE)) + stats::dchisq(y, df_p, log=TRUE)
      log_prob <- log_diff_exp(stats::pnorm(z_U, log.p = TRUE), stats::pnorm(z_L, log.p = TRUE))
      result[valid] <- log_prob[valid] + stats::dchisq(y[valid], df_p, log = TRUE)
      result
    }

    h_max <- h_max_global

    safe_integrand <- function(y) {
      vals <- exp(log_int(y) - h_max)
      vals[!is.finite(vals)] <- 0
      return(vals)
    }

    res <- stats::integrate(safe_integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    return(exp(h_max + log(res)))
  }
}

#' @export
make_solver.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                              metric = "Cpk", target = NULL,
                                              cached_state = NULL, ...) {
  # Case 2: Conjugate mu (Normal/Jeffreys), non-conjugate sigma
  # Integrate mu analytically using pnorm, then 1D numerical over sigma

  ss <- .extract_suff_stats(data, cached_state)
  n <- ss$n; x_bar <- ss$x_bar; sse <- ss$SS
  smp <- .semi_mu_posterior(prior, n, x_bar, sse)
  k_n <- smp$k_n; mu_n <- smp$mu_n; sse_n <- smp$sse_n

  n_eff <- if (prior$k0 == 0) n - 1 else n

  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)

    integrand <- function(sigma) {
      log_lik <- -n_eff * log(sigma) - sse_n / (2 * sigma^2)

      # Log prior on sigma (non-conjugate)
      log_prior_sigma <- prior$log_dens_sigma(sigma)

      # P(mu in constraint region | sigma, data)
      sd_mu <- sigma / sqrt(k_n)
      mb <- constr$mu_b_fn_vec(sigma)
      z_U <- (mb$upper - mu_n) / sd_mu
      z_L <- (mb$lower - mu_n) / sd_mu

      # Safety check for impossible regions
      log_prob_mu <- ifelse(z_L >= z_U, -Inf,
                             log_diff_exp(stats::pnorm(z_U, log.p = TRUE),
                                          stats::pnorm(z_L, log.p = TRUE)))

      exp(log_lik + log_prior_sigma + log_prob_mu)
    }

    # Integration bounds for sigma
    sigma_upper <- if (is.infinite(s_max)) 20 * sqrt(sse_n / max(n, 1)) else s_max
    sigma_upper <- max(sigma_upper, 1e-6)

    stats::integrate(integrand, 1e-10, sigma_upper,
                     rel.tol = 1e-5, subdivisions = 200)$value
  }
}

#' @export
make_solver.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 cached_state = NULL, ...) {
  # Case 3: Non-conjugate mu, conjugate sigma (InvGamma/Jeffreys)
  # Integrate sigma analytically using Gamma functions, then 1D over mu

  # Use cached state for sufficient statistics (data may not be in scope)
  if (!is.null(cached_state)) {
    n <- cached_state$n
    x_bar <- cached_state$x_bar
    sse <- cached_state$sse
  } else {
    n <- length(data)
    x_bar <- mean(data)
    sse <- sum((data - x_bar)^2)
  }
  alpha0 <- prior$alpha0
  beta0 <- prior$beta0

  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)

    integrand <- function(mu) {
      sse_mu <- sse + n * (mu - x_bar)^2
      alpha_n <- alpha0 + n / 2
      beta_n <- beta0 + sse_mu / 2
      log_marginal <- lgamma(alpha_n) - lgamma(alpha0) +
                      alpha0 * log(beta0) - alpha_n * log(beta_n)
      log_prior_mu <- prior$log_dens_mu(mu)

      if (is.infinite(s_max)) {
        log_prob_sigma <- 0
      } else {
        prob_sigma <- 1 - stats::pgamma(1 / s_max^2, shape = alpha_n, rate = beta_n)
        log_prob_sigma <- log(pmax(prob_sigma, 1e-300))
      }

      exp(log_marginal + log_prior_mu + log_prob_sigma)
    }

    mu_sd <- sqrt(sse / max(n - 1, 1))
    mu_low <- x_bar - 10 * mu_sd
    mu_high <- x_bar + 10 * mu_sd

    stats::integrate(integrand, mu_low, mu_high,
                     rel.tol = 1e-5, subdivisions = 200)$value
  }
}

# ==============================================================================
# Density Solver Factory: Direct PDF via Contour Integration
# ==============================================================================

#' Create density solver for direct PDF computation via contour integration
#'
#' Instead of computing P(metric > c) and differentiating, this computes the
#' marginal PDF p(c) directly by integrating along the contour lines where
#' metric = c. This is faster and avoids finite-difference noise.
#'
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate)
#' @param metric Capability metric (currently "Cpk", "Cpu", "Cpl", "Cp" supported)
#' @param target Target value for Cpm/Cpc
#' @param ... Additional arguments
#' @return Function pdf(c) returning the marginal posterior density at c
#' @keywords internal
make_density_solver <- function(data, LSL, USL, prior, metric = "Cpk", target = NULL, ...) {
  UseMethod("make_density_solver", prior)
}

#' @export
#' @export
make_density_solver.PriorConjugate <- function(data, LSL, USL, prior,
                                                metric = "Cpk", target = NULL, cached_state = NULL, ...) {
  ss   <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  M <- (LSL + USL) / 2
  tol <- USL - LSL

  log_Z_sigma <- lgamma(alpha_n) - log(2) - alpha_n * log(beta_n)

  # Cp: no integration needed (contour is a single sigma point)
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_c <- tol / (6 * c)
      if (sigma_c <= 0) return(0)
      log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_c) - beta_n / sigma_c^2
      log_jacobian <- log(tol / 6) - 2 * log(c)
      exp(log_p_sigma + log_jacobian - log_Z_sigma)
    })
  }

  # Precompute sigma grid for metrics requiring contour integration.
  # Replaces per-call stats::integrate with a fixed trapezoidal rule.
  sigma_mode_approx <- sqrt(beta_n / max(alpha_n, 1))
  sigma_lo <- max(1e-10, sigma_mode_approx * 0.02)
  sigma_hi <- sigma_mode_approx * 10
  n_sigma <- 1024L
  sigma_grid <- seq(sigma_lo, sigma_hi, length.out = n_sigma)
  d_sigma <- sigma_grid[2] - sigma_grid[1]
  lps_grid <- -(2 * alpha_n + 1) * log(sigma_grid) - beta_n / sigma_grid^2
  log_jac_grid <- log(3) + log(sigma_grid)
  sd_mu_grid <- sigma_grid / sqrt(k_n)

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / (6 * c)
      if (sigma_max <= sigma_lo) return(0)

      # Include sigma_max as the boundary point for smooth trapezoidal integration.
      idx <- which(sigma_grid < sigma_max)
      if (length(idx) == 0L) return(0)
      sig_pts <- c(sigma_grid[idx], min(sigma_max, sigma_hi))
      np <- length(sig_pts)
      ds <- diff(sig_pts)

      # Trapezoidal weights
      w_trap <- numeric(np)
      w_trap[1L] <- ds[1L] / 2
      w_trap[np] <- ds[np - 1L] / 2
      if (np > 2L) w_trap[2L:(np - 1L)] <- (ds[1L:(np - 2L)] + ds[2L:(np - 1L)]) / 2

      sd_mu <- sig_pts / sqrt(k_n)
      mu_L <- LSL + 3 * c * sig_pts
      mu_U <- USL - 3 * c * sig_pts
      lps <- c(lps_grid[idx], -(2 * alpha_n + 1) * log(sig_pts[np]) - beta_n / sig_pts[np]^2)
      lj  <- c(log_jac_grid[idx], log(3) + log(sig_pts[np]))
      log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
      log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)
      contrib_L <- exp(lps + log_p_mu_L + lj - log_Z_sigma)
      contrib_U <- exp(lps + log_p_mu_U + lj - log_Z_sigma)
      contrib_L[!is.finite(contrib_L)] <- 0
      contrib_U[!is.finite(contrib_U)] <- 0
      sum((contrib_L + contrib_U) * w_trap)
    })
  }

  if (metric == "Cpu") {
    return(function(c) {
      if (c <= 0) return(0)
      mu_U <- USL - 3 * c * sigma_grid
      log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu_grid, log = TRUE)
      vals <- exp(lps_grid + log_p_mu_U + log_jac_grid - log_Z_sigma)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  if (metric == "Cpl") {
    return(function(c) {
      if (c <= 0) return(0)
      mu_L <- LSL + 3 * c * sigma_grid
      log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu_grid, log = TRUE)
      vals <- exp(lps_grid + log_p_mu_L + log_jac_grid - log_Z_sigma)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  if (metric %in% c("Cpm", "Cpc")) {
    T_val <- if (metric == "Cpc") M else target
    if (is.null(T_val)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)
      R <- tol / (6 * c)
      if (R <= 0) return(0)
      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      valid <- sigma_v > 1e-10
      if (!any(valid)) return(0)
      sv <- sigma_v[valid]
      mu_v <- T_val + R * cos(theta_grid[valid])
      sd_mu <- sv / sqrt(k_n)

      log_p_sigma <- -(2 * alpha_n + 1) * log(sv) - beta_n / sv^2
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)
      vals <- exp(log_p_sigma + log_p_mu + log_jacobian - log_Z_sigma)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' @export
make_density_solver.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                                      metric = "Cpk", target = NULL,
                                                      cached_state = NULL, ...) {
  # Case 2: Contour integration with semi-analytical mu integration
  # Uses the same contour approach as PriorConjugate but with
  # log p(sigma | data) = log_lik(sigma) + log_prior_sigma(sigma)

  ss <- .extract_suff_stats(data, cached_state)
  n <- ss$n; x_bar <- ss$x_bar; sse <- ss$SS
  smp <- .semi_mu_posterior(prior, n, x_bar, sse)
  k_n <- smp$k_n; mu_n <- smp$mu_n; sse_n <- smp$sse_n

  M <- (LSL + USL) / 2
  tol <- USL - LSL

  # The sigma exponent in the marginal likelihood depends on whether the mu
  # prior contributes a sigma^{-1} factor. For k0 > 0 (NIG mu prior), the
  # N(mu0, sigma^2/k0) prior has sigma^{-1}, balanced by the conditional.
  # For k0 = 0 (flat mu prior), there is no sigma^{-1} from the prior, so
  # we use n-1 instead of n to avoid double-counting.
  n_eff <- if (prior$k0 == 0) n - 1 else n

  # Compute normalization constant Z for the marginal posterior of sigma
  log_Z <- .compute_semi_mu_log_Z(n_eff, sse_n, prior$log_dens_sigma)

  # Vectorized log marginal posterior of sigma (unnormalized)
  log_p_sigma <- function(sigma_v) {
    -n_eff * log(sigma_v) - sse_n / (2 * sigma_v^2) + prior$log_dens_sigma(sigma_v)
  }

  # Precompute fine sigma grid for numerical integration.
  # stats::integrate can miss the narrow peak in the contour integrand
  # (width ~ sigma / (3*c*sqrt(n))) when the integration range is wide.
  # A fixed grid with sufficient density guarantees peak detection.
  sigma_mode_approx <- sqrt(sse_n / max(n_eff + 1, 1))
  sigma_lo <- max(1e-10, sigma_mode_approx * 0.05)
  sigma_hi <- sigma_mode_approx * 8
  n_sigma <- 1024L
  sigma_grid <- seq(sigma_lo, sigma_hi, length.out = n_sigma)
  d_sigma <- sigma_grid[2] - sigma_grid[1]
  lps_grid <- log_p_sigma(sigma_grid)
  log_jac_grid <- log(3) + log(sigma_grid)

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / (6 * c)
      if (sigma_max <= 0) return(0)

      valid <- sigma_grid > 0 & sigma_grid < sigma_max
      if (!any(valid)) return(0)
      sg <- sigma_grid[valid]
      sd_mu <- sg / sqrt(k_n)

      mu_L <- LSL + 3 * c * sg
      mu_U <- USL - 3 * c * sg

      lps <- lps_grid[valid]
      log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
      log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)
      lj <- log_jac_grid[valid]

      contrib_L <- exp(lps + log_p_mu_L + lj - log_Z)
      contrib_U <- exp(lps + log_p_mu_U + lj - log_Z)
      contrib_L[!is.finite(contrib_L)] <- 0
      contrib_U[!is.finite(contrib_U)] <- 0

      sum(contrib_L + contrib_U) * d_sigma
    })
  }

  if (metric == "Cpu") {
    return(function(c) {
      if (c <= 0) return(0)

      mu_U <- USL - 3 * c * sigma_grid
      sd_mu <- sigma_grid / sqrt(k_n)

      log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)
      vals <- exp(lps_grid + log_p_mu_U + log_jac_grid - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  if (metric == "Cpl") {
    return(function(c) {
      if (c <= 0) return(0)

      mu_L <- LSL + 3 * c * sigma_grid
      sd_mu <- sigma_grid / sqrt(k_n)

      log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
      vals <- exp(lps_grid + log_p_mu_L + log_jac_grid - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  # Cp: sigma is fixed on the contour, no integration needed
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_c <- tol / (6 * c)
      if (sigma_c <= 0) return(0)

      lps <- log_p_sigma(sigma_c)
      log_jacobian <- log(tol / 6) - 2 * log(c)
      exp(lps + log_jacobian - log_Z)
    })
  }

  if (metric %in% c("Cpm", "Cpc")) {
    T_val <- if (metric == "Cpc") M else target
    if (is.null(T_val)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)

      R <- tol / (6 * c)
      if (R <= 0) return(0)

      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      mu_v <- T_val + R * cos(theta_grid)
      sd_mu <- sigma_v / sqrt(k_n)

      lps <- log_p_sigma(sigma_v)
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)

      vals <- exp(lps + log_p_mu + log_jacobian - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' Helper to compute log normalization constant for semi-conjugate mu prior
#' @keywords internal
.compute_semi_mu_log_Z <- function(n_eff, sse_n, log_dens_sigma) {
  integrand <- function(sigma) {
    log_lik <- -n_eff * log(sigma) - sse_n / (2 * sigma^2)
    log_prior <- log_dens_sigma(sigma)
    exp(log_lik + log_prior)
  }
  Z <- stats::integrate(integrand, 1e-10, Inf, rel.tol = 1e-5)$value
  log(max(Z, 1e-300))
}

#' @export
make_density_solver.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                         metric = "Cpk", target = NULL,
                                                         cached_state = NULL, ...) {
  # Case 3: Contour integration with semi-analytical sigma integration
  # For each mu on the contour, p(sigma|mu,data) is Inverse-Gamma

  # Use cached state for sufficient statistics (data may not be in scope)
  if (!is.null(cached_state)) {
    n <- cached_state$n
    x_bar <- cached_state$x_bar
    sse <- cached_state$sse
  } else {
    n <- length(data)
    x_bar <- mean(data)
    sse <- sum((data - x_bar)^2)
  }
  alpha0 <- prior$alpha0
  beta0 <- prior$beta0
  alpha_n <- alpha0 + n / 2

  # In the NIG model, alpha0=-0.5 with beta0=0 gives a FLAT prior on sigma,

  # but the true Jeffreys prior is sigma^{-1}. The missing sigma^{-1} factor
  # is equivalent to shifting alpha_n by +0.5 in the sigma exponent and Z.
  # We keep alpha0=-0.5 (since lgamma(-0.5) is finite) and only adjust alpha_n.
  jeffreys_adj <- if (alpha0 == -0.5 && beta0 == 0) 0.5 else 0
  alpha_n_eff <- alpha_n + jeffreys_adj

  M <- (LSL + USL) / 2
  tol <- USL - LSL

  # Compute normalization constant Z
  log_Z <- .compute_semi_sigma_log_Z(n, x_bar, sse, alpha0, beta0,
                                      prior$log_dens_mu, jeffreys_adj)

  # Constant for the unnormalized sigma kernel
  alpha_0_times_logbeta0 <- if (beta0 == 0) 0 else alpha0 * log(beta0)
  log_C <- log(2) + alpha_0_times_logbeta0 - lgamma(alpha0)

  # Precompute fine sigma grid for numerical integration.
  # The contour integrand has a narrow peak in sigma from the likelihood
  # factor exp(-n*(mu-xbar)^2/(2*sigma^2)) which stats::integrate can miss.
  sigma_mode_approx <- sqrt(sse / max(n - 1, 1))
  sigma_lo <- max(1e-10, sigma_mode_approx * 0.05)
  sigma_hi <- sigma_mode_approx * 8
  n_sigma <- 1024L
  sigma_grid <- seq(sigma_lo, sigma_hi, length.out = n_sigma)
  d_sigma <- sigma_grid[2] - sigma_grid[1]
  log_sigma_grid <- log(sigma_grid)
  log_jac_grid <- log(3) + log_sigma_grid

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / (6 * c)
      if (sigma_max <= 0) return(0)

      valid <- sigma_grid > 0 & sigma_grid < sigma_max
      if (!any(valid)) return(0)
      sg <- sigma_grid[valid]

      mu_L <- LSL + 3 * c * sg
      mu_U <- USL - 3 * c * sg

      sse_L <- sse + n * (mu_L - x_bar)^2
      sse_U <- sse + n * (mu_U - x_bar)^2
      beta_n_L <- beta0 + sse_L / 2
      beta_n_U <- beta0 + sse_U / 2

      lsg <- log_sigma_grid[valid]
      log_sk_L <- log_C - (2 * alpha_n_eff + 1) * lsg - beta_n_L / sg^2
      log_sk_U <- log_C - (2 * alpha_n_eff + 1) * lsg - beta_n_U / sg^2

      log_pm_L <- prior$log_dens_mu(mu_L)
      log_pm_U <- prior$log_dens_mu(mu_U)
      lj <- log_jac_grid[valid]

      contrib_L <- exp(log_pm_L + log_sk_L + lj - log_Z)
      contrib_U <- exp(log_pm_U + log_sk_U + lj - log_Z)
      contrib_L[!is.finite(contrib_L)] <- 0
      contrib_U[!is.finite(contrib_U)] <- 0

      sum(contrib_L + contrib_U) * d_sigma
    })
  }

  if (metric == "Cpu") {
    return(function(c) {
      if (c <= 0) return(0)

      mu_U <- USL - 3 * c * sigma_grid
      sse_U <- sse + n * (mu_U - x_bar)^2
      beta_n_U <- beta0 + sse_U / 2

      log_sk <- log_C - (2 * alpha_n_eff + 1) * log_sigma_grid - beta_n_U / sigma_grid^2
      log_pm <- prior$log_dens_mu(mu_U)

      vals <- exp(log_pm + log_sk + log_jac_grid - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  if (metric == "Cpl") {
    return(function(c) {
      if (c <= 0) return(0)

      mu_L <- LSL + 3 * c * sigma_grid
      sse_L <- sse + n * (mu_L - x_bar)^2
      beta_n_L <- beta0 + sse_L / 2

      log_sk <- log_C - (2 * alpha_n_eff + 1) * log_sigma_grid - beta_n_L / sigma_grid^2
      log_pm <- prior$log_dens_mu(mu_L)

      vals <- exp(log_pm + log_sk + log_jac_grid - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_sigma
    })
  }

  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_c <- tol / (6 * c)
      if (sigma_c <= 0) return(0)

      log_jacobian <- log(tol / 6) - 2 * log(c)

      integrand <- function(mu) {
        sse_mu <- sse + n * (mu - x_bar)^2
        beta_n <- beta0 + sse_mu / 2
        log_sigma_kernel <- log_C - (2 * alpha_n_eff + 1) * log(sigma_c) - beta_n / sigma_c^2
        log_prior_mu <- prior$log_dens_mu(mu)

        exp(log_prior_mu + log_sigma_kernel + log_jacobian - log_Z)
      }

      mu_lower <- x_bar - 20 * sqrt(sse / max(n, 1))
      mu_upper <- x_bar + 20 * sqrt(sse / max(n, 1))
      stats::integrate(integrand, mu_lower, mu_upper, rel.tol = 1e-5)$value
    })
  }

  # Cpm/Cpc: Semicircular contours
  if (metric %in% c("Cpm", "Cpc")) {
    T_val <- if (metric == "Cpc") M else target
    if (is.null(T_val)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)

      R <- tol / (6 * c)
      if (R <= 0) return(0)

      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      mu_v <- T_val + R * cos(theta_grid)

      sse_mu <- sse + n * (mu_v - x_bar)^2
      beta_n <- beta0 + sse_mu / 2
      log_sk <- log_C - (2 * alpha_n_eff + 1) * log(sigma_v) - beta_n / sigma_v^2
      log_pm <- prior$log_dens_mu(mu_v)

      vals <- exp(log_pm + log_sk + log_jacobian - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' Helper to compute log normalization constant for semi-conjugate sigma prior
#' @keywords internal
.compute_semi_sigma_log_Z <- function(n, x_bar, sse, alpha0, beta0, log_dens_mu,
                                      jeffreys_adj = 0) {
  alpha_n <- alpha0 + n / 2
  alpha_n_eff <- alpha_n + jeffreys_adj

  # for the Jeffreys prior case, beta0 = 0, so we need to handle that carefully
  alpha_0_times_logbeta0 <- if (beta0 == 0) 0 else alpha0 * log(beta0)


  integrand <- function(mu) {
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    log_marginal <- lgamma(alpha_n_eff) - lgamma(alpha0) + alpha_0_times_logbeta0 - alpha_n_eff * log(beta_n)
    log_prior <- log_dens_mu(mu)
    exp(log_marginal + log_prior)
  }

  mu_lower <- x_bar - 20 * sqrt(sse / max(n, 1))
  mu_upper <- x_bar + 20 * sqrt(sse / max(n, 1))
  Z <- stats::integrate(integrand, mu_lower, mu_upper, rel.tol = 1e-5)$value
  log(max(Z, 1e-300))
}

#' @export
make_density_solver.PriorGeneric <- function(data, LSL, USL, prior,
                                              metric = "Cpk", target = NULL,
                                              cached_state = NULL, ...) {
  # For generic priors, we still use 1D integration along contours,
  # but evaluate the joint posterior numerically using vectorized operations

  if (is.null(cached_state)) {
    cached_state <- precompute_generic_state(data, prior)
  }

  n <- cached_state$n
  x_bar <- cached_state$x_bar
  h_max <- cached_state$h_max
  Z <- cached_state$Z
  uni_s <- cached_state$uni_s
  map_mu <- cached_state$map_mu
  mu_scale <- cached_state$mu_scale
  log_post_vec <- cached_state$log_post_vec

  M <- (LSL + USL) / 2
  tol <- USL - LSL

  # Return pdf(c) function for Cpk
  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_max <- tol / (6 * c)
      if (sigma_max <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0 & sigma < sigma_max
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]

        # Contour points
        mu_L <- LSL + 3 * c * sigma_v
        mu_U <- USL - 3 * c * sigma_v

        # Jacobian |dmu/dc| = 3*sigma
        log_jacobian <- log(3) + log(sigma_v)

        # Evaluate joint posterior at contour points (vectorized)
        log_p_L <- log_post_vec(mu_L, sigma_v)
        log_p_U <- log_post_vec(mu_U, sigma_v)

        # Contributions from both contours (normalized by Z)
        # Z = ∫ exp(log_p - h_max) dmu dsigma, so we divide by Z
        contrib_L <- exp(log_p_L + log_jacobian - h_max) / Z
        contrib_U <- exp(log_p_U + log_jacobian - h_max) / Z

        contrib_L[!is.finite(contrib_L)] <- 0
        contrib_U[!is.finite(contrib_U)] <- 0

        result[valid] <- contrib_L + contrib_U
        result
      }

      stats::integrate(integrand, 0, sigma_max, rel.tol = 1e-5, subdivisions = 200)$value
    })
  }

  # Cpu: single contour
  # For Cpu = c, contour is mu = USL - 3*c*sigma
  # Need sigma_upper that ensures mu is within prior support (or has meaningful density)
  if (metric == "Cpu") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        mu_U <- USL - 3 * c * sigma_v
        log_jacobian <- log(3) + log(sigma_v)

        log_p <- log_post_vec(mu_U, sigma_v)
        contrib <- exp(log_p + log_jacobian - h_max) / Z
        contrib[!is.finite(contrib)] <- 0

        result[valid] <- contrib
        result
      }

      mu_center <- if (n > 0) x_bar else map_mu
      mu_sd <- if (n > 0) sqrt(uni_s^2 / n) else if (!is.null(mu_scale)) mu_scale else uni_s
      mu_lower_bound <- mu_center - 6 * mu_sd
      sigma_upper_mu <- (USL - mu_lower_bound) / (3 * c)
      sigma_upper <- min(uni_s, max(sigma_upper_mu, 1e-6))

      tryCatch({
        stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-5, subdivisions = 200)$value
      }, error = function(e) {
        # Fallback: return 0 if integration fails (e.g., at extreme c values)
        0
      })
    })
  }

  # Cpl: single contour
  # For Cpl = c, contour is mu = LSL + 3*c*sigma
  if (metric == "Cpl") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        mu_L <- LSL + 3 * c * sigma_v
        log_jacobian <- log(3) + log(sigma_v)

        log_p <- log_post_vec(mu_L, sigma_v)
        contrib <- exp(log_p + log_jacobian - h_max) / Z
        contrib[!is.finite(contrib)] <- 0

        result[valid] <- contrib
        result
      }

      mu_center <- if (n > 0) x_bar else map_mu
      mu_sd <- if (n > 0) sqrt(uni_s^2 / n) else if (!is.null(mu_scale)) mu_scale else uni_s
      mu_upper_bound <- mu_center + 6 * mu_sd
      sigma_upper_mu <- (mu_upper_bound - LSL) / (3 * c)
      sigma_upper <- min(uni_s, max(sigma_upper_mu, 1e-6))

      tryCatch({
        stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-5, subdivisions = 200)$value
      }, error = function(e) {
        0
      })
    })
  }

  # Cp: single point
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_c <- tol / (6 * c)
      if (sigma_c <= 0) return(0)

      # At sigma_c, mu can be anything - integrate over mu
      # p(Cp = c) = ∫ p(mu, sigma_c) |dsigma/dc| dmu
      # |dsigma/dc| = tol / (6 * c^2)
      log_jacobian <- log(tol / 6) - 2 * log(c)

      # Integrate over mu
      mu_integrand <- function(mu) {
        sigma_vec <- rep(sigma_c, length(mu))
        log_p <- log_post_vec(mu, sigma_vec)
        exp(log_p - h_max) / Z
      }

      if (n > 0) {
        sd_mu_cond <- sigma_c / sqrt(n)
        mu_width <- 15 * sd_mu_cond
        mu_lower <- x_bar - mu_width
        mu_upper <- x_bar + mu_width
      } else {
        mc <- map_mu
        ms <- mu_scale %||% uni_s
        mu_lower <- mc - 6 * ms
        mu_upper <- mc + 6 * ms
      }

      mu_integral <- stats::integrate(mu_integrand, mu_lower, mu_upper, rel.tol = 1e-5)$value
      mu_integral * exp(log_jacobian)
    })
  }

  # Cpm/Cpc: Semicircular contours
  if (metric %in% c("Cpm", "Cpc")) {
    T_val <- if (metric == "Cpc") M else target
    if (is.null(T_val)) stop("Target required for Cpm")

    return(function(c) {
      if (c <= 0) return(0)

      # Contour radius R = tol / (6*c)
      R <- tol / (6 * c)
      if (R <= 0) return(0)

      # Jacobian determinant for polar coordinates wrt c
      # |J| = R^2 / c
      log_jacobian <- 2 * log(R) - log(c)

      integrand <- function(theta) {
        # theta in [0, pi] covers sigma > 0
        sigma_v <- R * sin(theta)

        # Avoid singular boundary at sigma=0
        valid <- sigma_v > 1e-10
        result <- rep(0, length(theta))
        if (!any(valid)) return(result)

        sigma_v <- sigma_v[valid]
        mu_v <- T_val + R * cos(theta[valid])

        # Log joint posterior
        log_p <- log_post_vec(mu_v, sigma_v)

        # Combined density (normalized by Z)
        log_contrib <- log_p + log_jacobian - h_max

        result[valid] <- exp(log_contrib) / Z
        result
      }

      stats::integrate(integrand, 0, pi, rel.tol = 1e-5, subdivisions = 200)$value
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}


#' @export
make_solver.PriorGeneric <- function(data, LSL, USL, prior, metric = "Cpk",
                                      target = NULL, cached_state = NULL, ...) {

  # Use cached state if available (massive speedup for multiple metrics)
  if (!is.null(cached_state)) {
    log_post <- cached_state$log_post
    h_max <- cached_state$h_max
    uni_s <- cached_state$uni_s
    map_mu <- cached_state$map_mu
    Z <- cached_state$Z
    int_2d <- cached_state$int_2d
  } else {
    # Compute from scratch
    n <- length(data)
    x_bar <- mean(data)
    sse <- sum((data - x_bar)^2)

    # Scalar log posterior (for optim)
    log_post <- function(mu, sigma) {
      if (sigma <= 0) return(-Inf)
      -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2) +
        prior$log_dens(mu, sigma)
    }

    # Vectorized log posterior (for cubature)
    log_post_vec <- function(mu, sigma) {
      log_lik <- rep(-Inf, length(mu))
      valid <- sigma > 0
      if (any(valid)) {
        log_lik[valid] <- -n * log(sigma[valid]) -
                          (sse + n * (mu[valid] - x_bar)^2) / (2 * sigma[valid]^2)
      }
      log_prior <- prior$log_dens(mu, sigma)
      log_lik + log_prior
    }

    # Find MAP for integration bounds
    init_sd <- sqrt(sse / (n - 1))
    opt <- stats::optim(c(x_bar, init_sd), function(p) -log_post(p[1], p[2]))
    map_mu <- opt$par[1]
    map_sig <- opt$par[2]
    h_max <- -opt$value

    # Tighter bounds: 5 sigma from MAP
    uni_s <- map_sig * 5

    # cubature 2D integration with fully vectorized interface
    int_2d <- function(s_lim, m_fn) {
      s_top <- if (is.infinite(s_lim)) uni_s else min(s_lim, uni_s)

      # Fully vectorized integrand: x is 2 x n matrix
      integrand <- function(x) {
        mu <- x[1, ]
        sigma <- x[2, ]

        # Compute log posterior (vectorized)
        log_p <- log_post_vec(mu, sigma)

        # Apply constraint mask from m_fn (always returns list(lower, upper))
        mb <- m_fn(sigma)
        valid <- (mu >= mb$lower) & (mu <= mb$upper)

        result <- rep(0, length(mu))
        result[valid] <- exp(log_p[valid] - h_max)
        result[!is.finite(result)] <- 0
        matrix(result, nrow = 1)
      }


      # Tight mu bounds around MAP
      mu_lower <- map_mu - 5 * uni_s
      mu_upper <- map_mu + 5 * uni_s

      result <- cubature::pcubature(
        integrand,
        lowerLimit = c(mu_lower, 1e-10),
        upperLimit = c(mu_upper, s_top),
        tol = 1e-3,
        vectorInterface = TRUE
      )
      result$integral
    }

    Z <- int_2d(Inf, function(s) list(lower = rep(-Inf, length(s)), upper = rep(Inf, length(s))))
  }

  # Return closure
  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    num <- int_2d(constr$s_max_fn(), constr$mu_b_fn_vec)
    return(num / Z)
  }
}


# ==============================================================================
# Pre-computation for Generic Priors
# ==============================================================================

#' Precompute expensive posterior state for Generic priors
#' @param data Numeric vector of observations
#' @param prior PriorGeneric object
#' @return cached_state object to pass to make_solver
#' @keywords internal
precompute_generic_state <- function(data, prior) {
  n <- length(data)
  if (n > 0) {
    x_bar <- mean(data)
    sse <- sum((data - x_bar)^2)
    init_sd <- sqrt(sse / (n - 1))
    init_mu <- x_bar
  } else {
    x_bar <- 0
    sse <- 0
    bt <- prior$bayestools_priors
    if (!is.null(bt)) {
      mu_init  <- .extract_prior_init(bt$mu)
      sig_init <- .extract_prior_init(bt$sigma)
      init_mu <- mu_init$value
      init_sd <- sig_init$value
    } else {
      init_mu <- 0
      init_sd <- 1
    }
  }

  # Scalar log posterior (used by optim)
  log_post <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    log_lik <- if (n > 0) {
      -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2)
    } else {
      0
    }
    log_lik + prior$log_dens(mu, sigma)
  }

  # Vectorized log posterior (used by cubature)
  log_post_vec <- function(mu, sigma) {
    log_lik <- rep(0, length(mu))
    if (n > 0) {
      log_lik[] <- -Inf
      valid <- sigma > 0
      if (any(valid)) {
        log_lik[valid] <- -n * log(sigma[valid]) -
                          (sse + n * (mu[valid] - x_bar)^2) / (2 * sigma[valid]^2)
      }
    }
    log_prior <- prior$log_dens(mu, sigma)
    log_lik + log_prior
  }

  # Find MAP using constrained optimization (sigma > 0)
  # Use L-BFGS-B with reasonable bounds
  opt <- tryCatch({
    stats::optim(c(init_mu, init_sd), function(p) -log_post(p[1], p[2]),
                 method = "L-BFGS-B",
                 lower = c(-1e6, 1e-6),
                 upper = c(1e6, 1e6))
  }, error = function(e) {
    # Fallback to Nelder-Mead if L-BFGS-B fails
    stats::optim(c(init_mu, max(init_sd, 0.1)), function(p) {
      if (p[2] <= 0) return(1e10)
      -log_post(p[1], p[2])
    })
  })
  map_mu <- opt$par[1]
  map_sig <- max(opt$par[2], 1e-3)  # Ensure positive
  h_max <- -opt$value

  # Compute integration bounds based on posterior concentration
  # For Normal likelihood with data:
  #   - Posterior of sigma is concentrated around map_sig
  #   - Effective SD of sigma is roughly map_sig / sqrt(2*n)
  #   - Posterior of mu|sigma is N(x_bar, sigma^2/n)
  # Use bounds that capture ~6 SDs of the posterior to get >99.99% of mass
  if (n > 0) {
    # Sigma bounds: use multiplier based on posterior concentration
    # sigma_sd ~ map_sig / sqrt(2n), so 6 sigma_sd ~ 6 * map_sig / sqrt(2n)
    # For n=30, this is about 0.77 * map_sig, so bounds [0.5 * map_sig, 3 * map_sig] should work
    sigma_lower <- max(0.2 * map_sig, 1e-6)
    sigma_upper <- 4 * map_sig

    # Mu bounds: at sigma_upper, SD of mu|sigma is sigma_upper/sqrt(n)
    # Use 6 SDs for safety
    mu_sd_at_sigma_upper <- sigma_upper / sqrt(n)
    mu_lower <- map_mu - 6 * mu_sd_at_sigma_upper
    mu_upper <- map_mu + 6 * mu_sd_at_sigma_upper

    # For other functions (like make_solver), uni_s is used as a reference scale
    uni_s <- sigma_upper
  } else {
    # No data (prior-only): derive bounds from prior quantiles
    bt <- prior$bayestools_priors
    if (!is.null(bt) && inherits(bt$sigma, "prior")) {
      sig_q_lo <- BayesTools::quant(bt$sigma, 0.0001)
      sig_q_hi <- BayesTools::quant(bt$sigma, 0.9999)
      sigma_lower <- max(1e-6, sig_q_lo)
      sigma_upper <- sig_q_hi
    } else if (map_sig > 0.1 && map_sig < 100) {
      sigma_lower <- max(0.01, map_sig / 20)
      sigma_upper <- min(100, map_sig * 20)
    } else {
      sigma_lower <- 0.01
      sigma_upper <- 100
    }
    uni_s <- sigma_upper
    if (!is.null(bt)) {
      mu_init_info <- .extract_prior_init(bt$mu)
      mu_scale <- mu_init_info$scale
    } else {
      mu_scale <- uni_s
    }
    mu_lower <- map_mu - 6 * mu_scale
    mu_upper <- map_mu + 6 * mu_scale
  }

  # cubature 2D integration with fully vectorized interface
  int_2d <- function(s_lim, m_fn) {
    s_top <- if (is.infinite(s_lim)) sigma_upper else min(s_lim, sigma_upper)

    # Fully vectorized integrand: x is 2 x n matrix
    integrand <- function(x) {
      mu <- x[1, ]
      sigma <- x[2, ]

      # Compute log posterior (vectorized)
      log_p <- log_post_vec(mu, sigma)

      # Apply constraint mask from m_fn (always returns list(lower, upper))
      mb <- m_fn(sigma)
      valid <- (mu >= mb$lower) & (mu <= mb$upper)

      result <- rep(0, length(mu))
      result[valid] <- exp(log_p[valid] - h_max)
      result[!is.finite(result)] <- 0
      matrix(result, nrow = 1)
    }

    result <- cubature::pcubature(
      integrand,
      lowerLimit = c(mu_lower, max(sigma_lower, 1e-10)),
      upperLimit = c(mu_upper, s_top),
      tol = 1e-4,
      vectorInterface = TRUE
    )
    result$integral
  }

  Z <- int_2d(Inf, function(s) list(lower = rep(-Inf, length(s)), upper = rep(Inf, length(s))))

  # mu_scale: prior SD for mu (used by density solver for bounds when n=0)
  if (!exists("mu_scale")) mu_scale <- NULL

  list(log_post = log_post, log_post_vec = log_post_vec, h_max = h_max,
       uni_s = uni_s, map_mu = map_mu, mu_scale = mu_scale,
       int_2d = int_2d, Z = Z, n = n,
       x_bar = x_bar, sse = sse)
}


# ==============================================================================
# Main Analysis Functions
# ==============================================================================

#' Compute Posterior Probability that a Capability Index is in a Region
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param bounds Vector c(lower, upper) defining the interval
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param cached_state Pre-computed state from precompute_generic_state
#' @return Probability that metric is in (\code{bounds[1]}, \code{bounds[2]})
#' @keywords internal
compute_cpk_prob_integration <- function(data, LSL, USL, bounds, prior,
                                          metric = "Cpk", target = NULL,
                                          cached_state = NULL) {

  # Check if we can use density solver
  density_metrics <- c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")
  can_use_density <- (inherits(prior, "PriorConjugate") ||
                      inherits(prior, "PriorGeneric") ||
                      inherits(prior, "PriorSemiConjugateMu") ||
                      inherits(prior, "PriorSemiConjugateSigma")) &&
                     metric %in% density_metrics

  if (can_use_density) {
    # Use density solver PDF and integrate over bounds
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                   cached_state = cached_state)
    pdf_vec <- function(x) vapply(x, pdf_fn, numeric(1L))

    # Clip infinite bounds: capability indices are >= 0 and rarely exceed 10
    lower <- max(min(bounds), 0)
    upper <- if (is.finite(max(bounds))) max(bounds) else 10

    prob <- tryCatch(
      stats::integrate(pdf_vec, lower, upper, rel.tol = 1e-4)$value,
      error = function(e) NA_real_
    )
    if (!is.na(prob)) return(prob)

    # Fallback: try survival function approach
    S <- tryCatch({
      if (!is.null(cached_state))
        make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
      else
        make_solver(data, LSL, USL, prior, metric, target)
    }, error = function(e) NULL)
    if (!is.null(S)) {
      p_lower <- tryCatch(S(lower), error = function(e) NA_real_)
      p_upper <- tryCatch(S(upper), error = function(e) NA_real_)
      if (!is.na(p_lower) && !is.na(p_upper)) return(p_lower - p_upper)
    }
    return(NA_real_)

  } else {
    # Fallback: Survival function S(c) = P(Index > c)
    if (!is.null(cached_state)) {
      S <- make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
    } else {
      S <- make_solver(data, LSL, USL, prior, metric, target)
    }

    # P(lower < Index < upper) = P(Index > lower) - P(Index > upper)
    p_lower <- S(min(bounds))
    p_upper <- S(max(bounds))
    return(p_lower - p_upper)
  }
}

#' Compute negative-value mean correction for prior-only mode.
#' For Cpu, Cpl, Cpk, the metric can be negative when mu is outside \eqn{[LSL, USL]}.
#' This returns \eqn{E[metric * I(metric <= 0)]}, a negative correction to add to
#' \eqn{area * E[metric | metric > 0]} to get the full unconditional mean.
#' @keywords internal
.compute_negative_mean_correction <- function(metric, LSL, USL, bayestools_priors) {
  if (is.null(bayestools_priors)) return(0)
  if (metric %in% c("Cp", "Cpm", "Cpc")) return(0)

  log_dens_mu <- .make_prior_log_dens_fn(bayestools_priors$mu)
  log_dens_sigma <- .make_prior_log_dens_fn(bayestools_priors$sigma)

  # E[1/sigma] from the sigma prior
  E_inv_sigma <- tryCatch({
    integrand <- function(s) exp(log_dens_sigma(s)) / s
    stats::integrate(Vectorize(integrand), 1e-10, Inf,
                     rel.tol = 1e-6, subdivisions = 500)$value
  }, error = function(e) Inf)
  if (!is.finite(E_inv_sigma)) return(0)

  tol <- USL - LSL

  if (metric == "Cpu" || metric == "Cpk") {
    # E[(USL - mu) * I(mu > USL)]: negative since USL - mu < 0 when mu > USL
    E_Cpu_neg <- tryCatch({
      integrand <- function(mu) (USL - mu) * exp(log_dens_mu(mu))
      stats::integrate(Vectorize(integrand), USL, USL + 20 * tol,
                       rel.tol = 1e-6)$value
    }, error = function(e) 0)
    corr_Cpu <- (1 / 3) * E_Cpu_neg * E_inv_sigma
  }

  if (metric == "Cpl" || metric == "Cpk") {
    # E[(mu - LSL) * I(mu < LSL)]: negative since mu - LSL < 0 when mu < LSL
    E_Cpl_neg <- tryCatch({
      integrand <- function(mu) (mu - LSL) * exp(log_dens_mu(mu))
      stats::integrate(Vectorize(integrand), LSL - 20 * tol, LSL,
                       rel.tol = 1e-6)$value
    }, error = function(e) 0)
    corr_Cpl <- (1 / 3) * E_Cpl_neg * E_inv_sigma
  }

  switch(metric,
    "Cpu" = corr_Cpu,
    "Cpl" = corr_Cpl,
    "Cpk" = corr_Cpu + corr_Cpl,
    0
  )
}

#' Analyze Capability with Automatic Grid Detection
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param n_grid Number of grid points for density evaluation
#' @param cached_state Pre-computed state from precompute_generic_state
#' @param use_density_solver Logical, whether to use direct density solver (default TRUE for conjugate priors)
#' @return List with metric name, grid data.frame, and stats vector
#' @keywords internal
analyze_capability_integration <- function(data, LSL, USL, prior,
                                            metric = "Cpk", target = NULL,
                                            n_grid = 128,
                                            cached_state = NULL,
                                            use_density_solver = TRUE,
                                            mc_samples = NULL,
                                            divergence_info = NULL) {

  # Default: no divergence
  if (is.null(divergence_info))
    divergence_info <- list(mean_divergent = FALSE, sd_divergent = FALSE,
                            alpha = Inf, reason = NULL)

  # Check if we're in prior-only mode (no data) with PriorGeneric
  # In this case, use quantile-based grid bounds to avoid slow moment computation
  is_prior_only <- if (!is.null(cached_state$n)) cached_state$n == 0L else length(data) == 0
  has_bayestools_priors <- inherits(prior, "PriorGeneric") &&
                           !is.null(prior$bayestools_priors) &&
                           inherits(prior$bayestools_priors$sigma, "prior")

  if (is_prior_only && has_bayestools_priors) {
    # Fast path: compute grid bounds directly from sigma quantiles
    bounds <- .compute_metric_grid_bounds_from_quantiles(
      prior$bayestools_priors$sigma, LSL, USL, metric, target
    )
    x_start <- bounds$x_start
    x_end <- bounds$x_end
  } else if (divergence_info$mean_divergent) {
    # Divergent mean: skip moment computation, use solver-based grid bounds
    x_start <- 0
    x_end <- 3
    try({
      S_fn <- if (!is.null(cached_state)) {
        make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
      } else {
        make_solver(data, LSL, USL, prior, metric, target)
      }
      for (try_limit in c(10, 20, 50, 100, 500, 1000)) {
        if (S_fn(try_limit) < 0.001) {
          x_end <- try_limit
          break
        }
      }
    }, silent = TRUE)
  } else {
    # Standard path: compute posterior moments for grid bounds
    moments <- compute_metric_moments(data, LSL, USL, prior, metric = metric, target = target,
                                       use_analytic = TRUE, cached_state = cached_state)

    # Grid bounds via normal approximation (clip to natural bounds)
    # Handle NaN moments by falling back to reasonable defaults
    if (is.na(moments$mean) || is.na(moments$sd) || moments$sd <= 0) {
      # Fallback: use a wide grid centered around 1.0
      x_start <- 0
      x_end <- 3
    } else {
      x_start <- max(0, moments$mean - 3 * moments$sd)
      x_end <- moments$mean + 3 * moments$sd
    }

    # If bounds seem excessive (divergent mean), refine using quantiles via solver
    if (x_end > 50) {
       # Use S(c) = P(Index > c) to find effective upper bound (q99.9)
       try({
         S_fn <- if (!is.null(cached_state)) {
           make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
         } else {
           make_solver(data, LSL, USL, prior, metric, target)
         }
         # Search for upper bound where prob < 0.001
         for (try_limit in c(10, 20, 50, 100, 1000)) {
           if (S_fn(try_limit) < 0.001) {
             x_end <- try_limit
             break
           }
         }
       }, silent = TRUE)
    }
  }

  # Ensure valid grid (minimum width)
  if (x_end <= x_start || !is.finite(x_end)) {
    x_end <- max(x_start + 3, 3)
  }

  # Prior-only mode with broad distributions needs more grid resolution
  if (is_prior_only) n_grid <- max(n_grid, 1024L)

  # For prior-only mode with BayesTools priors, use direct Monte Carlo sampling.
  # This avoids grid resolution issues with heavy-tailed priors.
  if (is_prior_only && has_bayestools_priors && !is.null(mc_samples)) {
    mu_samples  <- mc_samples$mu
    sig_samples <- mc_samples$sig
    n_mc <- length(mu_samples)

    M <- (LSL + USL) / 2
    tol <- USL - LSL
    if (is.null(target)) target <- M

    metric_samples <- switch(metric,
      "Cp"  = tol / (6 * sig_samples),
      "Cpu" = (USL - mu_samples) / (3 * sig_samples),
      "Cpl" = (mu_samples - LSL) / (3 * sig_samples),
      "Cpk" = pmin((USL - mu_samples), (mu_samples - LSL)) / (3 * sig_samples),
      "Cpm" = pmin(USL - target, target - LSL) /
              (3 * sqrt(sig_samples^2 + (mu_samples - target)^2)),
      "Cpc" = {
        z <- (mu_samples - target) / sig_samples
        E_abs_dev <- sig_samples * sqrt(2 / pi) * exp(-0.5 * z^2) +
                     abs(mu_samples - target) * (1 - 2 * stats::pnorm(-abs(z)))
        tol / (6 * sqrt(pi / 2) * E_abs_dev)
      },
      stop("Unknown metric: ", metric)
    )

    post_mean   <- mean(metric_samples)
    post_median <- stats::median(metric_samples)
    post_sd     <- stats::sd(metric_samples)
    q2.5  <- unname(stats::quantile(metric_samples, 0.025))
    q97.5 <- unname(stats::quantile(metric_samples, 0.975))

    # HDI from samples
    sorted_samples <- sort(metric_samples)
    n_ci <- floor(0.95 * n_mc)
    ci_widths <- sorted_samples[(n_ci + 1):n_mc] - sorted_samples[1:(n_mc - n_ci)]
    best_ci <- which.min(ci_widths)
    hdi_lo <- sorted_samples[best_ci]
    hdi_hi <- sorted_samples[best_ci + n_ci]

    # Override with analytic Inf where divergent
    if (divergence_info$mean_divergent) {
      post_mean <- Inf
      post_sd   <- Inf
    } else if (divergence_info$sd_divergent) {
      post_sd <- Inf
    }

    return(list(
      metric = metric,
      samples = metric_samples,
      area = 1,
      stats = c(Mean = post_mean, Median = post_median, SD = post_sd,
                Q2.5 = q2.5, Q97.5 = q97.5,
                HDI_Lo = hdi_lo, HDI_Hi = hdi_hi),
      divergence_info = divergence_info
    ))
  }

  # Choose method: density solver (direct PDF) vs survival function + finite diff
  # Density solver works for both PriorConjugate and PriorGeneric with supported metrics
  density_metrics <- c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")
  can_use_density <- use_density_solver &&
                     (inherits(prior, "PriorConjugate") ||
                      inherits(prior, "PriorGeneric") ||
                      inherits(prior, "PriorSemiConjugateMu") ||
                      inherits(prior, "PriorSemiConjugateSigma")) &&
                     metric %in% density_metrics

  if (can_use_density) {
    # Direct PDF computation via contour integration
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                   cached_state = cached_state)

    # Adaptive grid: extend if density is non-negligible at edges
    max_extend <- if (is_prior_only) 5L else 3L
    for (attempt in seq_len(max_extend)) {
      grid_x <- seq(x_start, x_end, length.out = n_grid)
      mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2
      pdf_vals <- sapply(mid_x, pdf_fn)
      pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0

      peak <- max(pdf_vals)
      if (peak == 0) break
      edge_threshold <- 0.01 * peak
      need_extend_left  <- pdf_vals[1] > edge_threshold && x_start > 0
      need_extend_right <- pdf_vals[length(pdf_vals)] > edge_threshold

      if (!need_extend_left && !need_extend_right) break
      width <- x_end - x_start
      if (need_extend_left)  x_start <- max(0, x_start - width)
      if (need_extend_right) x_end   <- x_end + width
    }

    # For broad prior-only distributions (x_end/x_start > 20), use a composite
    # grid: log-spaced for the heavy tail to capture the full distribution.
    if (is_prior_only && x_end > 0 && x_start > 0 && x_end / x_start > 20) {
      log_grid <- exp(seq(log(max(x_start, 1e-4)), log(x_end), length.out = n_grid))
      grid_x <- log_grid
      mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2
      dx <- diff(grid_x)
      pdf_vals <- sapply(mid_x, pdf_fn)
      pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0
    }

    # Handle NaN/negative in PDF
    pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0

    # Normalize area to 1.0
    dx <- diff(grid_x)
    area <- sum(pdf_vals * dx)
    if (area > 0) pdf_vals <- pdf_vals / area

  } else {
    # Fallback: Survival function + finite differences
    if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
      S <- make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
    } else {
      S <- make_solver(data, LSL, USL, prior, metric, target)
    }

    # Evaluate grid
    grid_x <- seq(x_start, x_end, length.out = n_grid)
    S_vals <- sapply(grid_x, S)

    # Handle NaN in S_vals
    S_vals[!is.finite(S_vals)] <- 0

    # Compute PDF via finite differences
    pdf_vals <- -diff(S_vals) / diff(grid_x)
    mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2
    dx <- diff(grid_x)

    # Handle NaN/negative in PDF
    pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0

    # Normalize area to 1.0
    area <- sum(pdf_vals * dx)
    if (area > 0) pdf_vals <- pdf_vals / area
  }

  # Compute statistics
  # For prior-only mode, the density only covers c > 0 and area = P(metric > 0)
  # which can be < 1 when the mu prior has mass outside [LSL, USL].
  # MCMC includes negative values, so we must account for this.
  # For data mode, area < 1 is just grid truncation; use standard normalization.
  if (is_prior_only && area < 0.999) {
    neg_correction <- 0
    bt <- prior$bayestools_priors
    if (is.null(bt) && !is.null(cached_state)) bt <- cached_state$bayestools_priors
    neg_correction <- .compute_negative_mean_correction(metric, LSL, USL, bt)
    post_mean <- area * sum(mid_x * pdf_vals * dx) + neg_correction
    post_var <- area * sum((mid_x^2) * pdf_vals * dx) - post_mean^2
    post_sd <- sqrt(max(0, post_var))
  } else {
    area <- 1
    post_mean <- sum(mid_x * pdf_vals * dx)
    post_var <- sum((mid_x^2) * pdf_vals * dx) - post_mean^2
    post_sd <- sqrt(max(0, post_var))
  }

  # Quantiles (reconstruct CDF from grid, adjusted for mass at c <= 0)
  cdf_cond <- cumsum(pdf_vals * dx)
  cdf_adj <- (1 - area) + area * cdf_cond
  get_q <- function(q) {
    if (q <= 1 - area) return(0)
    mid_x[which.min(abs(cdf_adj - q))]
  }

  q2.5 <- get_q(0.025)
  q97.5 <- get_q(0.975)
  median_val <- get_q(0.5)

  # HDI (Highest Density Interval)
  sorted_idx <- order(pdf_vals, decreasing = TRUE)
  sorted_mass <- pdf_vals[sorted_idx] * dx[sorted_idx]
  cum_mass <- cumsum(sorted_mass)
  cutoff_idx <- which(cum_mass >= 0.95)[1]

  # Handle NA cutoff (e.g., if all pdf_vals are zero)
  if (is.na(cutoff_idx)) {
    cutoff_idx <- length(sorted_idx)
  }
  hdi_indices <- sorted_idx[1:cutoff_idx]

  # Override with analytic Inf where divergent
  if (divergence_info$mean_divergent) {
    post_mean <- Inf
    post_sd   <- Inf
  } else if (divergence_info$sd_divergent) {
    post_sd <- Inf
  }

  list(
    metric = metric,
    grid = data.frame(x = mid_x, density = pdf_vals),
    area = area,
    stats = c(Mean = post_mean, Median = median_val, SD = post_sd,
              Q2.5 = q2.5, Q97.5 = q97.5,
              HDI_Lo = min(mid_x[hdi_indices]),
              HDI_Hi = max(mid_x[hdi_indices])),
    divergence_info = divergence_info
  )
}


# ==============================================================================
# BayesTools Prior Conversion
# ==============================================================================

# ==============================================================================
# Analytic Moment Divergence Detection
# ==============================================================================

#' Extract the effective shape parameter alpha controlling tail behavior at sigma -> 0
#'
#' For a prior pi(sigma), alpha characterizes the density near zero:
#' pi(sigma) ~ sigma^(alpha - 1) as sigma -> 0.
#' Finite alpha indicates a potential divergence of \eqn{E[\sigma^{-k}]} for \eqn{k \ge \alpha}.
#' alpha = Inf means the prior is bounded away from zero or decays superexponentially.
#'
#' @param prior_sigma BayesTools prior object or string ("Jeffreys_sigma")
#' @return Numeric scalar: the effective alpha (possibly Inf)
#' @keywords internal
.extract_alpha_parameter <- function(prior_sigma) {
  if (identical(prior_sigma, "Jeffreys_sigma"))
    return(0)

  if (inherits(prior_sigma, "PriorConjugate"))
    return(prior_sigma$alpha0)

  if (!inherits(prior_sigma, "prior"))
    stop("prior_sigma must be a string or BayesTools::prior object")

  dist   <- prior_sigma[["distribution"]]
  params <- prior_sigma[["parameters"]]
  trunc  <- prior_sigma[["truncation"]]
  lower  <- trunc[["lower"]] %||% -Inf
  upper  <- trunc[["upper"]] %||%  Inf

  # Truncation bounded away from zero overrides everything
  if (is.finite(lower) && lower > 0)
    return(Inf)

  alpha <- switch(dist,
    "gamma"    = params[["shape"]],
    "exp"      = 1,
    "invgamma" = Inf,
    "lognormal"= Inf,
    "normal"   = Inf,
    "t"        = 1,
    "cauchy"   = 1,
    "point"    = Inf,
    "uniform"  = if (params[["a"]] <= 0) 1 else Inf,
    "beta"     = if (params[["a"]] <= 0 || (is.finite(lower) && lower <= 0)) params[["alpha"]] else Inf,
    Inf  # conservative default for unknown families
  )

  alpha
}

#' Check whether posterior moments of a capability metric diverge
#'
#' Uses the analytic decision rules from the divergence analysis: for metrics
#' scaling as \eqn{\sigma^{-1}} (Cp, Cpu, Cpl, Cpk), \eqn{E[C^k] < \infty} iff \eqn{\alpha > k}.
#' For metrics involving \eqn{\sqrt{\sigma^2 + (\mu - T)^2}} (Cpm, Cpc), the singularity
#' is regularised and \eqn{E[C^k] < \infty} iff \eqn{\alpha > k - 1}.
#'
#' @param metric One of "Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc"
#' @param alpha_sigma Effective alpha from .extract_alpha_parameter()
#' @param prior_sigma_label Human-readable label for the sigma prior (for messages)
#' @return List with mean_divergent, sd_divergent, alpha, reason
#' @keywords internal
.check_moment_divergence <- function(metric, alpha_sigma,
                                     prior_sigma_label = "sigma prior") {
  if (metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    mean_threshold <- 1
    var_threshold  <- 2
  } else if (metric %in% c("Cpm", "Cpc")) {
    mean_threshold <- 0
    var_threshold  <- 1
  } else {
    return(list(mean_divergent = FALSE, sd_divergent = FALSE,
                alpha = alpha_sigma, reason = NULL))
  }

  mean_div <- is.finite(alpha_sigma) && alpha_sigma <= mean_threshold
  sd_div   <- is.finite(alpha_sigma) && alpha_sigma <= var_threshold

  reason <- NULL
  if (mean_div) {
    reason <- sprintf(
      "%s has alpha=%.3g; %s requires alpha>%g for a finite mean (and alpha>%g for finite variance)",
      prior_sigma_label, alpha_sigma, metric, mean_threshold, var_threshold
    )
  } else if (sd_div) {
    reason <- sprintf(
      "%s has alpha=%.3g; %s requires alpha>%g for finite variance",
      prior_sigma_label, alpha_sigma, metric, var_threshold
    )
  }

  list(
    mean_divergent = mean_div,
    sd_divergent   = sd_div,
    alpha          = alpha_sigma,
    reason         = reason
  )
}

#' Build a human-readable label for a BayesTools prior
#' @param prior_sigma BayesTools prior object or string
#' @return Character string
#' @keywords internal
.prior_sigma_label <- function(prior_sigma) {
  if (identical(prior_sigma, "Jeffreys_sigma"))
    return("Jeffreys(sigma)")
  if (!inherits(prior_sigma, "prior"))
    return("unknown prior")

  dist   <- prior_sigma[["distribution"]]
  params <- prior_sigma[["parameters"]]
  param_str <- paste(vapply(params, function(p) format(p, digits = 3), character(1)), collapse = ", ")
  sprintf("%s(%s)", dist, param_str)
}

#' Convert BayesTools priors to integration prior format
#' @param prior_mu Prior for mu (string or BayesTools prior)
#' @param prior_sigma Prior for sigma (string or BayesTools prior)
#' @return List with $prior (integration prior object), $case (1-4), and $is_conjugate (logical)
#' @keywords internal

#' Classify mu prior as conjugate or non-conjugate
#' @param prior Prior specification for mu
#' @return List with is_conjugate, mu0, k0 (if conjugate), or log_dens_fn (if not)
#' @keywords internal
.classify_prior_mu <- function(prior) {
  # Jeffreys (flat) is conjugate with k0 = 0
  if (identical(prior, "Jeffreys_mu")) {
    return(list(is_conjugate = TRUE, mu0 = 0, k0 = 0))
  }

  # PriorConjugate carries the full NIG hyperparameters
  if (inherits(prior, "PriorConjugate")) {
    return(list(is_conjugate = TRUE, mu0 = prior$mu0, k0 = prior$k0))
  }

  # BayesTools Normal(mean, sd) prior on mu is sigma-INDEPENDENT: N(mean, sd^2).
  # The conjugate NIG model requires N(mu0, sigma^2/k0), which is sigma-DEPENDENT.
  # These are different models, so Normal mu priors are NOT conjugate here.

  # Not conjugate - build log-density function
  list(is_conjugate = FALSE, log_dens_fn = .make_prior_log_dens_fn(prior))
}

#' Classify sigma prior as conjugate or non-conjugate
#' @param prior Prior specification for sigma
#' @return List with is_conjugate, alpha0, beta0 (if conjugate), or log_dens_fn (if not)
#' @keywords internal
.classify_prior_sigma <- function(prior) {
  # Jeffreys (1/sigma) is conjugate with alpha0 = -0.5, beta0 = 0
  if (identical(prior, "Jeffreys_sigma")) {
    return(list(is_conjugate = TRUE, alpha0 = -0.5, beta0 = 0))
  }

  # PriorConjugate carries the full NIG hyperparameters
  if (inherits(prior, "PriorConjugate")) {
    return(list(is_conjugate = TRUE, alpha0 = prior$alpha0, beta0 = prior$beta0))
  }

  # InvGamma on sigma is NOT conjugate with the Normal likelihood because
  # BayesTools places InvGamma on sigma (not sigma^2). The conjugate model
  # requires InvGamma on sigma^2, which has a different functional form.
  # Gamma on sigma is also NOT conjugate.
  # All other distributions are non-conjugate.
  list(is_conjugate = FALSE, log_dens_fn = .make_prior_log_dens_fn(prior))
}

.bayestools_to_integration_prior <- function(prior_mu, prior_sigma) {

  # Classify each prior
  mu_info <- .classify_prior_mu(prior_mu)
  sigma_info <- .classify_prior_sigma(prior_sigma)
  bayestools_priors <- list(mu = prior_mu, sigma = prior_sigma)

  # Case 1: Full Conjugate (Normal/Jeffreys on mu AND InvGamma/Jeffreys on sigma)
  if (mu_info$is_conjugate && sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_conjugate(mu_info$mu0, mu_info$k0,
                                      sigma_info$alpha0, sigma_info$beta0),
      case = 1L,
      is_conjugate = TRUE
    ))
  }

  # Case 2: Semi-Conjugate Mu (conjugate mu, non-conjugate sigma)
  if (mu_info$is_conjugate && !sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_semi_mu(mu_info$mu0, mu_info$k0,
                                    sigma_info$log_dens_fn, bayestools_priors),
      case = 2L,
      is_conjugate = FALSE
    ))
  }

  # Case 3: Semi-Conjugate Sigma (non-conjugate mu, conjugate sigma)
  if (!mu_info$is_conjugate && sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_semi_sigma(sigma_info$alpha0, sigma_info$beta0,
                                       mu_info$log_dens_fn, bayestools_priors),
      case = 3L,
      is_conjugate = FALSE
    ))
  }

  # Case 4: Generic (both non-conjugate)
  log_dens_fn <- function(mu, sigma) {
    mu_info$log_dens_fn(mu) + sigma_info$log_dens_fn(sigma)
  }
  list(
    prior = create_prior_generic(log_dens_fn, bayestools_priors = bayestools_priors),
    case = 4L,
    is_conjugate = FALSE
  )
}

#' Extract a reasonable starting value and scale from a BayesTools prior
#' @keywords internal
.extract_prior_init <- function(bt_prior) {
  if (is.character(bt_prior)) {
    return(list(value = 1, scale = 10))
  }
  if (!inherits(bt_prior, "prior")) return(list(value = 1, scale = 1))

  dist   <- bt_prior[["distribution"]]
  params <- bt_prior[["parameters"]]
  trunc  <- bt_prior[["truncation"]]
  lower  <- trunc[["lower"]] %||% -Inf
  upper  <- trunc[["upper"]] %||%  Inf

  info <- switch(dist,
    "normal"    = list(value = params[["mean"]], scale = params[["sd"]]),
    "t"         = list(value = params[["location"]], scale = params[["scale"]]),
    "uniform"   = list(value = (params[["a"]] + params[["b"]]) / 2,
                       scale = (params[["b"]] - params[["a"]]) / 4),
    "gamma"     = { sh <- params[["shape"]]; rt <- params[["rate"]]
                    list(value = max((sh - 1) / rt, 0.5 / rt),
                         scale = sqrt(sh) / rt) },
    "invgamma"  = { a <- params[["shape"]]; b <- params[["scale"]]
                    list(value = b / (a + 1), scale = b / a) },
    "lognormal" = { ml <- params[["meanlog"]]; sl <- params[["sdlog"]]
                    list(value = exp(ml - sl^2),
                         scale = exp(ml) * sqrt(exp(sl^2) - 1)) },
    "exp"       = { rt <- params[["rate"]]
                    list(value = 1 / rt, scale = 1 / rt) },
    list(value = 1, scale = 1)
  )

  # Clamp value to truncation bounds
  if (is.finite(lower) && info$value < lower) {
    span <- if (is.finite(upper)) upper - lower else info$scale
    info$value <- lower + 0.1 * min(info$scale, span)
  }
  if (is.finite(upper) && info$value > upper) {
    span <- if (is.finite(lower)) upper - lower else info$scale
    info$value <- upper - 0.1 * min(info$scale, span)
  }
  info
}

#' Create a log-density function for a BayesTools prior
#' @param prior Prior specification (string or BayesTools prior object)
#' @return Function(x) returning vectorized log-density
#' @keywords internal
.make_prior_log_dens_fn <- function(prior) {

  # Handle string priors (Jeffreys)
  if (is.character(prior)) {
    if (prior == "Jeffreys_mu") {
      return(function(x) rep(0, length(x)))  # Improper flat prior
    } else if (prior == "Jeffreys_sigma") {
      return(function(x) ifelse(x <= 0, -Inf, -log(x)))
    }
    stop("Unknown string prior: ", prior)
  }

  # Handle BayesTools prior objects
  if (!inherits(prior, "prior")) {
    stop("prior must be a string or BayesTools::prior object")
  }

  dist <- prior[["distribution"]]
  params <- prior[["parameters"]]
  trunc <- prior[["truncation"]]

  # Check truncation bounds
  lower <- if (!is.null(trunc[["lower"]])) trunc[["lower"]] else -Inf
  upper <- if (!is.null(trunc[["upper"]])) trunc[["upper"]] else Inf

  # Pre-compute normalization constant for truncation
  log_norm <- 0
  if (is.finite(lower) || is.finite(upper)) {
    log_norm <- .compute_truncation_norm(dist, params, lower, upper)
  }


  # Construct efficient closure based on distribution
  base_dens_fn <- switch(dist,
    "point" = {
      loc <- params[["location"]]
      function(x) ifelse(x == loc, 0, -Inf)
    },
    "normal" = {
      mean <- params[["mean"]]
      sd <- params[["sd"]]
      function(x) stats::dnorm(x, mean = mean, sd = sd, log = TRUE)
    },
    "lognormal" = {
      meanlog <- params[["meanlog"]]
      sdlog <- params[["sdlog"]]
      function(x) ifelse(x <= 0, -Inf, stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE))
    },
    "t" = {
      loc <- params[["location"]]
      scale <- params[["scale"]]
      df <- params[["df"]]
      log_scale <- log(scale)
      function(x) {
        z <- (x - loc) / scale
        stats::dt(z, df = df, log = TRUE) - log_scale
      }
    },
    "gamma" = {
      shape <- params[["shape"]]
      rate <- params[["rate"]]
      # function(x) ifelse(x <= 0, -Inf, stats::dgamma(x, shape = shape, rate = rate, log = TRUE))
      function(x) stats::dgamma(x, shape = shape, rate = rate, log = TRUE)
    },
    "invgamma" = {
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      log_beta <- log(beta)
      lgamma_alpha <- lgamma(alpha)
      function(x) {
        ifelse(x <= 0, -Inf,
               alpha * log_beta - lgamma_alpha - (alpha + 1) * log(x) - beta / x)
      }
    },
    "uniform" = {
      a <- params[["a"]]
      b <- params[["b"]]
      log_width <- log(b - a)
      function(x) ifelse(x >= a & x <= b, -log_width, -Inf)
    },
    "beta" = {
      alpha <- params[["alpha"]]
      beta <- params[["beta"]]
      function(x) ifelse(x <= 0 | x >= 1, -Inf, stats::dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE))
    },
    "exp" = {
      rate <- params[["rate"]]
      function(x) ifelse(x < 0, -Inf, stats::dexp(x, rate = rate, log = TRUE))
    },
    stop("Unsupported prior distribution: ", dist)
  )

  # Return function that handles truncation and normalization
  function(x) {
    # Quick bounds check
    in_bounds <- x >= lower & x <= upper

    # If all valid, fast path
    if (all(in_bounds)) {
      return(base_dens_fn(x) - log_norm)
    }

    # If none valid, fast return
    if (!any(in_bounds)) {
      return(rep(-Inf, length(x)))
    }

    # Mixed case
    res <- rep(-Inf, length(x))
    res[in_bounds] <- base_dens_fn(x[in_bounds]) - log_norm
    res
  }
}


#' Compute log normalizing constant for truncated distribution
#' @keywords internal
.compute_truncation_norm <- function(dist, params, lower, upper) {

  cdf_fn <- switch(dist,
    "normal" = function(x) stats::pnorm(x, mean = params[["mean"]], sd = params[["sd"]]),
    "lognormal" = function(x) stats::plnorm(x, meanlog = params[["meanlog"]],
                                      sdlog = params[["sdlog"]]),
    "t" = function(x) stats::pt((x - params[["location"]]) / params[["scale"]],
                          df = params[["df"]]),
    "gamma" = function(x) stats::pgamma(x, shape = params[["shape"]], rate = params[["rate"]]),
    "invgamma" = function(x) {
      # CDF of inverse gamma
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      1 - stats::pgamma(beta / x, shape = alpha)
    },
    "uniform" = function(x) stats::punif(x, min = params[["a"]], max = params[["b"]]),
    "beta" = function(x) stats::pbeta(x, shape1 = params[["alpha"]], shape2 = params[["beta"]]),
    "exp" = function(x) stats::pexp(x, rate = params[["rate"]]),
    function(x) NA
  )

  p_upper <- if (is.infinite(upper)) 1 else cdf_fn(upper)
  p_lower <- if (is.infinite(lower)) 0 else cdf_fn(lower)

  log(p_upper - p_lower)
}

# ==============================================================================
# Integration Fit Function (called from bpc)
# ==============================================================================

#' Fit using numerical integration method
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param target Target value
#' @param prior_mu Prior for mu
#' @param prior_sigma Prior for sigma
#' @param sigma Number of standard deviations for capability metrics
#' @return List with metrics and integration results
#' @keywords internal
.bpc_fit_integration <- function(data, LSL, USL, target, prior_mu, prior_sigma,
                                   sigma = 3, sample_priors = FALSE,
                                   cached_state = NULL) {

  # If sampling from priors, ignore data (likelihood becomes flat/unity effectively)
  # Ideally we should pass sample_priors down, but clearing data works if the
  # functions handle empty data correctly.
  # However, the conjugate update uses n and sufficiency stats.
  # If we set data to empty, n=0, and the posterior parameters equal prior parameters.
  if (sample_priors) {
    data <- numeric(0)
  }

  # Convert priors to integration format (4-case classification)
  prior_info <- .bayestools_to_integration_prior(prior_mu, prior_sigma)
  prior <- prior_info$prior
  is_conjugate <- prior_info$is_conjugate
  case <- prior_info$case

  # Pre-compute state (sufficient statistics) unless already provided
  if (is.null(cached_state)) {
    n <- length(data)
    if (n > 0) {
      x_bar <- mean(data)
      sse <- sum((data - x_bar)^2)
    } else {
      x_bar <- 0
      sse <- 0
    }

    if (case == 4L) {
      cached_state <- precompute_generic_state(data, prior)
    } else {
      cached_state <- list(n = n, x_bar = x_bar, sse = sse)
    }
  }

  n <- cached_state$n %||% 0L

  # Analyze all metrics
  # Match order of MCMC results for consistency (Cp, Cpu, Cpl, Cpk, Cpc, Cpm)
  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")

  # Pre-compute analytic divergence info for each metric.
  # With data (n > 0), the likelihood provides superexponential decay at sigma = 0,
  # so all moments are finite regardless of the prior.
  is_prior_only <- n == 0
  alpha_sigma <- .extract_alpha_parameter(prior_sigma)
  sigma_label <- .prior_sigma_label(prior_sigma)
  divergence_map <- if (is_prior_only) {
    setNames(lapply(metrics, function(m)
      .check_moment_divergence(m, alpha_sigma, sigma_label)
    ), metrics)
  } else {
    no_div <- list(mean_divergent = FALSE, sd_divergent = FALSE,
                   alpha = alpha_sigma, reason = NULL)
    setNames(rep(list(no_div), length(metrics)), metrics)
  }

  # For prior-only mode with BayesTools priors, pre-generate shared samples
  # so all metrics use the same (mu, sigma) draws for consistency.
  mc_samples <- NULL
  has_bt_priors <- inherits(prior, "PriorGeneric") &&
                   !is.null(prior$bayestools_priors) &&
                   inherits(prior$bayestools_priors$sigma, "prior")
  if (is_prior_only && has_bt_priors) {
    bt <- prior$bayestools_priors
    n_mc <- 2000000L
    mc_samples <- list(
      mu  = BayesTools::rng(bt$mu, n_mc),
      sig = BayesTools::rng(bt$sigma, n_mc)
    )
  }

  results <- lapply(metrics, function(m) {
    analyze_capability_integration(data, LSL, USL, prior, metric = m,
                                    target = target, cached_state = cached_state,
                                    mc_samples = mc_samples,
                                    divergence_info = divergence_map[[m]])
  })
  names(results) <- metrics

  # Extract coefficients (posterior means)
  coefficients <- sapply(results, function(r) r$stats["Mean"])
  names(coefficients) <- metrics

  # Build metrics list compatible with bpc structure
  # For integration, we store stats rather than samples
  metrics_list <- lapply(metrics, function(m) {
    results[[m]]$stats
  })
  names(metrics_list) <- metrics
  class(metrics_list) <- "capability_metrics"
  attr(metrics_list, "LSL") <- LSL
  attr(metrics_list, "USL") <- USL
  attr(metrics_list, "target") <- target
  attr(metrics_list, "method") <- "integration"

  list(
    results = results,
    metrics = metrics_list,
    coefficients = coefficients,
    prior = prior,
    is_conjugate = is_conjugate,
    case = case,
    cached_state = cached_state,
    divergence = divergence_map
  )
}

