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
#' @keywords internal
create_prior_conjugate <- function(mu0 = 0, k0 = 0, alpha0 = -0.5, beta0 = 0) {
  structure(list(mu0 = mu0, k0 = k0, alpha0 = alpha0, beta0 = beta0),
            class = "PriorConjugate")
}

#' Create a generic prior with custom log-density function
#' @param log_dens_fn Function(mu, sigma) returning log prior density
#' @return PriorGeneric object
#' @keywords internal
create_prior_generic <- function(log_dens_fn) {
  structure(list(log_dens = log_dens_fn), class = "PriorGeneric")
}

# ==============================================================================
# Metric Constraints
# ==============================================================================

#' Get constraint functions for a capability metric
#' @param metric One of "Cp", "Cpk", "Cpm", "Cpc", "CpU", "CpL"
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
      if (metric == "CpU") {
        return(list(lower = rep(-Inf, n), upper = USL - 3 * c * s))
      }
      if (metric == "CpL") {
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
#' @param metric One of "Cp", "Cpk", "CpU", "CpL", "Cpm", "Cpc"
#' @return Metric value(s) (same length as mu/sigma)
#' @keywords internal
compute_metric_value <- function(mu, sigma, LSL, USL, target, metric) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  switch(metric,
    "Cp"  = tol / (6 * sigma),
    "CpU" = (USL - mu) / (3 * sigma),
    "CpL" = (mu - LSL) / (3 * sigma),
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

#' @export
compute_metric_moments.PriorConjugate <- function(data, LSL, USL, prior,
                                                   metric = "Cpk", target = NULL,
                                                   use_analytic = TRUE,
                                                   cached_state = NULL) {
  n <- length(data)
  if (n > 0) {
    x_bar <- mean(data)
    SS <- sum((data - x_bar)^2)

    # Posterior hyperparameters
    k_n <- prior$k0 + n
    mu_n <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
    alpha_n <- prior$alpha0 + n / 2
    beta_n <- prior$beta0 + 0.5 * SS + (prior$k0 * n * (x_bar - prior$mu0)^2) / (2 * k_n)
  } else {
    # No data: posterior = prior
    k_n <- prior$k0
    mu_n <- prior$mu0
    alpha_n <- prior$alpha0
    beta_n <- prior$beta0
  }

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
  E_inv_sigma <- sqrt(2) * gamma((df_p + 1) / 2) / gamma(df_p / 2) / sqrt(2 * beta_n)
  E_inv_sigma2 <- df_p / (2 * beta_n)  # E[y/(2*beta)] = E[y]/(2*beta) = df/(2*beta)

  if (metric == "Cp") {
    # Cp = tol / (6*sigma)
    E1 <- (tol / 6) * E_inv_sigma
    # E[Cp^2] = (tol/6)^2 * E[1/sigma^2]
    E2 <- (tol / 6)^2 * E_inv_sigma2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "CpU") {
    # CpU = (USL - mu) / (3*sigma)
    # E[CpU] = E[(USL - mu)/sigma] / 3
    # mu|sigma ~ N(mu_n, sigma^2/k_n), so E[mu|sigma] = mu_n
    # E[(USL - mu)/sigma] = (USL - mu_n) * E[1/sigma]
    E1 <- (USL - mu_n) / 3 * E_inv_sigma

    # E[CpU^2] = E[(USL - mu)^2 / sigma^2] / 9
    # (USL - mu)^2 = (USL - mu_n)^2 - 2*(USL - mu_n)*(mu - mu_n) + (mu - mu_n)^2
    # E[(mu - mu_n)^2 | sigma] = sigma^2/k_n
    # E[(USL - mu)^2 / sigma^2] = (USL - mu_n)^2 * E[1/sigma^2] + E[1/k_n] = (USL - mu_n)^2 * E[1/sigma^2] + 1/k_n
    E2 <- ((USL - mu_n)^2 * E_inv_sigma2 + 1 / k_n) / 9
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "CpL") {
    # CpL = (mu - LSL) / (3*sigma)
    E1 <- (mu_n - LSL) / 3 * E_inv_sigma
    E2 <- ((mu_n - LSL)^2 * E_inv_sigma2 + 1 / k_n) / 9
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpk") {
    # Cpk = min(CpU, CpL) = (tol/2 - |mu - mid|) / (3*sigma)
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

    integrand_cpm <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      # E[1/sqrt(sigma^2 + (mu - T)^2) | sigma]
      # mu ~ N(mu_n, tau^2), so (mu - T) ~ N(delta_T, tau^2)
      # E[1/sqrt(sigma^2 + X^2)] where X ~ N(delta_T, tau^2)
      # Use Gauss-Hermite style quadrature over standardized variable
      z_pts <- c(-3, -2, -1, 0, 1, 2, 3)
      w_pts <- stats::dnorm(z_pts)
      w_pts <- w_pts / sum(w_pts)
      x_pts <- delta_T + tau * z_pts
      inner <- sum(w_pts / sqrt(sigma^2 + x_pts^2))
      inner * stats::dchisq(y, df_p)
    }
    E_inv_sqrt <- stats::integrate(Vectorize(integrand_cpm), 1e-6, y_upper,
                            rel.tol = 1e-4, subdivisions = 200)$value
    E1 <- (tol / 6) * E_inv_sqrt

    integrand_cpm2 <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      z_pts <- c(-3, -2, -1, 0, 1, 2, 3)
      w_pts <- stats::dnorm(z_pts)
      w_pts <- w_pts / sum(w_pts)
      x_pts <- delta_T + tau * z_pts
      inner <- sum(w_pts / (sigma^2 + x_pts^2))
      inner * stats::dchisq(y, df_p)
    }
    E_inv <- stats::integrate(Vectorize(integrand_cpm2), 1e-6, y_upper,
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
                                        target = NULL, ...) {
  n <- length(data)
  if (n > 0) {
    x_bar <- mean(data)
    SS <- sum((data - x_bar)^2)

    # Posterior hyperparameters (Normal-Inverse-Gamma conjugate update)
    k_n <- prior$k0 + n
    mu_n <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
    alpha_n <- prior$alpha0 + n / 2
    beta_n <- prior$beta0 + 0.5 * SS + (prior$k0 * n * (x_bar - prior$mu0)^2) / (2 * k_n)
  } else {
    k_n <- prior$k0
    mu_n <- prior$mu0
    alpha_n <- prior$alpha0
    beta_n <- prior$beta0
  }
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
#' @param metric Capability metric (currently "Cpk", "CpU", "CpL", "Cp" supported)
#' @param target Target value for Cpm/Cpc
#' @param ... Additional arguments
#' @return Function pdf(c) returning the marginal posterior density at c
#' @keywords internal
make_density_solver <- function(data, LSL, USL, prior, metric = "Cpk", target = NULL, ...) {
  UseMethod("make_density_solver", prior)
}

#' @export
make_density_solver.PriorConjugate <- function(data, LSL, USL, prior,
                                                metric = "Cpk", target = NULL, ...) {
  n <- length(data)
  if (n > 0) {
    x_bar <- mean(data)
    SS <- sum((data - x_bar)^2)

    # Posterior hyperparameters (Normal-Inverse-Gamma conjugate update)
    k_n <- prior$k0 + n
    mu_n <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
    alpha_n <- prior$alpha0 + n / 2
    beta_n <- prior$beta0 + 0.5 * SS + (prior$k0 * n * (x_bar - prior$mu0)^2) / (2 * k_n)
  } else {
    k_n <- prior$k0
    mu_n <- prior$mu0
    alpha_n <- prior$alpha0
    beta_n <- prior$beta0
  }

  M <- (LSL + USL) / 2  # Midpoint
  tol <- USL - LSL      # Tolerance

  # Log-normalizing constant for sigma marginal
  # The marginal of sigma from NIG is:
  # p(sigma) = 2 * beta^alpha / Gamma(alpha) * sigma^(-2*alpha - 1) * exp(-beta/sigma^2)
  # So the normalizing constant 1/Z = 2 * beta^alpha / Gamma(alpha)
  # And log(Z) = lgamma(alpha) - log(2) - alpha * log(beta)
  log_Z_sigma <- lgamma(alpha_n) - log(2) - alpha_n * log(beta_n)

  # Return pdf(c) function for Cpk
  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_max <- tol / (6 * c)
      if (sigma_max <= 0) return(0)

      # Contour points: mu_L = LSL + 3*c*sigma (CpL = c), mu_U = USL - 3*c*sigma (CpU = c)
      # Both contours valid for sigma < sigma_max
      # At sigma = sigma_max, both contours meet at mu = M

      integrand <- function(sigma) {
        # Guard against boundary issues
        valid <- sigma > 0 & sigma < sigma_max
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        sd_mu <- sigma_v / sqrt(k_n)

        # Contour points
        mu_L <- LSL + 3 * c * sigma_v  # CpL = c contour
        mu_U <- USL - 3 * c * sigma_v  # CpU = c contour

        # Log marginal posterior of sigma (unnormalized)
        # p(sigma) ∝ sigma^(-2*alpha_n - 1) * exp(-beta_n / sigma^2)
        log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_v) - beta_n / sigma_v^2

        # Conditional densities p(mu|sigma) at contour points
        log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
        log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)

        # Jacobian |dmu/dc| = 3*sigma
        log_jacobian <- log(3) + log(sigma_v)

        # Contributions from both contours
        contrib_L <- exp(log_p_sigma + log_p_mu_L + log_jacobian - log_Z_sigma)
        contrib_U <- exp(log_p_sigma + log_p_mu_U + log_jacobian - log_Z_sigma)

        # Handle NaN/Inf
        contrib_L[!is.finite(contrib_L)] <- 0
        contrib_U[!is.finite(contrib_U)] <- 0

        result[valid] <- contrib_L + contrib_U
        result
      }

      stats::integrate(integrand, 0, sigma_max, rel.tol = 1e-6, subdivisions = 200)$value
    })
  }

  # CpU: single contour mu_U = USL - 3*c*sigma
  if (metric == "CpU") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        sd_mu <- sigma_v / sqrt(k_n)
        mu_U <- USL - 3 * c * sigma_v

        log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_v) - beta_n / sigma_v^2
        log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)
        log_jacobian <- log(3) + log(sigma_v)

        contrib <- exp(log_p_sigma + log_p_mu_U + log_jacobian - log_Z_sigma)
        contrib[!is.finite(contrib)] <- 0
        result[valid] <- contrib
        result
      }

      # Integrate over all sigma > 0 (use large upper bound)
      sigma_upper <- sqrt(beta_n / alpha_n) * 10  # ~10x posterior mode
      stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-6, subdivisions = 200)$value
    })
  }

  # CpL: single contour mu_L = LSL + 3*c*sigma
  if (metric == "CpL") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        sd_mu <- sigma_v / sqrt(k_n)
        mu_L <- LSL + 3 * c * sigma_v

        log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_v) - beta_n / sigma_v^2
        log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
        log_jacobian <- log(3) + log(sigma_v)

        contrib <- exp(log_p_sigma + log_p_mu_L + log_jacobian - log_Z_sigma)
        contrib[!is.finite(contrib)] <- 0
        result[valid] <- contrib
        result
      }

      sigma_upper <- sqrt(beta_n / alpha_n) * 10
      stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-6, subdivisions = 200)$value
    })
  }

  # Cp: single contour sigma = tol / (6*c)
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_c <- tol / (6 * c)
      if (sigma_c <= 0) return(0)

      # For Cp = tol/(6*sigma), at c the contour is sigma = tol/(6*c)
      # Jacobian: |dsigma/dc| = tol / (6 * c^2)

      log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_c) - beta_n / sigma_c^2
      log_jacobian <- log(tol / 6) - 2 * log(c)

      exp(log_p_sigma + log_jacobian - log_Z_sigma)
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
        sd_mu <- sigma_v / sqrt(k_n)

        # Log marginal posterior of sigma
        log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_v) - beta_n / sigma_v^2

        # Conditional densities p(mu|sigma)
        log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)

        # Combined density
        log_contrib <- log_p_sigma + log_p_mu + log_jacobian - log_Z_sigma

        result[valid] <- exp(log_contrib)
        result
      }

      stats::integrate(integrand, 0, pi, rel.tol = 1e-5, subdivisions = 200)$value
    })
  }

  stop("Unsupported metric for density solver: ", metric)
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

  # CpU: single contour
  if (metric == "CpU") {
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

      sigma_upper <- uni_s
      stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-5, subdivisions = 200)$value
    })
  }

  # CpL: single contour
  if (metric == "CpL") {
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

      sigma_upper <- uni_s
      stats::integrate(integrand, 0, sigma_upper, rel.tol = 1e-5, subdivisions = 200)$value
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

      # Use dynamic bounds based on sigma_c to capture the likelihood peak
      # p(mu|sigma) is roughly N(x_bar, sigma^2/n)
      # We use +/- 15 SDs to be safe
      sd_mu_cond <- sigma_c / sqrt(n)
      mu_width <- 15 * sd_mu_cond

      mu_lower <- x_bar - mu_width
      mu_upper <- x_bar + mu_width

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
    # No data: use default initial values
    x_bar <- 0
    sse <- 0
    init_sd <- 1
    init_mu <- 0
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

  opt <- stats::optim(c(init_mu, init_sd), function(p) -log_post(p[1], p[2]))
  map_mu <- opt$par[1]
  map_sig <- opt$par[2]
  h_max <- -opt$value

  # Tighter bounds: 5 sigma from MAP
  # If n=0, map_sig is determined by prior. If prior is flat or wide, 5 sigma might be too small
  # or too large depending on the prior.
  # For safety, ensure uni_s is at least a reasonable minimum and scale up if n=0
  uni_s <- map_sig * (if (n > 0) 5 else 10)
  uni_s <- max(uni_s, 20.0) # Ensure at least some width

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

  list(log_post = log_post, log_post_vec = log_post_vec, h_max = h_max,
       uni_s = uni_s, map_mu = map_mu, int_2d = int_2d, Z = Z, n = n,
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
  density_metrics <- c("Cpk", "Cp", "CpU", "CpL", "Cpm", "Cpc")
  can_use_density <- (inherits(prior, "PriorConjugate") || inherits(prior, "PriorGeneric")) &&
                     metric %in% density_metrics

  if (can_use_density) {
    # Use density solver PDF and integrate over bounds
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                   cached_state = cached_state)
    pdf_vec <- Vectorize(pdf_fn)

    # Integrate PDF from bounds[1] to bounds[2]
    prob <- stats::integrate(pdf_vec, min(bounds), max(bounds), rel.tol = 1e-4)$value
    return(prob)

  } else {
    # Fallback: Survival function S(c) = P(Index > c)
    if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
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
                                            use_density_solver = TRUE) {

  # Compute posterior moments to determine grid bounds
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
       S_fn <- if (!is.null(cached_state) && inherits(prior, "PriorGeneric")) {
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

  # Ensure valid grid (minimum width)
  if (x_end <= x_start || !is.finite(x_end)) {
    x_end <- max(x_start + 3, 3)
  }

  # Choose method: density solver (direct PDF) vs survival function + finite diff
  # Density solver works for both PriorConjugate and PriorGeneric with supported metrics
  density_metrics <- c("Cpk", "Cp", "CpU", "CpL", "Cpm", "Cpc")
  can_use_density <- use_density_solver &&
                     (inherits(prior, "PriorConjugate") || inherits(prior, "PriorGeneric")) &&
                     metric %in% density_metrics

  if (can_use_density) {
    # Direct PDF computation via contour integration
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                   cached_state = cached_state)

    # Evaluate on grid (use midpoints for consistency)
    grid_x <- seq(x_start, x_end, length.out = n_grid)
    mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2
    pdf_vals <- sapply(mid_x, pdf_fn)

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
  post_mean <- sum(mid_x * pdf_vals * dx)
  post_var <- sum((mid_x^2) * pdf_vals * dx) - post_mean^2
  post_sd <- sqrt(max(0, post_var))

  # Quantiles (reconstruct CDF from grid)
  cdf_vals <- cumsum(pdf_vals * dx)
  get_q <- function(q) mid_x[which.min(abs(cdf_vals - q))]

  q2.5 <- get_q(0.025)
  q97.5 <- get_q(0.975)
  median_val <- get_q(0.5)

  # HDI (Highest Density Interval)
  sorted_idx <- order(pdf_vals, decreasing = TRUE)
  sorted_mass <- pdf_vals[sorted_idx] * dx[1]
  cum_mass <- cumsum(sorted_mass)
  cutoff_idx <- which(cum_mass >= 0.95)[1]

  # Handle NA cutoff (e.g., if all pdf_vals are zero)
  if (is.na(cutoff_idx)) {
    cutoff_idx <- length(sorted_idx)
  }
  hdi_indices <- sorted_idx[1:cutoff_idx]

  list(
    metric = metric,
    grid = data.frame(x = mid_x, density = pdf_vals),
    stats = c(Mean = post_mean, Median = median_val, SD = post_sd,
              Q2.5 = q2.5, Q97.5 = q97.5,
              HDI_Lo = min(mid_x[hdi_indices]),
              HDI_Hi = max(mid_x[hdi_indices]))
  )
}


# ==============================================================================
# BayesTools Prior Conversion
# ==============================================================================

#' Convert BayesTools priors to integration prior format
#' @param prior_mu Prior for mu (string or BayesTools prior)
#' @param prior_sigma Prior for sigma (string or BayesTools prior)
#' @return List with $prior (integration prior object) and $is_conjugate (logical)
#' @keywords internal

.bayestools_to_integration_prior <- function(prior_mu, prior_sigma) {

  # Check for conjugate case: Jeffreys priors
  if (identical(prior_mu, "Jeffreys_mu") && identical(prior_sigma, "Jeffreys_sigma")) {
    return(list(
      prior = create_prior_conjugate(),
      is_conjugate = TRUE
    ))
  }

  # Compile efficient log-density functions for mu and sigma
  # This avoids repeated switch/lookup overhead on every evaluation
  log_dens_mu_fn <- .make_prior_log_dens_fn(prior_mu)
  log_dens_sigma_fn <- .make_prior_log_dens_fn(prior_sigma)

  # Non-conjugate: combine the compiled functions
  log_dens_fn <- function(mu, sigma) {
    # If param_name was needed for Jeffreys, it's handled inside the factory now
    log_dens_mu_fn(mu) + log_dens_sigma_fn(sigma)
  }

  list(
    prior = create_prior_generic(log_dens_fn),
    is_conjugate = FALSE
  )
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
                                  sigma = 3, sample_priors = FALSE) {

  # If sampling from priors, ignore data (likelihood becomes flat/unity effectively)
  # Ideally we should pass sample_priors down, but clearing data works if the
  # functions handle empty data correctly.
  # However, the conjugate update uses n and sufficiency stats.
  # If we set data to empty, n=0, and the posterior parameters equal prior parameters.
  if (sample_priors) {
    data <- numeric(0)
  }

  # Convert priors to integration format
  prior_info <- .bayestools_to_integration_prior(prior_mu, prior_sigma)
  prior <- prior_info$prior
  is_conjugate <- prior_info$is_conjugate

  # Pre-compute state for non-conjugate priors
  cached_state <- NULL
  if (!is_conjugate) {
    cached_state <- precompute_generic_state(data, prior)
  }

  # Analyze all metrics
  metrics <- c("Cp", "Cpk", "Cpm", "CpU", "CpL", "Cpc")
  results <- lapply(metrics, function(m) {
    analyze_capability_integration(data, LSL, USL, prior, metric = m,
                                    target = target, cached_state = cached_state)
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
    cached_state = cached_state
  )
}
