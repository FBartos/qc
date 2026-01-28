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
#' @return List with s_max_fn, mu_b_fn (scalar), and mu_b_fn_vec (vectorized)
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
    # Scalar version (for PriorGeneric)
    mu_b_fn = function(s) {
      if (metric == "Cp") return(c(-Inf, Inf))
      if (metric == "Cpk") return(c(LSL + 3 * c * s, USL - 3 * c * s))
      if (metric == "CpU") return(c(-Inf, USL - 3 * c * s))
      if (metric == "CpL") return(c(LSL + 3 * c * s, Inf))
      if (metric %in% c("Cpm", "Cpc")) {
        T_val <- if (metric == "Cpc") mid else target
        R <- tol / (6 * c)
        if (s >= R) return(c(0, -1))  # Empty interval
        w <- sqrt(R^2 - s^2)
        return(c(T_val - w, T_val + w))
      }
      stop("Unknown metric: ", metric)
    },
    # Vectorized version (for PriorConjugate): takes vector of sigma, returns list(lower, upper)
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
  x_bar <- mean(data)
  SS <- sum((data - x_bar)^2)

  # Posterior hyperparameters (Normal-Inverse-Gamma conjugate update)
  k_n <- prior$k0 + n
  mu_n <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
  alpha_n <- prior$alpha0 + n / 2
  beta_n <- prior$beta0 + 0.5 * SS + (prior$k0 * n * (x_bar - prior$mu0)^2) / (2 * k_n)
  df_p <- 2 * alpha_n

  # Pre-compute global h_max (chi-square mode density for numerical stability)
  y_mode <- max(df_p - 2, 1e-6)
  h_max_global <- dchisq(y_mode, df_p, log = TRUE)

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

      # log_diff_exp(pnorm(z_U, log.p=TRUE), pnorm(z_L, log.p=TRUE)) + dchisq(y, df_p, log=TRUE)
      log_prob <- log_diff_exp(pnorm(z_U, log.p = TRUE), pnorm(z_L, log.p = TRUE))
      result[valid] <- log_prob[valid] + dchisq(y[valid], df_p, log = TRUE)
      result
    }

    h_max <- h_max_global

    safe_integrand <- function(y) {
      vals <- exp(log_int(y) - h_max)
      vals[!is.finite(vals)] <- 0
      return(vals)
    }

    res <- integrate(safe_integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    return(exp(h_max + log(res)))
  }
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

    # Log posterior (scalar)
    log_post <- function(mu, sigma) {
      if (sigma <= 0) return(-Inf)
      -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2) +
        prior$log_dens(mu, sigma)
    }

    # Find MAP for integration bounds
    init_sd <- sqrt(sse / (n - 1))
    opt <- optim(c(x_bar, init_sd), function(p) -log_post(p[1], p[2]))
    map_mu <- opt$par[1]
    map_sig <- opt$par[2]
    h_max <- -opt$value

    # Tighter bounds: 5 sigma from MAP
    uni_s <- map_sig * 5

    # cubature 2D integration with vectorized interface
    int_2d <- function(s_lim, m_fn) {
      s_top <- if (is.infinite(s_lim)) uni_s else min(s_lim, uni_s)

      # Vectorized integrand: x is 2 x n matrix
      integrand <- function(x) {
        mu <- x[1, ]
        sigma <- x[2, ]

        result <- numeric(length(mu))
        for (i in seq_along(mu)) {
          if (sigma[i] <= 0) next
          mb <- m_fn(sigma[i])
          if (mu[i] < mb[1] || mu[i] > mb[2]) next
          result[i] <- exp(log_post(mu[i], sigma[i]) - h_max)
        }
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

    Z <- int_2d(Inf, function(s) c(-Inf, Inf))
  }

  # Return closure
  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    num <- int_2d(constr$s_max_fn(), constr$mu_b_fn)
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
  x_bar <- mean(data)
  sse <- sum((data - x_bar)^2)

  log_post <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2) +
      prior$log_dens(mu, sigma)
  }

  init_sd <- sqrt(sse / (n - 1))
  opt <- optim(c(x_bar, init_sd), function(p) -log_post(p[1], p[2]))
  map_mu <- opt$par[1]
  map_sig <- opt$par[2]
  h_max <- -opt$value

  # Tighter bounds: 5 sigma from MAP
  uni_s <- map_sig * 5

  # cubature 2D integration with vectorized interface
  int_2d <- function(s_lim, m_fn) {
    s_top <- if (is.infinite(s_lim)) uni_s else min(s_lim, uni_s)

    # Vectorized integrand: x is 2 x n matrix
    integrand <- function(x) {
      mu <- x[1, ]
      sigma <- x[2, ]

      result <- numeric(length(mu))
      for (i in seq_along(mu)) {
        if (sigma[i] <= 0) next
        mb <- m_fn(sigma[i])
        if (mu[i] < mb[1] || mu[i] > mb[2]) next
        result[i] <- exp(log_post(mu[i], sigma[i]) - h_max)
      }
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

  Z <- int_2d(Inf, function(s) c(-Inf, Inf))

  list(log_post = log_post, h_max = h_max, uni_s = uni_s, map_mu = map_mu,
       int_2d = int_2d, Z = Z)
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
  # Create solver function P(Index > c)
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

#' Analyze Capability with Automatic Grid Detection
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param n_grid Number of grid points for density evaluation
#' @param alpha_tail Probability mass to leave in tails for grid detection
#' @param cached_state Pre-computed state from precompute_generic_state
#' @return List with metric name, grid data.frame, and stats vector
#' @keywords internal
analyze_capability_integration <- function(data, LSL, USL, prior,
                                            metric = "Cpk", target = NULL,
                                            n_grid = 128, alpha_tail = 0.0001,
                                            cached_state = NULL) {

  # Create Solver Function: P(Index > c)
  if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
    S <- make_solver(data, LSL, USL, prior, metric, target, cached_state = cached_state)
  } else {
    S <- make_solver(data, LSL, USL, prior, metric, target)
  }

  # Heuristic search for bracket interval [0, max_c]
  max_c <- 3.0
  while (S(max_c) > alpha_tail) {
    max_c <- max_c * 2
    if (max_c > 100) break
  }

  # Find grid bounds using root finding (inverse CDF)
  get_quantile <- function(target_prob) {
    tryCatch({
      uniroot(function(c) S(c) - target_prob,
              interval = c(0, max_c),
              extendInt = "downX",
              tol = 1e-4)$root
    }, error = function(e) NA)
  }

  x_start <- get_quantile(1 - alpha_tail)
  x_end <- get_quantile(alpha_tail)

  if (is.na(x_start)) x_start <- 0
  if (is.na(x_end)) x_end <- max_c

  # Evaluate grid
  grid_x <- seq(x_start, x_end, length.out = n_grid)
  S_vals <- sapply(grid_x, S)

  # Compute PDF via finite differences
  pdf_vals <- -diff(S_vals) / diff(grid_x)
  mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2

  # Normalize area to 1.0
  area <- sum(pdf_vals * diff(grid_x))
  if (area > 0) pdf_vals <- pdf_vals / area

  # Compute statistics
  post_mean <- sum(mid_x * pdf_vals * diff(grid_x))
  post_var <- sum((mid_x^2) * pdf_vals * diff(grid_x)) - post_mean^2
  post_sd <- sqrt(max(0, post_var))

  # Quantiles (reconstruct CDF from grid)
  cdf_vals <- cumsum(pdf_vals * diff(grid_x))
  get_q <- function(q) mid_x[which.min(abs(cdf_vals - q))]

  q2.5 <- get_q(0.025)
  q97.5 <- get_q(0.975)
  median_val <- get_q(0.5)

  # HDI (Highest Density Interval)
  sorted_idx <- order(pdf_vals, decreasing = TRUE)
  sorted_mass <- pdf_vals[sorted_idx] * diff(grid_x)[1]
  cum_mass <- cumsum(sorted_mass)
  cutoff_idx <- which(cum_mass >= 0.95)[1]
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

  # Non-conjugate: build log_dens function from BayesTools priors
  log_dens_fn <- function(mu, sigma) {
    log_prior_mu <- .evaluate_prior_log_dens(prior_mu, mu, "mu")
    log_prior_sigma <- .evaluate_prior_log_dens(prior_sigma, sigma, "sigma")
    log_prior_mu + log_prior_sigma
  }

  list(
    prior = create_prior_generic(log_dens_fn),
    is_conjugate = FALSE
  )
}

#' Evaluate log-density for a BayesTools prior at a given value
#' @param prior Prior specification (string or BayesTools prior object)
#' @param x Value(s) at which to evaluate (can be vector)
#' @param param_name Parameter name ("mu" or "sigma") for Jeffreys handling
#' @return Log prior density at x (vectorized)
#' @keywords internal
.evaluate_prior_log_dens <- function(prior, x, param_name = NULL) {

  # Handle string priors (Jeffreys) - vectorized
  if (is.character(prior)) {
    if (prior == "Jeffreys_mu") {
      return(rep(0, length(x)))  # Improper flat prior on mu: log(1) = 0
    } else if (prior == "Jeffreys_sigma") {
      return(ifelse(x <= 0, -Inf, -log(x)))  # log(1/sigma)
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

  # Initialize output with -Inf for out-of-bounds values
  log_dens <- rep(-Inf, length(x))
  in_bounds <- x >= lower & x <= upper

  if (!any(in_bounds)) return(log_dens)

  x_valid <- x[in_bounds]

  # Evaluate log density based on distribution type (vectorized)
  log_dens_valid <- switch(dist,
    "point" = {
      ifelse(x_valid == params[["location"]], 0, -Inf)
    },
    "normal" = {
      dnorm(x_valid, mean = params[["mean"]], sd = params[["sd"]], log = TRUE)
    },
    "lognormal" = {
      ifelse(x_valid <= 0, -Inf,
             dlnorm(x_valid, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log = TRUE))
    },
    "t" = {
      # Location-scale t
      z <- (x_valid - params[["location"]]) / params[["scale"]]
      dt(z, df = params[["df"]], log = TRUE) - log(params[["scale"]])
    },
    "gamma" = {
      ifelse(x_valid <= 0, -Inf,
             dgamma(x_valid, shape = params[["shape"]], rate = params[["rate"]], log = TRUE))
    },
    "invgamma" = {
      # Inverse gamma: shape, scale parameterization
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      ifelse(x_valid <= 0, -Inf,
             alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x_valid) - beta / x_valid)
    },
    "uniform" = {
      a <- params[["a"]]
      b <- params[["b"]]
      ifelse(x_valid >= a & x_valid <= b, -log(b - a), -Inf)
    },
    "beta" = {
      ifelse(x_valid <= 0 | x_valid >= 1, -Inf,
             dbeta(x_valid, shape1 = params[["alpha"]], shape2 = params[["beta"]], log = TRUE))
    },
    "exp" = {
      ifelse(x_valid < 0, -Inf,
             dexp(x_valid, rate = params[["rate"]], log = TRUE))
    },
    stop("Unsupported prior distribution: ", dist)
  )

  # Adjust for truncation normalization (if truncated)
  if (is.finite(lower) || is.finite(upper)) {
    log_norm <- .compute_truncation_norm(dist, params, lower, upper)
    log_dens_valid <- log_dens_valid - log_norm
  }

  log_dens[in_bounds] <- log_dens_valid
  log_dens
}

#' Compute log normalizing constant for truncated distribution
#' @keywords internal
.compute_truncation_norm <- function(dist, params, lower, upper) {

  cdf_fn <- switch(dist,
    "normal" = function(x) pnorm(x, mean = params[["mean"]], sd = params[["sd"]]),
    "lognormal" = function(x) plnorm(x, meanlog = params[["meanlog"]],
                                      sdlog = params[["sdlog"]]),
    "t" = function(x) pt((x - params[["location"]]) / params[["scale"]],
                          df = params[["df"]]),
    "gamma" = function(x) pgamma(x, shape = params[["shape"]], rate = params[["rate"]]),
    "invgamma" = function(x) {
      # CDF of inverse gamma
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      1 - pgamma(beta / x, shape = alpha)
    },
    "uniform" = function(x) punif(x, min = params[["a"]], max = params[["b"]]),
    "beta" = function(x) pbeta(x, shape1 = params[["alpha"]], shape2 = params[["beta"]]),
    "exp" = function(x) pexp(x, rate = params[["rate"]]),
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
                                  sigma = 3) {

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

  list(
    results = results,
    metrics = metrics_list,
    coefficients = coefficients,
    prior = prior,
    is_conjugate = is_conjugate,
    cached_state = cached_state
  )
}
