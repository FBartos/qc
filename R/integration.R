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

#' Create a unit information prior for a normal model
#'
#' Constructs a Normal-Inverse-Gamma (NIG) conjugate prior that carries
#' approximately one unit of Fisher information.  The prior parameters are
#' derived from the observed data:
#' \itemize{
#'   \item \eqn{\mu \mid \sigma^2 \sim \mathrm{Normal}(\bar{x},\, \sigma^2 / 1)}
#'   \item \eqn{\sigma^2 \sim \mathrm{Inv-Gamma}(1/2,\, s^2/2)}
#' }
#' which corresponds to \code{create_prior_conjugate(mu0 = xbar, k0 = 1,
#' alpha0 = 0.5, beta0 = s2 / 2)}.
#'
#' The returned \code{PriorConjugate} object can be passed as \code{prior_mu}
#' to \code{\link{bpc}} when \code{method = "integration"}.
#'
#' @param x Numeric vector of observations.  \code{NA} values are removed.
#'   At least 2 finite observations are required.
#' @return A \code{PriorConjugate} object.
#' @export
create_prior_unit_information <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L)
    stop("At least 2 finite observations are required to compute a unit information prior.")
  xbar <- mean(x)
  s2   <- stats::var(x)
  create_prior_conjugate(mu0 = xbar, k0 = 1, alpha0 = 0.5, beta0 = s2 / 2)
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

.is_improper_conjugate_posterior <- function(k_n, alpha_n, beta_n) {
  !is.finite(k_n) || !is.finite(alpha_n) || !is.finite(beta_n) ||
    k_n <= 0 || alpha_n <= 0 || beta_n < 0
}

.is_degenerate_conjugate_posterior <- function(beta_n) {
  is.finite(beta_n) && beta_n == 0
}

.scalar_almost_equal <- function(x, y) {
  tol <- sqrt(.Machine$double.eps) * max(1, abs(x), abs(y))
  abs(x - y) <= tol
}

.degenerate_conjugate_metric_distribution <- function(mu_n, k_n, LSL, USL,
                                                      target, metric,
                                                      sigma_level = 3) {
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  if (!is.finite(k_n) || k_n <= 0) {
    return(NULL)
  }

  boundary_sd <- 1 / (sigma_level * sqrt(k_n))

  if (metric == "Cp") {
    return(list(type = "pos_inf"))
  }

  if (metric == "Cpu") {
    if (.scalar_almost_equal(mu_n, USL)) {
      return(list(type = "normal", mean = 0, sd = boundary_sd))
    }
    return(list(type = if (mu_n < USL) "pos_inf" else "neg_inf"))
  }

  if (metric == "Cpl") {
    if (.scalar_almost_equal(mu_n, LSL)) {
      return(list(type = "normal", mean = 0, sd = boundary_sd))
    }
    return(list(type = if (mu_n > LSL) "pos_inf" else "neg_inf"))
  }

  if (metric == "Cpk") {
    if (.scalar_almost_equal(mu_n, LSL) || .scalar_almost_equal(mu_n, USL)) {
      return(list(type = "normal", mean = 0, sd = boundary_sd))
    }
    if (mu_n > LSL && mu_n < USL) {
      return(list(type = "pos_inf"))
    }
    return(list(type = "neg_inf"))
  }

  if (metric == "Cpm") {
    delta <- abs(mu_n - target)
    if (.scalar_almost_equal(delta, 0)) {
      return(list(type = "pos_inf"))
    }
    value <- (USL - LSL) / ((2 * sigma_level) * delta)
    return(list(type = "point", value = value))
  }

  if (metric == "Cpc") {
    delta <- abs(mu_n - target)
    if (.scalar_almost_equal(delta, 0)) {
      return(list(type = "pos_inf"))
    }
    value <- (USL - LSL) / ((2 * sigma_level) * sqrt(pi / 2) * delta)
    return(list(type = "point", value = value))
  }

  NULL
}

.degenerate_metric_moments <- function(dist) {
  switch(dist$type,
    "point" = list(mean = dist$value, sd = 0),
    "normal" = list(mean = dist$mean, sd = dist$sd),
    "pos_inf" = list(mean = Inf, sd = Inf),
    "neg_inf" = list(mean = -Inf, sd = Inf),
    stop("Unknown degenerate metric distribution type: ", dist$type)
  )
}

.degenerate_metric_quantiles <- function(dist, probs) {
  probs <- pmin(pmax(probs, 0), 1)
  switch(dist$type,
    "point" = rep(dist$value, length(probs)),
    "normal" = stats::qnorm(probs, mean = dist$mean, sd = dist$sd),
    "pos_inf" = rep(Inf, length(probs)),
    "neg_inf" = rep(-Inf, length(probs)),
    stop("Unknown degenerate metric distribution type: ", dist$type)
  )
}

.degenerate_metric_interval <- function(dist, ci, ci_level) {
  ci_level <- max(min(ci_level, 1), 0)
  h <- (1 - ci_level) / 2

  switch(dist$type,
    "point" = rep(dist$value, 2),
    "normal" = {
      if (ci == "central" || ci == "HPD") {
        stats::qnorm(c(h, 1 - h), mean = dist$mean, sd = dist$sd)
      } else {
        stop("Unknown ci for degenerate metric distribution.")
      }
    },
    "pos_inf" = rep(Inf, 2),
    "neg_inf" = rep(-Inf, 2),
    stop("Unknown degenerate metric distribution type: ", dist$type)
  )
}

.degenerate_metric_solver <- function(dist) {
  force(dist)

  function(c) {
    c <- as.numeric(c)
    switch(dist$type,
      "point" = as.numeric(c < dist$value),
      "normal" = stats::pnorm(c, mean = dist$mean, sd = dist$sd,
                               lower.tail = FALSE),
      "pos_inf" = as.numeric(c < Inf),
      "neg_inf" = rep(0, length(c)),
      stop("Unknown degenerate metric distribution type: ", dist$type)
    )
  }
}

.degenerate_metric_density <- function(dist) {
  force(dist)

  function(c) {
    c <- as.numeric(c)
    switch(dist$type,
      "normal" = stats::dnorm(c, mean = dist$mean, sd = dist$sd),
      "point" = rep(0, length(c)),
      "pos_inf" = rep(0, length(c)),
      "neg_inf" = rep(0, length(c)),
      stop("Unknown degenerate metric distribution type: ", dist$type)
    )
  }
}

.degenerate_metric_prob <- function(dist, bounds) {
  lower <- min(bounds)
  upper <- max(bounds)

  if (!is.finite(lower) && !is.finite(upper)) {
    return(1)
  }

  switch(dist$type,
    "point" = as.numeric(lower < dist$value && dist$value < upper),
    "normal" = {
      stats::pnorm(upper, mean = dist$mean, sd = dist$sd) -
        stats::pnorm(lower, mean = dist$mean, sd = dist$sd)
    },
    "pos_inf" = as.numeric(is.infinite(upper) && upper > 0),
    "neg_inf" = as.numeric(is.infinite(lower) && lower < 0),
    stop("Unknown degenerate metric distribution type: ", dist$type)
  )
}

.degenerate_metric_grid <- function(dist, n_grid = 512L, metric_can_be_negative = FALSE) {
  n_grid <- max(as.integer(n_grid), 64L)

  if (dist$type == "normal") {
    probs <- seq(0.001, 0.999, length.out = n_grid)
    x <- stats::qnorm(probs, mean = dist$mean, sd = dist$sd)
    density <- stats::dnorm(x, mean = dist$mean, sd = dist$sd)
    return(data.frame(x = x, density = density))
  }

  if (dist$type == "point") {
    plot_sd <- max(1e-6, 0.01 * max(1, abs(dist$value)))
    x <- seq(dist$value - 4 * plot_sd, dist$value + 4 * plot_sd, length.out = n_grid)
    density <- stats::dnorm(x, mean = dist$value, sd = plot_sd)
    return(data.frame(x = x, density = density))
  }

  if (metric_can_be_negative || identical(dist$type, "neg_inf")) {
    x <- seq(-1, 1, length.out = n_grid)
  } else {
    x <- seq(0, 1, length.out = n_grid)
  }

  data.frame(x = x, density = rep(0, length(x)))
}

.analyze_degenerate_metric_distribution <- function(metric, dist, n_grid,
                                                    divergence_info = NULL,
                                                    metric_can_be_negative = FALSE) {
  if (is.null(divergence_info)) {
    divergence_info <- list(mean_divergent = FALSE, sd_divergent = FALSE,
                            alpha = Inf, reason = NULL)
  }

  moments <- .degenerate_metric_moments(dist)
  central <- .degenerate_metric_interval(dist, "central", 0.95)
  hdi <- .degenerate_metric_interval(dist, "HPD", 0.95)
  median_val <- .degenerate_metric_quantiles(dist, 0.5)

  list(
    metric = metric,
    grid = .degenerate_metric_grid(dist, n_grid,
                                   metric_can_be_negative = metric_can_be_negative),
    area = if (dist$type %in% c("point", "normal")) 1 else 0,
    stats = c(Mean = moments$mean, Median = median_val, SD = moments$sd,
              Q2.5 = central[1], Q97.5 = central[2],
              HDI_Lo = hdi[1], HDI_Hi = hdi[2]),
    divergence_info = divergence_info,
    degenerate = dist
  )
}

.semi_mu_posterior <- function(prior, n, x_bar, sse) {
  k_n   <- prior$k0 + n
  mu_n  <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
  sse_n <- sse + prior$k0 * n * (x_bar - prior$mu0)^2 / k_n
  list(k_n = k_n, mu_n = mu_n, sse_n = sse_n)
}

.cpc_lookup <- local({
  z <- seq(0, 12, length.out = 4097L)
  g <- sqrt(2 / pi) * exp(-0.5 * z^2) + z * (2 * stats::pnorm(z) - 1)
  r <- z / g
  r[1L] <- 0
  list(
    z = z,
    g = g,
    r = r,
    g0 = g[1L],
    g_max = g[length(g)],
    r_max = r[length(r)]
  )
})

.cpc_g <- function(z) {
  sqrt(2 / pi) * exp(-0.5 * z^2) + z * (2 * stats::pnorm(z) - 1)
}

.cpc_E_abs_dev_normal <- function(mu, sigma, target) {
  delta <- abs(mu - target)
  z <- delta / sigma
  sigma * sqrt(2 / pi) * exp(-0.5 * z^2) +
    delta * (1 - 2 * stats::pnorm(-z))
}

.cpc_mu_width_from_sigma <- function(sigma, c, tol, sigma_level = 3) {
  A <- tol / ((2 * sigma_level) * sqrt(pi / 2))
  K <- A / (c * sigma)
  width <- rep(NA_real_, length(K))
  feasible <- is.finite(K) & K >= .cpc_lookup$g0
  if (!any(feasible)) return(width)

  z <- numeric(sum(feasible))
  Kf <- K[feasible]
  use_interp <- Kf <= .cpc_lookup$g_max
  if (any(use_interp)) {
    z[use_interp] <- stats::approx(
      x = .cpc_lookup$g,
      y = .cpc_lookup$z,
      xout = Kf[use_interp],
      ties = "ordered",
      rule = 2
    )$y
  }
  if (any(!use_interp)) {
    # For z > 8, g(z) is numerically indistinguishable from z.
    z[!use_interp] <- Kf[!use_interp]
  }

  width[feasible] <- sigma[feasible] * z
  width
}

.cpc_sigma_limit_from_delta <- function(delta, c, tol, sigma_level = 3) {
  delta <- abs(delta)
  sigma_max <- tol / ((2 * sigma_level) * c)
  A <- tol / ((2 * sigma_level) * sqrt(pi / 2))

  limit <- numeric(length(delta))
  zero_delta <- delta == 0
  limit[zero_delta] <- sigma_max

  q <- delta * c / A
  feasible <- !zero_delta & is.finite(q) & q < 1
  if (any(feasible)) {
    qf <- q[feasible]
    z <- stats::approx(
      x = .cpc_lookup$r,
      y = .cpc_lookup$z,
      xout = qf,
      ties = "ordered",
      rule = 2
    )$y
    limit[feasible] <- delta[feasible] / z
  }

  pmin(limit, sigma_max)
}

.cpc_contour_from_z <- function(z, c, tol, target, sigma_level = 3) {
  g_z <- .cpc_g(z)
  sigma <- tol / ((2 * sigma_level) * sqrt(pi / 2) * c * g_z)
  delta <- sigma * z
  list(
    sigma = sigma,
    mu_lower = target - delta,
    mu_upper = target + delta,
    log_jacobian = 2 * log(sigma) - log(c)
  )
}

.metric_can_be_negative <- function(metric) {
  metric %in% c("Cpu", "Cpl", "Cpk")
}

.metric_sigma_limit_from_mu <- function(metric, mu, c, LSL, USL, target,
                                        sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  switch(
    metric,
    "Cp" = rep(tol / ((2 * sigma_level) * c), length(mu)),
    "Cpu" = pmax(0, (USL - mu) / (sigma_level * c)),
    "Cpl" = pmax(0, (mu - LSL) / (sigma_level * c)),
    "Cpk" = pmax(0, pmin(USL - mu, mu - LSL) / (sigma_level * c)),
    "Cpm" = {
      sigma_cap <- tol / ((2 * sigma_level) * c)
      delta <- abs(mu - target)
      limit <- sqrt(pmax(0, sigma_cap^2 - delta^2))
      limit[delta >= sigma_cap] <- 0
      limit
    },
    "Cpc" = .cpc_sigma_limit_from_delta(abs(mu - target), c, tol, sigma_level),
    stop("Unknown metric: ", metric)
  )
}

.metric_sigma_region_from_mu <- function(metric, mu, c, LSL, USL, target,
                                         sigma_level = 3) {
  lower <- rep(0, length(mu))
  upper <- rep(Inf, length(mu))

  if (metric %in% c("Cp", "Cpm", "Cpc")) {
    if (c <= 0) {
      return(list(lower = lower, upper = upper))
    }
    upper <- .metric_sigma_limit_from_mu(
      metric, mu, c, LSL, USL, target,
      sigma_level = sigma_level
    )
    upper[!is.finite(upper) | upper <= 0] <- -Inf
    return(list(lower = lower, upper = upper))
  }

  numerator <- switch(
    metric,
    "Cpu" = USL - mu,
    "Cpl" = mu - LSL,
    "Cpk" = pmin(USL - mu, mu - LSL),
    stop("Unknown metric: ", metric)
  )

  if (c > 0) {
    upper <- numerator / (sigma_level * c)
    upper[numerator <= 0 | !is.finite(upper) | upper <= 0] <- -Inf
    return(list(lower = lower, upper = upper))
  }

  if (c < 0) {
    full_region <- numerator >= 0
    lower[full_region] <- 0
    upper[full_region] <- Inf

    needs_lower_bound <- !full_region
    lower[needs_lower_bound] <- numerator[needs_lower_bound] / (sigma_level * c)
    lower[needs_lower_bound & (!is.finite(lower) | lower <= 0)] <- Inf
    upper[needs_lower_bound] <- Inf
    return(list(lower = lower, upper = upper))
  }

  inside_support <- numerator > 0
  lower[inside_support] <- 0
  upper[inside_support] <- Inf
  lower[!inside_support] <- Inf
  upper[!inside_support] <- -Inf
  list(lower = lower, upper = upper)
}

.sigma_interval_prob_inv_gamma <- function(lower, upper, shape, rate) {
  prob <- numeric(length(lower))

  full_region <- lower <= 0 & is.infinite(upper)
  prob[full_region] <- 1

  upper_only <- lower <= 0 & is.finite(upper) & upper > 0
  if (any(upper_only)) {
    prob[upper_only] <- stats::pgamma(
      1 / upper[upper_only]^2,
      shape = shape,
      rate = rate[upper_only],
      lower.tail = FALSE
    )
  }

  lower_only <- is.finite(lower) & lower > 0 & is.infinite(upper)
  if (any(lower_only)) {
    prob[lower_only] <- stats::pgamma(
      1 / lower[lower_only]^2,
      shape = shape,
      rate = rate[lower_only],
      lower.tail = TRUE
    )
  }

  finite_band <- is.finite(lower) & lower > 0 &
    is.finite(upper) & upper > lower
  if (any(finite_band)) {
    prob[finite_band] <- stats::pgamma(
      1 / lower[finite_band]^2,
      shape = shape,
      rate = rate[finite_band],
      lower.tail = TRUE
    ) - stats::pgamma(
      1 / upper[finite_band]^2,
      shape = shape,
      rate = rate[finite_band],
      lower.tail = TRUE
    )
  }

  prob[!is.finite(prob) | prob < 0] <- 0
  pmin(prob, 1)
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
get_metric_constraints <- function(metric, c, LSL, USL, target,
                                   sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  list(
    s_max_fn = function() {
      if (metric %in% c("Cp", "Cpm", "Cpc"))
        return(if (c <= 0) Inf else tol / ((2 * sigma_level) * c))
      if (metric == "Cpk")
        return(if (c <= 0) Inf else tol / ((2 * sigma_level) * c))
      if (metric %in% c("Cpu", "Cpl"))
        return(Inf)
      return(Inf)
    },
    # Vectorized version: takes vector of sigma, returns list(lower, upper)
    mu_b_fn_vec = function(s) {
      n <- length(s)
      if (metric == "Cp") {
        return(list(lower = rep(-Inf, n), upper = rep(Inf, n)))
      }
      if (metric == "Cpk") {
        return(list(lower = LSL + sigma_level * c * s,
                    upper = USL - sigma_level * c * s))
      }
      if (metric == "Cpu") {
        return(list(lower = rep(-Inf, n),
                    upper = USL - sigma_level * c * s))
      }
      if (metric == "Cpl") {
        return(list(lower = LSL + sigma_level * c * s,
                    upper = rep(Inf, n)))
      }
      if (metric == "Cpm") {
        if (c <= 0) {
          return(list(lower = rep(-Inf, n), upper = rep(Inf, n)))
        }
        R <- tol / ((2 * sigma_level) * c)
        w <- sqrt(pmax(0, R^2 - s^2))
        lower <- target - w
        upper <- target + w
        invalid <- s >= R
        lower[invalid] <- 0
        upper[invalid] <- -1
        return(list(lower = lower, upper = upper))
      }
      if (metric == "Cpc") {
        if (c <= 0) {
          return(list(lower = rep(-Inf, n), upper = rep(Inf, n)))
        }
        w <- .cpc_mu_width_from_sigma(s, c, tol, sigma_level = sigma_level)
        lower <- target - w
        upper <- target + w
        invalid <- !is.finite(w)
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
compute_metric_value <- function(mu, sigma, LSL, USL, target, metric,
                                 sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  switch(metric,
    "Cp"  = tol / ((2 * sigma_level) * sigma),
    "Cpu" = (USL - mu) / (sigma_level * sigma),
    "Cpl" = (mu - LSL) / (sigma_level * sigma),
    "Cpk" = pmin((USL - mu) / (sigma_level * sigma),
                 (mu - LSL) / (sigma_level * sigma)),
    "Cpm" = tol / ((2 * sigma_level) * sqrt(sigma^2 + (mu - target)^2)),
    "Cpc" = tol / ((2 * sigma_level) * sqrt(pi / 2) *
                     .cpc_E_abs_dev_normal(mu, sigma, target)),
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
                                   cached_state = NULL,
                                   sigma_level = 3) {
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
                                                        p_high = 0.999,
                                                        sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  # Get sigma quantiles using BayesTools::quant
  sigma_low <- BayesTools::quant(prior_sigma, p_low)
  sigma_high <- BayesTools::quant(prior_sigma, p_high)

  # For metrics inversely related to sigma (Cp, Cpk, Cpu, Cpl, Cpm, Cpc),

# metric_high corresponds to sigma_low and vice versa
  metric_at_sigma_low <- switch(metric,
    "Cp" = tol / ((2 * sigma_level) * sigma_low),
    "Cpk" = tol / ((2 * sigma_level) * sigma_low),  # Upper bound (assumes mu = mid)
    "Cpu" = (USL - mid) / (sigma_level * sigma_low),
    "Cpl" = (mid - LSL) / (sigma_level * sigma_low),
    "Cpm" = tol / ((2 * sigma_level) * sigma_low),  # Upper bound (assumes mu = target)
    "Cpc" = tol / ((2 * sigma_level) * sigma_low),  # Upper bound approximation
    tol / ((2 * sigma_level) * sigma_low)  # Default
  )

  metric_at_sigma_high <- switch(metric,
    "Cp" = tol / ((2 * sigma_level) * sigma_high),
    "Cpk" = tol / ((2 * sigma_level) * sigma_high),
    "Cpu" = (USL - mid) / (sigma_level * sigma_high),
    "Cpl" = (mid - LSL) / (sigma_level * sigma_high),
    "Cpm" = tol / ((2 * sigma_level) * sigma_high),
    "Cpc" = tol / ((2 * sigma_level) * sigma_high),
    tol / ((2 * sigma_level) * sigma_high)
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
                                                   cached_state = NULL,
                                                   sigma_level = 3) {
  ss  <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  if (.is_degenerate_conjugate_posterior(beta_n)) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_moments(dist))
    }
  }

  if (!use_analytic) {
    # Numerical fallback: 2D integration
    return(.compute_moments_numerical_conjugate(
      mu_n, k_n, alpha_n, beta_n, LSL, USL, target, metric,
      sigma_level = sigma_level
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
    # Cp = tol / (2 * sigma_level * sigma)
    E1 <- (tol / (2 * sigma_level)) * E_inv_sigma
    E2 <- (tol / (2 * sigma_level))^2 * E_inv_sigma2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpu") {
    # Cpu = (USL - mu) / (sigma_level * sigma)
    # E[Cpu] = E[(USL - mu)/sigma] / sigma_level
    # mu|sigma ~ N(mu_n, sigma^2/k_n), so E[mu|sigma] = mu_n
    # E[(USL - mu)/sigma] = (USL - mu_n) * E[1/sigma]
    E1 <- (USL - mu_n) / sigma_level * E_inv_sigma

    # E[Cpu^2] = E[(USL - mu)^2 / sigma^2] / sigma_level^2
    # (USL - mu)^2 = (USL - mu_n)^2 - 2*(USL - mu_n)*(mu - mu_n) + (mu - mu_n)^2
    # E[(mu - mu_n)^2 | sigma] = sigma^2/k_n
    # E[(USL - mu)^2 / sigma^2] = (USL - mu_n)^2 * E[1/sigma^2] + E[1/k_n] = (USL - mu_n)^2 * E[1/sigma^2] + 1/k_n
    E2 <- ((USL - mu_n)^2 * E_inv_sigma2 + 1 / k_n) / sigma_level^2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpl") {
    # Cpl = (mu - LSL) / (sigma_level * sigma)
    E1 <- (mu_n - LSL) / sigma_level * E_inv_sigma
    E2 <- ((mu_n - LSL)^2 * E_inv_sigma2 + 1 / k_n) / sigma_level^2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpk") {
    # Cpk = min(Cpu, Cpl) = (tol/2 - |mu - mid|) / (sigma_level * sigma)
    # E[Cpk] = (tol/2) / sigma_level * E[1/sigma] -
    #          E[|mu - mid|/sigma] / sigma_level

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

    E1 <- (tol / 2) / sigma_level * E_inv_sigma - E_abs_div_sigma / sigma_level

    # E[Cpk^2] - more complex, use 1D numerical integration
    integrand_sq <- function(y) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      # E[(tol/2 - |mu - mid|)^2 | sigma] / (sigma_level^2 * sigma^2)
      # = E[(tol/2)^2 - tol*|mu-mid| + |mu-mid|^2 | sigma] / (9*sigma^2)
      abs_mean <- tau * sqrt(2 / pi) * exp(-delta^2 / (2 * tau^2)) +
                  delta * (1 - 2 * stats::pnorm(-delta / tau))
      # E[|X|^2] = E[X^2] = delta^2 + tau^2
      abs2_mean <- delta^2 + tau^2
      cpk2_given_sigma <- ((tol / 2)^2 - tol * abs_mean + abs2_mean) /
        (sigma_level^2 * sigma^2)
      cpk2_given_sigma * stats::dchisq(y, df_p)
    }
    E2 <- stats::integrate(integrand_sq, 0, Inf, rel.tol = 1e-6)$value

    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpm") {
    # Cpm = tol / (2 * sigma_level * sqrt(sigma^2 + (mu - T)^2))
    # No closed form - use 1D numerical integration over sigma (via chi-square)
    delta_T <- mu_n - target

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
    E1 <- (tol / (2 * sigma_level)) * E_inv_sqrt

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
    E2 <- (tol / (2 * sigma_level))^2 * E_inv

    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpc") {
    y_upper <- max(100, df_p * 10)

    gh_nodes <- sqrt(2) * c(
      -3.190994, -2.266581, -1.468554, -0.723551, 0,
       0.723551,  1.468554,  2.266581,  3.190994
    )
    gh_wts <- c(
      3.961e-05, 4.944e-03, 8.847e-02, 4.326e-01, 7.202e-01,
      4.326e-01, 8.847e-02, 4.944e-03, 3.961e-05
    )
    gh_wts <- gh_wts / sum(gh_wts)
    n_gh <- length(gh_nodes)

    integrand_cpc <- function(y, power) {
      sigma <- sqrt(2 * beta_n / y)
      tau <- sigma / sqrt(k_n)
      mu_mat <- tcrossprod(tau, gh_nodes) + mu_n
      sigma_rep <- rep(sigma, times = n_gh)
      metric_vals <- compute_metric_value(
        as.vector(mu_mat), sigma_rep, LSL, USL, target, metric,
        sigma_level = sigma_level
      )
      metric_mat <- matrix(metric_vals, nrow = length(sigma), ncol = n_gh)
      inner <- drop((metric_mat^power) %*% gh_wts)
      inner * stats::dchisq(y, df_p)
    }

    E1 <- stats::integrate(
      function(y) integrand_cpc(y, 1),
      1e-6, y_upper,
      rel.tol = 1e-4,
      subdivisions = 200
    )$value
    E2 <- stats::integrate(
      function(y) integrand_cpc(y, 2),
      1e-6, y_upper,
      rel.tol = 1e-4,
      subdivisions = 200
    )$value

    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }


  stop("Unknown metric: ", metric)
}

#' Numerical fallback for conjugate prior moments (2D integration)
#' @keywords internal
.compute_moments_numerical_conjugate <- function(mu_n, k_n, alpha_n, beta_n,
                                                  LSL, USL, target, metric,
                                                  sigma_level = 3) {
  if (.is_degenerate_conjugate_posterior(beta_n)) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_moments(dist))
    }
  }

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
    metric_val <- compute_metric_value(
      mu, sigma, LSL, USL, target, metric,
      sigma_level = sigma_level
    )

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
                                                         cached_state = NULL,
                                                         sigma_level = 3) {
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
      inv_sigma_scale <- 1 / (sigma_level * sigma)
      switch(metric,
        "Cp"  = list(E1 = spec_tol * inv_sigma_scale / 2,
                     E2 = (spec_tol * inv_sigma_scale / 2)^2),
        "Cpu" = list(E1 = cc * inv_sigma_scale,
                     E2 = (cc^2 + b^2) * inv_sigma_scale^2),
        "Cpl" = list(E1 = a * inv_sigma_scale,
                     E2 = (a^2 + b^2) * inv_sigma_scale^2),
        "Cpk" = {
          zs  <- (M - mu_n) / sd_mu
          Phi <- stats::pnorm(zs); phi <- stats::dnorm(zs)
          list(
            E1 = (a * Phi + cc * (1 - Phi) - 2 * b * phi) * inv_sigma_scale,
            E2 = ((a^2 + b^2) * Phi + (cc^2 + b^2) * (1 - Phi) -
                    2 * b * (a + cc) * phi) * inv_sigma_scale^2
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
      m_vals <- compute_metric_value(
        as.vector(t(mu_mat)), sigma_rep, LSL, USL, target, metric,
        sigma_level = sigma_level
      )
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
                                                            cached_state = NULL,
                                                            sigma_level = 3) {
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
        rep(mu_vec[i], n_gh), sigma_pts, LSL, USL, target, metric,
        sigma_level = sigma_level
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
                                                 cached_state = NULL,
                                                 sigma_level = 3) {
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
    metric_val <- compute_metric_value(
      mu, sigma, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
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
make_solver <- function(data, LSL, USL, prior, metric = "Cpk", target = NULL,
                        sigma_level = 3, ...) {
  UseMethod("make_solver", prior)
}

#' @export
make_solver.PriorConjugate <- function(data, LSL, USL, prior, metric = "Cpk",
                                        target = NULL, cached_state = NULL,
                                        sigma_level = 3, ...) {
  ss   <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n
  df_p <- 2 * alpha_n
  tol <- USL - LSL

  if (.is_degenerate_conjugate_posterior(beta_n)) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_solver(dist))
    }
  }

  # Analytic fast paths for Cp, Cpu, Cpl
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(1.0)
      stats::pchisq(8 * beta_n * sigma_level^2 * c^2 / tol^2,
                    df = df_p, lower.tail = FALSE)
    })
  }
  if (metric == "Cpu") {
    x_U <- (USL - mu_n) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      suppressWarnings(stats::pt(x_U, df = df_p,
                                 ncp = sigma_level * c * sqrt_k_n))
    })
  }
  if (metric == "Cpl") {
    x_L <- (mu_n - LSL) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      suppressWarnings(stats::pt(x_L, df = df_p,
                                 ncp = sigma_level * c * sqrt_k_n))
    })
  }

  # Pre-compute global h_max (chi-square mode density for numerical stability)
  y_mode <- max(df_p - 2, 1e-6)
  h_max_global <- stats::dchisq(y_mode, df_p, log = TRUE)

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)

    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
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
                                              cached_state = NULL,
                                              sigma_level = 3, ...) {
  # Case 2: Conjugate mu (Normal/Jeffreys), non-conjugate sigma
  # Integrate mu analytically using pnorm, then 1D numerical over sigma

  ss <- .extract_suff_stats(data, cached_state)
  n <- ss$n; x_bar <- ss$x_bar; sse <- ss$SS
  smp <- .semi_mu_posterior(prior, n, x_bar, sse)
  k_n <- smp$k_n; mu_n <- smp$mu_n; sse_n <- smp$sse_n

  n_eff <- if (prior$k0 == 0) n - 1 else n
  # make_solver() must return a survival probability, not raw posterior mass.
  log_Z <- .compute_semi_mu_log_Z(n_eff, sse_n, prior$log_dens_sigma)

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
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

      exp(log_lik + log_prior_sigma + log_prob_mu - log_Z)
    }

    sigma_upper <- if (is.infinite(s_max)) Inf else max(s_max, 1e-6)

    prob <- stats::integrate(integrand, 1e-10, sigma_upper,
                             rel.tol = 1e-5, subdivisions = 200)$value
    min(max(prob, 0), 1)
  }
}

#' @export
make_solver.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 cached_state = NULL,
                                                 sigma_level = 3, ...) {
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
  alpha_0_times_logbeta0 <- if (beta0 == 0) 0 else alpha0 * log(beta0)
  jeffreys_adj <- if (alpha0 == -0.5 && beta0 == 0) 0.5 else 0
  # Normalize to the same posterior mass used by the density solver.
  log_Z <- .compute_semi_sigma_log_Z(n, x_bar, sse, alpha0, beta0,
                                     prior$log_dens_mu, jeffreys_adj)
  mu_window <- 20 * sqrt(sse / max(n, 1))

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)

    integrand <- function(mu) {
      sse_mu <- sse + n * (mu - x_bar)^2
      alpha_n <- alpha0 + n / 2
      alpha_n_eff <- alpha_n + jeffreys_adj
      beta_n <- beta0 + sse_mu / 2
      log_marginal <- lgamma(alpha_n_eff) - lgamma(alpha0) +
        alpha_0_times_logbeta0 - alpha_n_eff * log(beta_n)
      log_prior_mu <- prior$log_dens_mu(mu)

      sigma_region <- .metric_sigma_region_from_mu(
        metric, mu, c, LSL, USL, target,
        sigma_level = sigma_level
      )
      prob_sigma <- .sigma_interval_prob_inv_gamma(
        sigma_region$lower, sigma_region$upper,
        shape = alpha_n_eff, rate = beta_n
      )
      log_prob_sigma <- rep(-Inf, length(mu))
      pos_prob <- prob_sigma > 0
      log_prob_sigma[pos_prob] <- log(prob_sigma[pos_prob])

      exp(log_marginal + log_prior_mu + log_prob_sigma - log_Z)
    }

    mu_low <- x_bar - mu_window
    mu_high <- x_bar + mu_window

    prob <- stats::integrate(integrand, mu_low, mu_high,
                             rel.tol = 1e-5, subdivisions = 200)$value
    min(max(prob, 0), 1)
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
make_density_solver <- function(data, LSL, USL, prior, metric = "Cpk",
                                target = NULL, sigma_level = 3, ...) {
  UseMethod("make_density_solver", prior)
}

#' @export
#' @export
make_density_solver.PriorConjugate <- function(data, LSL, USL, prior,
                                                metric = "Cpk", target = NULL,
                                                cached_state = NULL,
                                                sigma_level = 3, ...) {
  ss   <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  M <- (LSL + USL) / 2
  tol <- USL - LSL

  if (.is_degenerate_conjugate_posterior(beta_n)) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_density(dist))
    }
  }

  log_Z_sigma <- lgamma(alpha_n) - log(2) - alpha_n * log(beta_n)

  # Cp: no integration needed (contour is a single sigma point)
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)
      log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_c) - beta_n / sigma_c^2
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)
      exp(log_p_sigma + log_jacobian - log_Z_sigma)
    })
  }

  # Cpu/Cpl: analytic via non-central t finite-differencing
  if (metric == "Cpu") {
    x_U <- (USL - mu_n) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      if (c <= 0) return(0)
      h <- max(c * 1e-5, 1e-8)
      S_plus  <- suppressWarnings(stats::pt(x_U, df = 2 * alpha_n,
                                            ncp = sigma_level * (c + h) * sqrt_k_n))
      S_minus <- suppressWarnings(stats::pt(x_U, df = 2 * alpha_n,
                                            ncp = sigma_level * (c - h) * sqrt_k_n))
      max(0, -(S_plus - S_minus) / (2 * h))
    })
  }

  if (metric == "Cpl") {
    x_L <- (mu_n - LSL) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      if (c <= 0) return(0)
      h <- max(c * 1e-5, 1e-8)
      S_plus  <- suppressWarnings(stats::pt(x_L, df = 2 * alpha_n,
                                            ncp = sigma_level * (c + h) * sqrt_k_n))
      S_minus <- suppressWarnings(stats::pt(x_L, df = 2 * alpha_n,
                                            ncp = sigma_level * (c - h) * sqrt_k_n))
      max(0, -(S_plus - S_minus) / (2 * h))
    })
  }

  # Precompute sigma grid for metrics requiring contour integration.
  # Replaces per-call stats::integrate with a fixed trapezoidal rule.
  sigma_mode_approx <- sqrt(beta_n / max(alpha_n, 1))
  sigma_lo <- max(1e-10, sigma_mode_approx * 0.02)
  sigma_hi <- sigma_mode_approx * 10
  n_sigma <- 1024L
  sigma_grid <- seq(sigma_lo, sigma_hi, length.out = n_sigma)
  lps_grid <- -(2 * alpha_n + 1) * log(sigma_grid) - beta_n / sigma_grid^2
  log_jac_grid <- log(sigma_level) + log(sigma_grid)

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / ((2 * sigma_level) * c)
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
      mu_L <- LSL + sigma_level * c * sig_pts
      mu_U <- USL - sigma_level * c * sig_pts
      lps <- c(lps_grid[idx], -(2 * alpha_n + 1) * log(sig_pts[np]) - beta_n / sig_pts[np]^2)
      lj  <- c(log_jac_grid[idx], log(sigma_level) + log(sig_pts[np]))
      log_p_mu_L <- stats::dnorm(mu_L, mu_n, sd_mu, log = TRUE)
      log_p_mu_U <- stats::dnorm(mu_U, mu_n, sd_mu, log = TRUE)
      contrib_L <- exp(lps + log_p_mu_L + lj - log_Z_sigma)
      contrib_U <- exp(lps + log_p_mu_U + lj - log_Z_sigma)
      contrib_L[!is.finite(contrib_L)] <- 0
      contrib_U[!is.finite(contrib_U)] <- 0
      sum((contrib_L + contrib_U) * w_trap)
    })
  }

  if (metric == "Cpm") {
    if (is.null(target)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)
      R <- tol / ((2 * sigma_level) * c)
      if (R <= 0) return(0)
      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      valid <- sigma_v > 1e-10
      if (!any(valid)) return(0)
      sv <- sigma_v[valid]
      mu_v <- target + R * cos(theta_grid[valid])
      sd_mu <- sv / sqrt(k_n)

      log_p_sigma <- -(2 * alpha_n + 1) * log(sv) - beta_n / sv^2
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)
      vals <- exp(log_p_sigma + log_p_mu + log_jacobian - log_Z_sigma)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M

    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )
        sd_mu <- contour$sigma / sqrt(k_n)
        log_p_sigma <- -(2 * alpha_n + 1) * log(contour$sigma) -
          beta_n / contour$sigma^2
        log_p_mu_L <- stats::dnorm(contour$mu_lower, mu_n, sd_mu, log = TRUE)
        log_p_mu_U <- stats::dnorm(contour$mu_upper, mu_n, sd_mu, log = TRUE)

        vals <- exp(log_p_sigma + log_p_mu_L + contour$log_jacobian - log_Z_sigma) +
          exp(log_p_sigma + log_p_mu_U + contour$log_jacobian - log_Z_sigma)
        vals[!is.finite(vals)] <- 0
        vals
      }

      stats::integrate(
        integrand, 0, 20,
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' @export
make_density_solver.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                                      metric = "Cpk", target = NULL,
                                                      cached_state = NULL,
                                                      sigma_level = 3, ...) {
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

  sigma_support <- list(lower = 0, upper = Inf)
  bt_priors <- prior$bayestools_priors
  if (!is.null(bt_priors) && inherits(bt_priors$sigma, "prior")) {
    sigma_support <- .extract_prior_bounds(bt_priors$sigma)
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
  log_jac_grid <- log(sigma_level) + log(sigma_grid)

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / ((2 * sigma_level) * c)
      if (sigma_max <= 0) return(0)

      valid <- sigma_grid > 0 & sigma_grid < sigma_max
      if (!any(valid)) return(0)
      sg <- sigma_grid[valid]
      sd_mu <- sg / sqrt(k_n)

      mu_L <- LSL + sigma_level * c * sg
      mu_U <- USL - sigma_level * c * sg

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

      mu_U <- USL - sigma_level * c * sigma_grid
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

      mu_L <- LSL + sigma_level * c * sigma_grid
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

      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0 ||
          sigma_c < sigma_support$lower ||
          sigma_c > sigma_support$upper) {
        return(0)
      }

      lps <- log_p_sigma(sigma_c)
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)
      exp(lps + log_jacobian - log_Z)
    })
  }

  if (metric == "Cpm") {
    if (is.null(target)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)

      R <- tol / ((2 * sigma_level) * c)
      if (R <= 0) return(0)

      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      mu_v <- target + R * cos(theta_grid)
      sd_mu <- sigma_v / sqrt(k_n)

      lps <- log_p_sigma(sigma_v)
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)

      vals <- exp(lps + log_p_mu + log_jacobian - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M

    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )
        sd_mu <- contour$sigma / sqrt(k_n)
        lps <- log_p_sigma(contour$sigma)
        log_p_mu_L <- stats::dnorm(contour$mu_lower, mu_n, sd_mu, log = TRUE)
        log_p_mu_U <- stats::dnorm(contour$mu_upper, mu_n, sd_mu, log = TRUE)

        vals <- exp(lps + log_p_mu_L + contour$log_jacobian - log_Z) +
          exp(lps + log_p_mu_U + contour$log_jacobian - log_Z)
        vals[!is.finite(vals)] <- 0
        vals
      }

      stats::integrate(
        integrand, 0, 20,
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
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
                                                         cached_state = NULL,
                                                         sigma_level = 3, ...) {
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
  log_jac_grid <- log(sigma_level) + log_sigma_grid

  if (metric == "Cpk") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_max <- tol / ((2 * sigma_level) * c)
      if (sigma_max <= 0) return(0)

      valid <- sigma_grid > 0 & sigma_grid < sigma_max
      if (!any(valid)) return(0)
      sg <- sigma_grid[valid]

      mu_L <- LSL + sigma_level * c * sg
      mu_U <- USL - sigma_level * c * sg

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

      mu_U <- USL - sigma_level * c * sigma_grid
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

      mu_L <- LSL + sigma_level * c * sigma_grid
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
      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)

      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)

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

  if (metric == "Cpm") {
    if (is.null(target)) stop("Target required for Cpm")

    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)

      R <- tol / ((2 * sigma_level) * c)
      if (R <= 0) return(0)

      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      mu_v <- target + R * cos(theta_grid)

      sse_mu <- sse + n * (mu_v - x_bar)^2
      beta_n <- beta0 + sse_mu / 2
      log_sk <- log_C - (2 * alpha_n_eff + 1) * log(sigma_v) - beta_n / sigma_v^2
      log_pm <- prior$log_dens_mu(mu_v)

      vals <- exp(log_pm + log_sk + log_jacobian - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M

    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )

        sse_L <- sse + n * (contour$mu_lower - x_bar)^2
        sse_U <- sse + n * (contour$mu_upper - x_bar)^2
        beta_n_L <- beta0 + sse_L / 2
        beta_n_U <- beta0 + sse_U / 2

        log_sk_L <- log_C - (2 * alpha_n_eff + 1) * log(contour$sigma) -
          beta_n_L / contour$sigma^2
        log_sk_U <- log_C - (2 * alpha_n_eff + 1) * log(contour$sigma) -
          beta_n_U / contour$sigma^2
        log_pm_L <- prior$log_dens_mu(contour$mu_lower)
        log_pm_U <- prior$log_dens_mu(contour$mu_upper)

        vals <- exp(log_pm_L + log_sk_L + contour$log_jacobian - log_Z) +
          exp(log_pm_U + log_sk_U + contour$log_jacobian - log_Z)
        vals[!is.finite(vals)] <- 0
        vals
      }

      stats::integrate(
        integrand, 0, 20,
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
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
                                              cached_state = NULL,
                                              sigma_level = 3, ...) {
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

      sigma_max <- tol / ((2 * sigma_level) * c)
      if (sigma_max <= 0) return(0)
      sigma_upper <- min(sigma_max, uni_s)
      if (sigma_upper <= 1e-10) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0 & sigma < sigma_upper
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]

        # Contour points
        mu_L <- LSL + sigma_level * c * sigma_v
        mu_U <- USL - sigma_level * c * sigma_v

        # Jacobian |dmu/dc| = sigma_level * sigma
        log_jacobian <- log(sigma_level) + log(sigma_v)

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

      tryCatch({
        stats::integrate(integrand, 1e-10, sigma_upper, rel.tol = 1e-5, subdivisions = 200)$value
      }, error = function(e) {
        0
      })
    })
  }

  # Cpu: single contour
  # For Cpu = c, contour is mu = USL - sigma_level * c * sigma
  # Need sigma_upper that ensures mu is within prior support (or has meaningful density)
  if (metric == "Cpu") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        mu_U <- USL - sigma_level * c * sigma_v
        log_jacobian <- log(sigma_level) + log(sigma_v)

        log_p <- log_post_vec(mu_U, sigma_v)
        contrib <- exp(log_p + log_jacobian - h_max) / Z
        contrib[!is.finite(contrib)] <- 0

        result[valid] <- contrib
        result
      }

      mu_center <- if (n > 0) x_bar else map_mu
      mu_sd <- if (n > 0) sqrt(uni_s^2 / n) else if (!is.null(mu_scale)) mu_scale else uni_s
      mu_lower_bound <- mu_center - 6 * mu_sd
      sigma_upper_mu <- (USL - mu_lower_bound) / (sigma_level * c)
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
  # For Cpl = c, contour is mu = LSL + sigma_level * c * sigma
  if (metric == "Cpl") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(sigma) {
        valid <- sigma > 0
        result <- rep(0, length(sigma))
        if (!any(valid)) return(result)

        sigma_v <- sigma[valid]
        mu_L <- LSL + sigma_level * c * sigma_v
        log_jacobian <- log(sigma_level) + log(sigma_v)

        log_p <- log_post_vec(mu_L, sigma_v)
        contrib <- exp(log_p + log_jacobian - h_max) / Z
        contrib[!is.finite(contrib)] <- 0

        result[valid] <- contrib
        result
      }

      mu_center <- if (n > 0) x_bar else map_mu
      mu_sd <- if (n > 0) sqrt(uni_s^2 / n) else if (!is.null(mu_scale)) mu_scale else uni_s
      mu_upper_bound <- mu_center + 6 * mu_sd
      sigma_upper_mu <- (mu_upper_bound - LSL) / (sigma_level * c)
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

      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)

      # At sigma_c, mu can be anything - integrate over mu
      # p(Cp = c) = ∫ p(mu, sigma_c) |dsigma/dc| dmu
      # |dsigma/dc| = tol / (2 * sigma_level * c^2)
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)

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

      tryCatch({
        mu_integral <- stats::integrate(mu_integrand, mu_lower, mu_upper, rel.tol = 1e-5)$value
        mu_integral * exp(log_jacobian)
      }, error = function(e) {
        0
      })
    })
  }

  if (metric == "Cpm") {
    if (is.null(target)) stop("Target required for Cpm")

    return(function(c) {
      if (c <= 0) return(0)

      # Contour radius R = tol / (2 * sigma_level * c)
      R <- tol / ((2 * sigma_level) * c)
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
        mu_v <- target + R * cos(theta[valid])

        # Log joint posterior
        log_p <- log_post_vec(mu_v, sigma_v)

        # Combined density (normalized by Z)
        log_contrib <- log_p + log_jacobian - h_max

        result[valid] <- exp(log_contrib) / Z
        result
      }

      tryCatch({
        stats::integrate(integrand, 0, pi, rel.tol = 1e-5, subdivisions = 200)$value
      }, error = function(e) {
        0
      })
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M

    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )
        log_p_L <- log_post_vec(contour$mu_lower, contour$sigma)
        log_p_U <- log_post_vec(contour$mu_upper, contour$sigma)

        vals <- exp(log_p_L + contour$log_jacobian - h_max) / Z +
          exp(log_p_U + contour$log_jacobian - h_max) / Z
        vals[!is.finite(vals)] <- 0
        vals
      }

      tryCatch({
        stats::integrate(
          integrand, 0, 20,
          rel.tol = 1e-5,
          subdivisions = 400
        )$value
      }, error = function(e) {
        0
      })
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}


#' @export
make_solver.PriorGeneric <- function(data, LSL, USL, prior, metric = "Cpk",
                                      target = NULL, cached_state = NULL,
                                      sigma_level = 3, ...) {

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
    init <- .find_feasible_generic_init(x_bar, init_sd, prior, log_post)
    opt <- stats::optim(c(init$mu, init$sigma), function(p) -log_post(p[1], p[2]))
    map_mu <- opt$par[1]
    map_sig <- opt$par[2]
    h_max <- -opt$value

    # Tighter bounds intersected with prior support when available.
    uni_s <- map_sig * 5
    bt <- prior$bayestools_priors
    if (!is.null(bt) && inherits(bt$sigma, "prior")) {
      sig_bounds <- .extract_prior_bounds(bt$sigma)
      if (is.finite(sig_bounds$upper)) {
        uni_s <- min(uni_s, sig_bounds$upper)
      }
      uni_s <- max(uni_s, max(sig_bounds$lower, 1e-6) * 1.1)
    }

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
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
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
precompute_generic_state <- function(data, prior, cached_state = NULL) {
  if (!is.null(cached_state)) {
    n <- cached_state$n %||% 0L
    x_bar <- cached_state$x_bar %||% 0
    sse <- cached_state$sse %||% 0
  } else {
    n <- length(data)
    if (n > 0) {
      x_bar <- mean(data)
      sse <- sum((data - x_bar)^2)
    } else {
      x_bar <- 0
      sse <- 0
    }
  }

  if (n > 0) {
    init_sd <- max(sqrt(sse / max(n - 1, 1)), 1e-6)
    init_mu <- x_bar
  } else {
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

  init <- .find_feasible_generic_init(init_mu, init_sd, prior, log_post)
  init_mu <- init$mu
  init_sd <- init$sigma

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

    bt <- prior$bayestools_priors
    if (!is.null(bt) && inherits(bt$sigma, "prior")) {
      sig_bounds <- .extract_prior_bounds(bt$sigma)
      if (is.finite(sig_bounds$lower)) {
        sigma_lower <- max(sigma_lower, sig_bounds$lower + 1e-6)
      }
      if (is.finite(sig_bounds$upper)) {
        sigma_upper <- min(sigma_upper, sig_bounds$upper)
      }
    }
    sigma_upper <- max(sigma_upper, sigma_lower * 1.1)

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
                                          cached_state = NULL,
                                          sigma_level = 3) {
  lower_bound <- min(bounds)
  upper_bound <- max(bounds)

  interval_prob_from_survival <- function(S) {
    p_lower <- if (is.finite(lower_bound)) S(lower_bound) else 1
    p_upper <- if (is.finite(upper_bound)) S(upper_bound) else 0
    min(max(p_lower - p_upper, 0), 1)
  }

  if (inherits(prior, "PriorConjugate")) {
    ss <- .extract_suff_stats(data, cached_state)
    post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
    if (.is_degenerate_conjugate_posterior(post$beta_n)) {
      dist <- .degenerate_conjugate_metric_distribution(
        post$mu_n, post$k_n, LSL, USL, target, metric,
        sigma_level = sigma_level
      )
      if (!is.null(dist)) {
        return(.degenerate_metric_prob(dist, bounds))
      }
    }
  }

  # Check if we can use density solver
  density_metrics <- c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")
  can_use_density <- (inherits(prior, "PriorConjugate") ||
                      inherits(prior, "PriorGeneric") ||
                      inherits(prior, "PriorSemiConjugateMu") ||
                      inherits(prior, "PriorSemiConjugateSigma")) &&
                     metric %in% density_metrics &&
                     !.metric_can_be_negative(metric)

  if (can_use_density) {
    # Use density solver PDF and integrate over bounds
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                  sigma_level = sigma_level,
                                  cached_state = cached_state)
    pdf_vec <- function(x) vapply(x, pdf_fn, numeric(1L))

    # Capability indices on this path have support on [0, Inf).
    lower <- max(lower_bound, 0)
    upper <- upper_bound
    if (upper <= lower) return(0)

    prob <- tryCatch(
      stats::integrate(pdf_vec, lower, upper, rel.tol = 1e-4)$value,
      error = function(e) NA_real_
    )
    if (!is.na(prob)) return(prob)

    # Fallback: try survival function approach
    S <- tryCatch({
      if (!is.null(cached_state))
        make_solver(data, LSL, USL, prior, metric, target,
                    sigma_level = sigma_level,
                    cached_state = cached_state)
      else
        make_solver(data, LSL, USL, prior, metric, target,
                    sigma_level = sigma_level)
    }, error = function(e) NULL)
    if (!is.null(S)) {
      prob <- tryCatch(interval_prob_from_survival(S), error = function(e) NA_real_)
      if (!is.na(prob)) return(prob)
    }
    return(NA_real_)

  } else {
    # Fallback: Survival function S(c) = P(Index > c)
    if (!is.null(cached_state)) {
      S <- make_solver(data, LSL, USL, prior, metric, target,
                       sigma_level = sigma_level,
                       cached_state = cached_state)
    } else {
      S <- make_solver(data, LSL, USL, prior, metric, target,
                       sigma_level = sigma_level)
    }

    # P(lower < Index < upper) = P(Index > lower) - P(Index > upper)
    return(interval_prob_from_survival(S))
  }
}

#' Compute negative-value mean correction for prior-only mode.
#' For Cpu, Cpl, Cpk, the metric can be negative when mu is outside \eqn{[LSL, USL]}.
#' This returns \eqn{E[metric * I(metric <= 0)]}, a negative correction to add to
#' \eqn{area * E[metric | metric > 0]} to get the full unconditional mean.
#' @keywords internal
.compute_negative_mean_correction <- function(metric, LSL, USL, bayestools_priors,
                                              sigma_level = 3) {
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
    corr_Cpu <- (1 / sigma_level) * E_Cpu_neg * E_inv_sigma
  }

  if (metric == "Cpl" || metric == "Cpk") {
    # E[(mu - LSL) * I(mu < LSL)]: negative since mu - LSL < 0 when mu < LSL
    E_Cpl_neg <- tryCatch({
      integrand <- function(mu) (mu - LSL) * exp(log_dens_mu(mu))
      stats::integrate(Vectorize(integrand), LSL - 20 * tol, LSL,
                       rel.tol = 1e-6)$value
    }, error = function(e) 0)
    corr_Cpl <- (1 / sigma_level) * E_Cpl_neg * E_inv_sigma
  }

  switch(metric,
    "Cpu" = corr_Cpu,
    "Cpl" = corr_Cpl,
    "Cpk" = corr_Cpu + corr_Cpl,
    0
  )
}

# Two-pass adaptive grid for density evaluation.
# Pass 1 (coarse): ~n_coarse uniform points to locate the significant region
#   and extend edges if density is non-negligible at boundaries.
# Pass 2 (fine):   concentrate remaining points in the significant region
#   (where pdf > 1% of peak), with sparse coverage in the tails.
# Returns a list with grid_x, mid_x, pdf_vals, dx, area — ready for
# downstream statistics and plotting.
.adaptive_density_grid <- function(pdf_fn, x_start, x_end, n_grid,
                                   is_prior_only = FALSE) {

  n_coarse <- min(32L, n_grid)

  # --- Pass 1: coarse grid + edge extension ---
  max_extend <- if (is_prior_only) 5L else 3L
  for (attempt in seq_len(max_extend)) {
    coarse_x   <- seq(x_start, x_end, length.out = n_coarse)
    coarse_mid <- (coarse_x[-1] + coarse_x[-n_coarse]) / 2
    coarse_pdf <- vapply(coarse_mid, pdf_fn, numeric(1))
    coarse_pdf[!is.finite(coarse_pdf) | coarse_pdf < 0] <- 0

    peak <- max(coarse_pdf)
    if (peak == 0) break
    edge_threshold <- 0.01 * peak
    need_left  <- coarse_pdf[1] > edge_threshold && x_start > 0
    need_right <- coarse_pdf[length(coarse_pdf)] > edge_threshold

    if (!need_left && !need_right) break
    width <- x_end - x_start
    if (need_left)  x_start <- max(0, x_start - width)
    if (need_right) x_end   <- x_end + width
  }

  if (peak == 0) {
    grid_x   <- seq(x_start, x_end, length.out = n_grid)
    mid_x    <- (grid_x[-1] + grid_x[-n_grid]) / 2
    pdf_vals <- rep(0, length(mid_x))
    dx       <- diff(grid_x)
    return(list(grid_x = grid_x, mid_x = mid_x, pdf_vals = pdf_vals,
                dx = dx, area = 0))
  }

  # For broad prior-only distributions, use log-spaced grid
  if (is_prior_only && x_end > 0 && x_start > 0 && x_end / x_start > 20) {
    n_grid_po <- max(n_grid, 1024L)
    grid_x   <- exp(seq(log(max(x_start, 1e-4)), log(x_end), length.out = n_grid_po))
    mid_x    <- (grid_x[-1] + grid_x[-n_grid_po]) / 2
    pdf_vals <- vapply(mid_x, pdf_fn, numeric(1))
    pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0
    dx   <- diff(grid_x)
    area <- sum(pdf_vals * dx)
    if (area > 0) pdf_vals <- pdf_vals / area
    return(list(grid_x = grid_x, mid_x = mid_x, pdf_vals = pdf_vals,
                dx = dx, area = area))
  }

  if (is_prior_only) n_grid <- max(n_grid, 1024L)

  # --- Pass 2: adaptive refinement ---
  # Identify significant region from coarse grid (pdf > 1% of peak)
  sig <- coarse_pdf > edge_threshold
  sig_lo <- coarse_mid[which.max(sig)]
  sig_hi <- coarse_mid[length(sig) - which.max(rev(sig)) + 1L]

  # Allocate ~75% of budget to the significant region, ~25% to tails
  n_fine <- max(16L, as.integer(0.75 * n_grid))
  n_tail <- n_grid - n_fine

  fine_grid <- seq(sig_lo, sig_hi, length.out = n_fine)

  # Tail points: split between left and right tails
  left_width  <- sig_lo - x_start
  right_width <- x_end - sig_hi
  total_tail  <- left_width + right_width
  if (total_tail > 0 && n_tail >= 2L) {
    n_left  <- max(1L, round(n_tail * left_width / total_tail))
    n_right <- max(1L, n_tail - n_left)
    left_grid  <- if (n_left  >= 2L && left_width  > 0) seq(x_start, sig_lo, length.out = n_left + 1L)[-(n_left + 1L)]  else numeric(0)
    right_grid <- if (n_right >= 2L && right_width > 0) seq(sig_hi, x_end, length.out = n_right + 1L)[-1L] else numeric(0)
  } else {
    left_grid <- numeric(0)
    right_grid <- numeric(0)
  }

  grid_x <- sort(unique(c(x_start, left_grid, fine_grid, right_grid, x_end)))
  mid_x    <- (grid_x[-1] + grid_x[-length(grid_x)]) / 2
  pdf_vals <- vapply(mid_x, pdf_fn, numeric(1))
  pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0
  dx   <- diff(grid_x)
  area <- sum(pdf_vals * dx)
  if (area > 0) pdf_vals <- pdf_vals / area

  list(grid_x = grid_x, mid_x = mid_x, pdf_vals = pdf_vals,
       dx = dx, area = area)
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
                                            n_grid = 512L,
                                            cached_state = NULL,
                                            use_density_solver = TRUE,
                                            mc_samples = NULL,
                                            divergence_info = NULL,
                                            sigma_level = 3) {

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
  metric_can_be_negative <- .metric_can_be_negative(metric)

  if (inherits(prior, "PriorConjugate")) {
    ss <- .extract_suff_stats(data, cached_state)
    post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
    if (.is_degenerate_conjugate_posterior(post$beta_n)) {
      dist <- .degenerate_conjugate_metric_distribution(
        post$mu_n, post$k_n, LSL, USL, target, metric,
        sigma_level = sigma_level
      )
      if (!is.null(dist)) {
        return(.analyze_degenerate_metric_distribution(
          metric, dist, n_grid,
          divergence_info = divergence_info,
          metric_can_be_negative = metric_can_be_negative
        ))
      }
    }
  }

  if (is_prior_only && has_bayestools_priors) {
    # Fast path: compute grid bounds directly from sigma quantiles
    bounds <- .compute_metric_grid_bounds_from_quantiles(
      prior$bayestools_priors$sigma, LSL, USL, metric, target,
      sigma_level = sigma_level
    )
    x_start <- bounds$x_start
    x_end <- bounds$x_end
  } else if (divergence_info$mean_divergent) {
    # Divergent mean: skip moment computation, use solver-based grid bounds
    x_start <- if (metric_can_be_negative) -3 else 0
    x_end <- 3
  } else {
    # Standard path: compute posterior moments for grid bounds
    moments <- compute_metric_moments(
      data, LSL, USL, prior, metric = metric, target = target,
      use_analytic = TRUE, cached_state = cached_state,
      sigma_level = sigma_level
    )

    # Grid bounds via normal approximation (clip to natural bounds)
    # Handle NaN moments by falling back to reasonable defaults
    if (is.na(moments$mean) || is.na(moments$sd) || moments$sd <= 0) {
      # Fallback: use a wide grid centered around 1.0
      x_start <- if (metric_can_be_negative) -3 else 0
      x_end <- 3
    } else {
      x_start <- moments$mean - 3 * moments$sd
      if (!metric_can_be_negative) {
        x_start <- max(0, x_start)
      }
      x_end <- moments$mean + 3 * moments$sd
    }
  }

  # Negative-support metrics need explicit lower-tail coverage because the
  # direct density backend only covers c > 0. We use survival quantiles to
  # capture the full support before building the grid.
  if (metric_can_be_negative || divergence_info$mean_divergent || x_end > 50) {
    try({
      S_fn <- if (!is.null(cached_state)) {
        make_solver(data, LSL, USL, prior, metric, target,
                    sigma_level = sigma_level,
                    cached_state = cached_state)
      } else {
        make_solver(data, LSL, USL, prior, metric, target,
                    sigma_level = sigma_level)
      }

      if (metric_can_be_negative) {
        x_start <- min(x_start, -1)
        for (try_limit in c(-1, -2, -5, -10, -20, -50, -100, -500, -1000)) {
          if (S_fn(try_limit) > 0.999) {
            x_start <- try_limit
            break
          }
        }
      }

      for (try_limit in c(1, 2, 5, 10, 20, 50, 100, 500, 1000)) {
        if (S_fn(try_limit) < 0.001) {
          x_end <- if (divergence_info$mean_divergent || x_end > 50) {
            try_limit
          } else {
            max(x_end, try_limit)
          }
          break
        }
      }
    }, silent = TRUE)
  }

  # Ensure valid grid (minimum width)
  if (x_end <= x_start || !is.finite(x_end)) {
    x_end <- max(x_start + 3, 3)
  }

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
      "Cp"  = tol / ((2 * sigma_level) * sig_samples),
      "Cpu" = (USL - mu_samples) / (sigma_level * sig_samples),
      "Cpl" = (mu_samples - LSL) / (sigma_level * sig_samples),
      "Cpk" = pmin((USL - mu_samples), (mu_samples - LSL)) /
        (sigma_level * sig_samples),
      "Cpm" = pmin(USL - target, target - LSL) /
              (sigma_level * sqrt(sig_samples^2 + (mu_samples - target)^2)),
      "Cpc" = {
        z <- (mu_samples - target) / sig_samples
        E_abs_dev <- sig_samples * sqrt(2 / pi) * exp(-0.5 * z^2) +
                     abs(mu_samples - target) * (1 - 2 * stats::pnorm(-abs(z)))
        tol / ((2 * sigma_level) * sqrt(pi / 2) * E_abs_dev)
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
                     metric %in% density_metrics &&
                     !metric_can_be_negative

  if (can_use_density) {
    # Direct PDF computation via contour integration
    pdf_fn <- make_density_solver(data, LSL, USL, prior, metric, target,
                                  sigma_level = sigma_level,
                                  cached_state = cached_state)

    result <- .adaptive_density_grid(pdf_fn, x_start, x_end, n_grid,
                                     is_prior_only = is_prior_only)
    grid_x   <- result$grid_x
    mid_x    <- result$mid_x
    pdf_vals <- result$pdf_vals
    dx       <- result$dx
    area     <- result$area

  } else {
    # Fallback: Survival function + finite differences
    if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
      S <- make_solver(data, LSL, USL, prior, metric, target,
                       sigma_level = sigma_level,
                       cached_state = cached_state)
    } else {
      S <- make_solver(data, LSL, USL, prior, metric, target,
                       sigma_level = sigma_level)
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
    neg_correction <- .compute_negative_mean_correction(
      metric, LSL, USL, bt, sigma_level = sigma_level
    )
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

#' Extract effective support bounds for a BayesTools prior
#' @keywords internal
.extract_prior_bounds <- function(bt_prior) {
  if (!inherits(bt_prior, "prior")) {
    return(list(lower = -Inf, upper = Inf))
  }

  dist <- bt_prior[["distribution"]]
  params <- bt_prior[["parameters"]]
  trunc <- bt_prior[["truncation"]]

  lower <- trunc[["lower"]] %||% -Inf
  upper <- trunc[["upper"]] %||% Inf

  natural_bounds <- switch(dist,
    "uniform" = list(lower = params[["a"]], upper = params[["b"]]),
    "gamma" = list(lower = 0, upper = Inf),
    "invgamma" = list(lower = 0, upper = Inf),
    "lognormal" = list(lower = 0, upper = Inf),
    "exp" = list(lower = 0, upper = Inf),
    "beta" = list(lower = 0, upper = 1),
    "point" = list(lower = params[["location"]], upper = params[["location"]]),
    list(lower = -Inf, upper = Inf)
  )

  list(
    lower = max(lower, natural_bounds$lower),
    upper = min(upper, natural_bounds$upper)
  )
}

#' Find finite optimizer initials for generic integration problems
#' @keywords internal
.find_feasible_generic_init <- function(init_mu, init_sd, prior, log_post) {
  candidate_pairs <- list(c(init_mu, max(init_sd, 1e-6)))

  bt <- prior$bayestools_priors
  if (!is.null(bt)) {
    mu_init  <- .extract_prior_init(bt$mu)
    sig_init <- .extract_prior_init(bt$sigma)

    mu_candidates <- unique(c(init_mu, mu_init$value, 0))
    sig_candidates <- unique(c(init_sd, sig_init$value, 0.1, 1))

    if (inherits(bt$sigma, "prior")) {
      sig_q <- tryCatch(
        as.numeric(BayesTools::quant(bt$sigma, c(0.25, 0.5, 0.75))),
        error = function(e) numeric(0)
      )
      sig_candidates <- c(sig_candidates, sig_q)
    }

    mu_candidates <- mu_candidates[is.finite(mu_candidates)]
    sig_candidates <- sig_candidates[is.finite(sig_candidates) & sig_candidates > 0]

    for (sigma in sig_candidates) {
      for (mu in mu_candidates) {
        candidate_pairs[[length(candidate_pairs) + 1L]] <- c(mu, sigma)
      }
    }
  }

  candidate_pairs <- c(candidate_pairs, list(c(0, max(init_sd, 1e-6)), c(0, 1)))

  for (candidate in candidate_pairs) {
    value <- tryCatch(log_post(candidate[1], candidate[2]), error = function(e) -Inf)
    if (is.finite(value)) {
      return(list(mu = candidate[1], sigma = candidate[2]))
    }
  }

  stop(
    "Could not find finite initial values for integration with the supplied priors. ",
    "Consider widening the prior support or using method = 'mcmc'."
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
                                   sigma = 3, sample_priors = FALSE,
                                   cached_state = NULL) {

  # If sampling from priors, ignore data (likelihood becomes flat/unity effectively)
  # Ideally we should pass sample_priors down, but clearing data works if the
  # functions handle empty data correctly.
  # However, the conjugate update uses n and sufficiency stats.
  # If we set data to empty, n=0, and the posterior parameters equal prior parameters.
  if (sample_priors) {
    data <- numeric(0)
    cached_state <- NULL
  }

  # Convert priors to integration format (4-case classification).
  # A PriorConjugate passed as prior_mu (e.g. from create_prior_unit_information)
  # is already in integration format and does not need conversion.
  prior_info <- if (inherits(prior_mu, "PriorConjugate")) {
    list(prior = prior_mu, is_conjugate = TRUE, case = 1L)
  } else {
    .bayestools_to_integration_prior(prior_mu, prior_sigma)
  }
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
  } else if (case == 4L && is.null(cached_state$log_post_vec)) {
    cached_state <- precompute_generic_state(data, prior, cached_state = cached_state)
  }

  n <- cached_state$n %||% 0L

  # Analyze all metrics
  # Match order of MCMC results for consistency (Cp, Cpu, Cpl, Cpk, Cpc, Cpm)
  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")

  # Pre-compute analytic divergence info for each metric.
  # With data (n > 0), the likelihood provides superexponential decay at sigma = 0,
  # so all moments are finite regardless of the prior.
  is_prior_only <- n == 0
  if (is_conjugate) {
    alpha_sigma <- prior$alpha0
    sigma_label <- sprintf("NIG(alpha0 = %s)", format(prior$alpha0, digits = 3))
  } else {
    alpha_sigma <- .extract_alpha_parameter(prior_sigma)
    sigma_label <- .prior_sigma_label(prior_sigma)
  }
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
                                   divergence_info = divergence_map[[m]],
                                   sigma_level = sigma)
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
  attr(metrics_list, "sigma") <- sigma
  attr(metrics_list, "method") <- "integration"

  list(
    results = results,
    metrics = metrics_list,
    coefficients = coefficients,
    sigma = sigma,
    prior = prior,
    is_conjugate = is_conjugate,
    case = case,
    cached_state = cached_state,
    divergence = divergence_map
  )
}

