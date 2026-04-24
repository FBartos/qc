# ==============================================================================
# Metric Geometry and Constraint Helpers
# ==============================================================================

# Helper function for log-space difference of exponentials
log_diff_exp <- function(x, y) {
  ifelse(x <= y, -Inf, x + log1p(-exp(y - x)))
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

.cpm_spec_distance <- function(LSL, USL, target) {
  if (is.null(target)) target <- (LSL + USL) / 2
  min(USL - target, target - LSL)
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
      sigma_cap <- .cpm_spec_distance(LSL, USL, target) / (sigma_level * c)
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

  full_region <- lower <= 0 & is.infinite(upper) & upper > 0
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

  lower_only <- is.finite(lower) & lower > 0 & is.infinite(upper) & upper > 0
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
#' @param sigma_level Number of process standard deviations used to define the
#'   capability metric scale.
#' @return List with s_max_fn and mu_b_fn_vec (vectorized)
#' @export
get_metric_constraints <- function(metric, c, LSL, USL, target,
                                   sigma_level = 3) {
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = FALSE,
    sigma_name = "sigma_level"
  )

  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)

  list(
    s_max_fn = function() {
      if (metric %in% c("Cp", "Cpc"))
        return(if (c <= 0) Inf else tol / ((2 * sigma_level) * c))
      if (metric == "Cpm")
        return(if (c <= 0) Inf else cpm_dist / (sigma_level * c))
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
        R <- .cpm_spec_distance(LSL, USL, target) / (sigma_level * c)
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
#' @param sigma_level Number of process standard deviations used to define the
#'   capability metric scale.
#' @return Metric value(s) (same length as mu/sigma)
#' @export
compute_metric_value <- function(mu, sigma, LSL, USL, target, metric,
                                 sigma_level = 3) {
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = FALSE,
    sigma_name = "sigma_level"
  )

  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  switch(metric,
    "Cp"  = tol / ((2 * sigma_level) * sigma),
    "Cpu" = (USL - mu) / (sigma_level * sigma),
    "Cpl" = (mu - LSL) / (sigma_level * sigma),
    "Cpk" = pmin((USL - mu) / (sigma_level * sigma),
                 (mu - LSL) / (sigma_level * sigma)),
    "Cpm" = .cpm_spec_distance(LSL, USL, target) /
      (sigma_level * sqrt(sigma^2 + (mu - target)^2)),
    "Cpc" = tol / ((2 * sigma_level) * sqrt(pi / 2) *
                     .cpc_E_abs_dev_normal(mu, sigma, target)),
    stop("Unknown metric: ", metric)
  )
}

# ==============================================================================
# Posterior Moments of Capability Metrics
# ==============================================================================

