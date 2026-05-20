# ==============================================================================
# Conjugate Posterior Helpers and Degenerate Metric Handling
# ==============================================================================

# ==============================================================================
# Sufficient Statistics & Posterior Update Helpers
# ==============================================================================

.extract_suff_stats <- function(data, cached_state) {
  suff_state <- .integration_require_suff_state(
    data = data,
    cached_state = cached_state
  )

  list(
    n = suff_state$n,
    x_bar = suff_state$x_bar,
    SS = suff_state$sse
  )
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

.format_conjugate_posterior_parameter <- function(x) {
  format(signif(x, 6), trim = TRUE, scientific = FALSE)
}

.conjugate_posterior_error_message <- function(k_n, alpha_n, beta_n,
                                               context = "The conjugate posterior") {
  sprintf(
    "%s is improper for the supplied data and prior (k_n = %s, alpha_n = %s, beta_n = %s).",
    context,
    .format_conjugate_posterior_parameter(k_n),
    .format_conjugate_posterior_parameter(alpha_n),
    .format_conjugate_posterior_parameter(beta_n)
  )
}

.compute_validated_conjugate_posterior <- function(prior, data, cached_state = NULL,
                                                   context = "The conjugate posterior") {
  ss <- .extract_suff_stats(data, cached_state)
  post <- .nig_posterior(prior, ss$n, ss$x_bar, ss$SS)

  if (.is_improper_conjugate_posterior(post$k_n, post$alpha_n, post$beta_n)) {
    stop(
      .conjugate_posterior_error_message(
        post$k_n, post$alpha_n, post$beta_n,
        context = context
      ),
      call. = FALSE
    )
  }

  list(
    ss = ss,
    post = post,
    is_degenerate = .is_degenerate_conjugate_posterior(post$beta_n)
  )
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
    value <- .cpm_spec_distance(LSL, USL, target) / (sigma_level * delta)
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
    # Match cut(..., include.lowest = TRUE): bins are left-open/right-closed
    # after the first interval, so boundary point masses are not dropped.
    "point" = as.numeric(lower < dist$value && dist$value <= upper),
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

.semi_mu_posterior <- function(prior, n, x_bar, sse,
                               context = "The semi-conjugate mu posterior") {
  k_n   <- prior$k0 + n
  if (!is.finite(k_n) || k_n <= 0) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (k_n = %s).",
        context,
        .format_conjugate_posterior_parameter(k_n)
      ),
      call. = FALSE
    )
  }
  mu_n  <- (prior$k0 * prior$mu0 + n * x_bar) / k_n
  sse_n <- sse + prior$k0 * n * (x_bar - prior$mu0)^2 / k_n
  list(k_n = k_n, mu_n = mu_n, sse_n = sse_n)
}
