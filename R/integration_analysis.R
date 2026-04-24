# ==============================================================================
# Integration Analysis and Probability Utilities
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
                                          sigma_level = 3,
                                          request = NULL) {
  request <- .as_qc_integration_request(
    request = request,
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  metric <- request$metric

  case_info <- .integration_metric_case(
    request = request,
    context = "probability"
  )

  if (!is.null(case_info$degenerate)) {
    return(.degenerate_metric_prob(case_info$degenerate, bounds))
  }

  backend <- .integration_backend_resolver(
    request = request,
    prefer_density = TRUE
  )

  if (backend$can_use_density) {
    prob <- tryCatch({
      pdf_vec <- function(x) vapply(x, backend$pdf_fn, numeric(1L))

      # Capability indices on this path have support on [0, Inf).
      lower <- max(min(bounds), 0)
      upper <- max(bounds)
      if (upper <= lower) return(0)

      stats::integrate(pdf_vec, lower, upper, rel.tol = 1e-4)$value
    }, error = function(e) NA_real_)
    if (!is.na(prob)) return(prob)

    backend$S <- tryCatch(
      .integration_make_solver(request = request),
      error = function(e) NULL
    )
    if (!is.null(backend$S)) {
      prob <- tryCatch(
        .integration_interval_prob_from_survival(backend$S, bounds),
        error = function(e) NA_real_
      )
      if (!is.na(prob)) return(prob)
    }

    return(NA_real_)
  }

  if (is.null(backend$S)) {
    return(NA_real_)
  }

  .integration_interval_prob_from_survival(backend$S, bounds)
}

.integration_extract_bayestools_priors <- function(prior, cached_state = NULL) {
  bt <- prior$bayestools_priors %||% NULL
  if (is.null(bt) && !is.null(cached_state)) {
    bt <- cached_state$bayestools_priors %||% NULL
  }
  bt
}

.integration_make_sigma_log_dens_fn <- function(sigma_prior) {
  if (inherits(sigma_prior, "PriorConjugate")) {
    alpha0 <- sigma_prior$alpha0
    beta0 <- sigma_prior$beta0

    if (is.finite(alpha0) && alpha0 > 0 && is.finite(beta0) && beta0 > 0) {
      log_const <- log(2) + alpha0 * log(beta0) - lgamma(alpha0)
      return(function(sigma) {
        ifelse(
          sigma <= 0,
          -Inf,
          log_const - (2 * alpha0 + 1) * log(sigma) - beta0 / sigma^2
        )
      })
    }

    if (identical(alpha0, -0.5) && identical(beta0, 0)) {
      return(function(sigma) ifelse(sigma <= 0, -Inf, -log(sigma)))
    }

    return(function(sigma) {
      ifelse(
        sigma <= 0,
        -Inf,
        -(2 * alpha0 + 1) * log(sigma) - beta0 / sigma^2
      )
    })
  }

  .make_prior_log_dens_fn(sigma_prior)
}

.integration_can_use_density <- function(prior, metric) {
  if (.metric_can_be_negative(metric)) {
    return(FALSE)
  }

  if (inherits(prior, "PriorSemiConjugateSigma") &&
      metric %in% c("Cpm", "Cpc")) {
    return(FALSE)
  }

  if (inherits(prior, "PriorGeneric") && identical(metric, "Cpc")) {
    return(FALSE)
  }

  inherits(prior, "PriorConjugate") ||
    inherits(prior, "PriorGeneric") ||
    inherits(prior, "PriorSemiConjugateMu") ||
    inherits(prior, "PriorSemiConjugateSigma")
}

.integration_divergent_mean_value <- function(metric, approx_mean,
                                              mu_prior = NULL,
                                              LSL = NULL,
                                              USL = NULL) {
  if (!.metric_can_be_negative(metric)) {
    return(Inf)
  }

  if (inherits(mu_prior, "prior") && !is.null(LSL) && !is.null(USL)) {
    mu_bounds <- .extract_prior_bounds(mu_prior)

    if (metric == "Cpu") {
      if (is.finite(mu_bounds$lower) && mu_bounds$lower >= USL) return(-Inf)
      if (is.finite(mu_bounds$upper) && mu_bounds$upper <= USL) return(Inf)
    }

    if (metric == "Cpl") {
      if (is.finite(mu_bounds$upper) && mu_bounds$upper <= LSL) return(-Inf)
      if (is.finite(mu_bounds$lower) && mu_bounds$lower >= LSL) return(Inf)
    }

    if (metric == "Cpk") {
      if (is.finite(mu_bounds$lower) && mu_bounds$lower >= USL) return(-Inf)
      if (is.finite(mu_bounds$upper) && mu_bounds$upper <= LSL) return(-Inf)
      if (is.finite(mu_bounds$lower) && is.finite(mu_bounds$upper) &&
          mu_bounds$lower >= LSL && mu_bounds$upper <= USL) {
        return(Inf)
      }
    }
  }

  if (is.finite(approx_mean) && approx_mean < 0) {
    return(-Inf)
  }

  Inf
}

.analyze_prior_only_mc_samples <- function(mc_samples, metric, LSL, USL, target,
                                           divergence_info, sigma_level = 3) {
  mu_samples <- mc_samples$mu
  sig_samples <- mc_samples$sig
  n_mc <- length(mu_samples)

  metric_samples <- compute_metric_value(
    mu = mu_samples,
    sigma = sig_samples,
    LSL = LSL,
    USL = USL,
    target = target,
    metric = metric,
    sigma_level = sigma_level
  )

  post_mean <- mean(metric_samples)
  post_median <- stats::median(metric_samples)
  post_sd <- stats::sd(metric_samples)
  q2.5 <- unname(stats::quantile(metric_samples, 0.025))
  q97.5 <- unname(stats::quantile(metric_samples, 0.975))

  sorted_samples <- sort(metric_samples)
  n_ci <- floor(0.95 * n_mc)
  ci_widths <- sorted_samples[(n_ci + 1):n_mc] - sorted_samples[1:(n_mc - n_ci)]
  best_ci <- which.min(ci_widths)
  hdi_lo <- sorted_samples[best_ci]
  hdi_hi <- sorted_samples[best_ci + n_ci]

  if (divergence_info$mean_divergent) {
    post_mean <- .integration_divergent_mean_value(metric, post_mean)
    post_sd <- Inf
  } else if (divergence_info$sd_divergent) {
    post_sd <- Inf
  }

  list(
    metric = metric,
    samples = metric_samples,
    area = 1,
    stats = c(Mean = post_mean, Median = post_median, SD = post_sd,
              Q2.5 = q2.5, Q97.5 = q97.5,
              HDI_Lo = hdi_lo, HDI_Hi = hdi_hi),
    divergence_info = divergence_info
  )
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
  log_dens_sigma <- .integration_make_sigma_log_dens_fn(bayestools_priors$sigma)

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

.integration_density_tail_masses <- function(pdf_fn, grid_x, grid_mass,
                                             support_lower = 0) {
  left_tail_mass <- 0

  if (length(grid_x) >= 2L && is.finite(support_lower) && grid_x[1] > support_lower) {
    pdf_vec <- function(x) vapply(x, pdf_fn, numeric(1L))
    left_tail_mass <- tryCatch(
      stats::integrate(pdf_vec, support_lower, grid_x[1], rel.tol = 1e-4)$value,
      error = function(e) NA_real_
    )

    if (!is.finite(left_tail_mass) || left_tail_mass < 0) {
      left_tail_mass <- 0
    }
  }

  list(
    mass_nonpositive = 0,
    left_tail_mass = left_tail_mass,
    right_tail_mass = max(0, 1 - left_tail_mass - grid_mass),
    support_lower = support_lower
  )
}

.integration_survival_tail_masses <- function(S, grid_x, grid_mass,
                                              support_lower = 0,
                                              metric_can_be_negative = FALSE) {
  clamp_prob <- function(value) {
    if (!is.finite(value)) {
      return(0)
    }
    max(min(value, 1), 0)
  }

  if (!length(grid_x)) {
    return(list(
      mass_nonpositive = 0,
      left_tail_mass = 0,
      right_tail_mass = max(0, 1 - grid_mass),
      support_lower = support_lower
    ))
  }

  start_prob <- clamp_prob(tryCatch(S(grid_x[1]), error = function(e) NA_real_))
  end_prob <- clamp_prob(tryCatch(S(grid_x[length(grid_x)]), error = function(e) NA_real_))

  mass_nonpositive <- 0
  left_tail_mass <- 0
  support_floor <- if (metric_can_be_negative) grid_x[1] else support_lower

  if (metric_can_be_negative && grid_x[1] > 0) {
    zero_prob <- clamp_prob(tryCatch(S(0), error = function(e) NA_real_))
    mass_nonpositive <- max(0, 1 - zero_prob)
    left_tail_mass <- max(0, zero_prob - start_prob)
    support_floor <- 0
  } else if (!metric_can_be_negative &&
             is.finite(support_lower) &&
             grid_x[1] > support_lower) {
    left_tail_mass <- max(0, 1 - start_prob)
    support_floor <- support_lower
  }

  right_tail_mass <- if (is.finite(grid_x[length(grid_x)])) {
    end_prob
  } else {
    0
  }
  known_mass <- mass_nonpositive + left_tail_mass + grid_mass + right_tail_mass
  if (is.finite(known_mass) && abs(known_mass - 1) > 1e-8) {
    right_tail_mass <- max(0, right_tail_mass + (1 - known_mass))
  }

  list(
    mass_nonpositive = mass_nonpositive,
    left_tail_mass = left_tail_mass,
    right_tail_mass = right_tail_mass,
    support_lower = support_floor
  )
}

.integration_moments_from_grid_entry <- function(entry) {
  grid_info <- .integration_density_grid_info(
    entry = entry,
    x = entry$grid,
    density = entry$grid$density,
    area = entry$area %||% 1
  )
  if (is.null(grid_info)) {
    return(list(mean = NA_real_, sd = NA_real_))
  }

  add_uniform_segment <- function(mass, lower, upper) {
    if (!is.finite(mass) || mass <= 0) {
      return(c(mean = 0, second = 0))
    }

    if (!is.finite(lower) || !is.finite(upper) || upper < lower) {
      point <- if (is.finite(lower)) lower else upper
      if (!is.finite(point)) {
        return(c(mean = Inf, second = Inf))
      }
      return(c(mean = mass * point, second = mass * point^2))
    }

    c(
      mean = mass * (lower + upper) / 2,
      second = mass * (lower^2 + lower * upper + upper^2) / 3
    )
  }

  mean_term <- sum(grid_info$x * grid_info$abs_mass)
  second_term <- sum((grid_info$x^2) * grid_info$abs_mass)

  if (grid_info$left_tail_mass > 0) {
    left_tail <- add_uniform_segment(
      grid_info$left_tail_mass,
      grid_info$support_lower,
      grid_info$left[1]
    )
    mean_term <- mean_term + left_tail[["mean"]]
    second_term <- second_term + left_tail[["second"]]
  }

  if (grid_info$right_tail_mass > 0) {
    right_tail <- add_uniform_segment(
      grid_info$right_tail_mass,
      tail(grid_info$right, 1),
      .integration_right_tail_endpoint(grid_info)
    )
    mean_term <- mean_term + right_tail[["mean"]]
    second_term <- second_term + right_tail[["second"]]
  }

  if (!is.finite(mean_term) || !is.finite(second_term)) {
    return(list(mean = mean_term, sd = Inf))
  }

  variance <- second_term - mean_term^2
  list(mean = mean_term, sd = sqrt(max(0, variance)))
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
                                            sigma_level = 3,
                                            request = NULL) {
  request <- .as_qc_integration_request(
    request = request,
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  data <- request$data
  LSL <- request$LSL
  USL <- request$USL
  prior <- request$prior
  metric <- request$metric
  target <- request$target
  cached_state <- request$cached_state
  sigma_level <- request$sigma_level

  case_info <- .integration_metric_case(
    request = request,
    context = "analysis"
  )

  # Default: no divergence
  if (is.null(divergence_info))
    divergence_info <- list(mean_divergent = FALSE, sd_divergent = FALSE,
                            alpha = Inf, reason = NULL)

  # Check if we're in prior-only mode (no data) with PriorGeneric
  # In this case, use quantile-based grid bounds to avoid slow moment computation
  is_prior_only <- .integration_is_prior_only(data = data, cached_state = cached_state)
  bt_priors <- .integration_extract_bayestools_priors(prior, cached_state)
  has_sigma_quantile_prior <- !is.null(bt_priors) &&
    inherits(bt_priors$sigma %||% NULL, "prior")
  metric_can_be_negative <- case_info$metric_can_be_negative
  moments <- NULL

  if (is_prior_only && !is.null(mc_samples)) {
    return(.analyze_prior_only_mc_samples(
      mc_samples = mc_samples,
      metric = metric,
      LSL = LSL,
      USL = USL,
      target = target,
      divergence_info = divergence_info,
      sigma_level = sigma_level
    ))
  }

  if (!is.null(case_info$degenerate)) {
    return(.analyze_degenerate_metric_distribution(
      metric,
      case_info$degenerate,
      n_grid,
      divergence_info = divergence_info,
      metric_can_be_negative = metric_can_be_negative
    ))
  }

  if (is_prior_only && has_sigma_quantile_prior) {
    # Fast path: compute grid bounds directly from sigma quantiles
    bounds <- .compute_metric_grid_bounds_from_quantiles(
      bt_priors$sigma, LSL, USL, metric, target,
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
    moments <- .integration_compute_metric_moments(
      request = request,
      use_analytic = TRUE
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
      S_fn <- .integration_make_solver(request = request)

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
  # Choose method: density solver (direct PDF) vs survival function + finite diff
  # Density solver works for both PriorConjugate and PriorGeneric with supported metrics
  backend <- .integration_backend_resolver(
    request = request,
    prefer_density = use_density_solver
  )

  if (backend$can_use_density) {
    result <- .adaptive_density_grid(
      backend$pdf_fn,
      x_start,
      x_end,
      n_grid,
      is_prior_only = is_prior_only
    )
    grid_x <- result$grid_x
    mid_x <- result$mid_x
    pdf_vals <- result$pdf_vals
    dx <- result$dx
    area <- result$area
  } else {
    # Fallback: Survival function + finite differences
    if (is.null(backend$S)) {
      stop("Failed to construct an integration backend for metric '", metric, "'.")
    }

    # Evaluate grid
    grid_x <- seq(x_start, x_end, length.out = n_grid)
    S_vals <- sapply(grid_x, backend$S)

    # Handle NaN in S_vals
    S_vals[!is.finite(S_vals)] <- 0

    # Compute PDF via finite differences
    pdf_vals <- -diff(S_vals) / diff(grid_x)
    mid_x <- (grid_x[-1] + grid_x[-n_grid]) / 2
    dx <- diff(grid_x)

    # Handle NaN/negative in PDF
    pdf_vals[!is.finite(pdf_vals) | pdf_vals < 0] <- 0

    # Capture the realized mass on the finite survival grid.
    area <- sum(pdf_vals * dx)
    if (area > 0) pdf_vals <- pdf_vals / area
  }

  grid_left <- grid_x[-length(grid_x)]
  grid_right <- grid_x[-1L]
  grid_df <- data.frame(
    x = mid_x,
    density = pdf_vals,
    x_left = grid_left,
    x_right = grid_right,
    width = dx
  )

  grid_entry <- .integration_grid_entry(
    metric = metric,
    grid_df = grid_df,
    area = area,
    backend = backend,
    metric_can_be_negative = metric_can_be_negative
  )

  # Reuse the moment engine whenever it already produced a stable answer for
  # this posterior. The grid/tail reconstruction is kept as a fallback for
  # cases where we skipped or could not trust the moment computation.
  moment_info <- .integration_moments_from_grid_entry(grid_entry)
  post_mean <- if (!is.null(moments) && is.finite(moments$mean)) {
    moments$mean
  } else {
    moment_info$mean
  }
  post_sd <- if (!is.null(moments) &&
                 is.finite(moments$sd) &&
                 moments$sd >= 0) {
    moments$sd
  } else {
    moment_info$sd
  }
  quantiles <- .integration_quantiles_from_density(
    x = grid_df,
    density = grid_df$density,
    probs = c(0.025, 0.5, 0.975),
    area = grid_entry$area,
    entry = grid_entry
  )
  hdi <- .integration_hdi_from_density(
    x = grid_df,
    density = grid_df$density,
    ci_level = 0.95,
    area = grid_entry$area,
    entry = grid_entry
  )

  # Override with analytic Inf where divergent
  if (divergence_info$mean_divergent) {
    post_mean <- .integration_divergent_mean_value(
      metric,
      post_mean,
      mu_prior = bt_priors$mu %||% NULL,
      LSL = LSL,
      USL = USL
    )
    post_sd <- Inf
  } else if (divergence_info$sd_divergent) {
    post_sd <- Inf
  }

  grid_entry$stats <- c(
    Mean = post_mean,
    Median = quantiles[2],
    SD = post_sd,
    Q2.5 = quantiles[1],
    Q97.5 = quantiles[3],
    HDI_Lo = hdi[1],
    HDI_Hi = hdi[2]
  )
  grid_entry$divergence_info <- divergence_info
  grid_entry
}


# ==============================================================================
# BayesTools Prior Conversion
# ==============================================================================

# ==============================================================================
# Analytic Moment Divergence Detection
# ==============================================================================

