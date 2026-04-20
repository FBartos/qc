.integration_result_stats <- function(x) {
  if (is.list(x) && !is.null(x$stats)) {
    return(x$stats)
  }

  x
}

.integration_grid_cell_widths <- function(x) {
  n <- length(x)

  if (n == 0L) {
    return(numeric(0))
  }

  if (n == 1L) {
    return(1)
  }

  edges <- c(
    x[1] - (x[2] - x[1]) / 2,
    (x[-n] + x[-1L]) / 2,
    x[n] + (x[n] - x[n - 1L]) / 2
  )

  diff(edges)
}

.integration_quantiles_from_density <- function(x, density, probs, area = 1) {
  dx <- .integration_grid_cell_widths(x)
  if (length(dx) == 0L) {
    return(rep(NA_real_, length(probs)))
  }

  mass <- density * dx
  total_mass <- sum(mass)
  if (!is.finite(total_mass) || total_mass <= 0) {
    return(rep(NA_real_, length(probs)))
  }

  cdf_cond <- cumsum(mass) / total_mass
  if (is.null(area) || !is.finite(area)) {
    area <- 1
  }
  area <- max(min(area, 1), 0)
  cdf_adj <- (1 - area) + area * cdf_cond

  vapply(probs, function(prob) {
    if (!is.finite(prob)) {
      return(NA_real_)
    }
    if (prob <= 1 - area) {
      return(0)
    }

    x[which.min(abs(cdf_adj - prob))]
  }, numeric(1))
}

.integration_hdi_from_density <- function(x, density, ci_level) {
  dx <- .integration_grid_cell_widths(x)
  if (length(dx) == 0L) {
    return(c(NA_real_, NA_real_))
  }

  mass <- density * dx
  total_mass <- sum(mass)
  if (!is.finite(total_mass) || total_mass <= 0) {
    return(c(NA_real_, NA_real_))
  }

  ci_level <- max(min(ci_level, 1), 0)
  mass <- mass / total_mass

  sorted_idx <- order(density, decreasing = TRUE)
  cum_mass <- cumsum(mass[sorted_idx])
  cutoff_idx <- which(cum_mass >= ci_level)[1]

  if (is.na(cutoff_idx)) {
    cutoff_idx <- length(sorted_idx)
  }

  hdi_idx <- sorted_idx[seq_len(cutoff_idx)]
  c(min(x[hdi_idx]), max(x[hdi_idx]))
}

.integration_hdi_from_samples <- function(x, ci_level) {
  x <- sort(x[is.finite(x)])
  n <- length(x)

  if (n == 0L) {
    return(c(NA_real_, NA_real_))
  }

  ci_level <- max(min(ci_level, 1), 0)
  n_ci <- floor(ci_level * n)

  if (n_ci >= n) {
    return(c(x[1], x[n]))
  }

  ci_widths <- x[(n_ci + 1):n] - x[seq_len(n - n_ci)]
  best_ci <- which.min(ci_widths)

  c(x[best_ci], x[best_ci + n_ci])
}

.integration_interval_probs_from_grid <- function(entry, interval_breaks) {
  if (!is.list(entry) || is.null(entry$grid)) {
    return(NULL)
  }

  g <- entry$grid
  if (!is.data.frame(g) || !all(c("x", "density") %in% names(g)) || nrow(g) < 2L) {
    return(NULL)
  }

  density <- g$density
  if (!any(is.finite(density) & density > 0)) {
    return(NULL)
  }

  area <- entry$area %||% 1
  if (!is.finite(area)) {
    area <- 1
  }
  area <- max(min(area, 1), 0)

  n_g <- nrow(g)
  dxx <- diff(g$x)
  avg_dens <- (density[-n_g] + density[-1L]) / 2
  cdf_cond <- c(0, cumsum(avg_dens * dxx))
  total_mass <- cdf_cond[length(cdf_cond)]
  if (!is.finite(total_mass) || total_mass <= 0) {
    return(NULL)
  }
  cdf_cond <- cdf_cond / total_mass

  n_intervals <- length(interval_breaks) - 1L
  probs <- numeric(n_intervals)
  grid_min <- min(g$x)
  grid_max <- max(g$x)

  for (j in seq_len(n_intervals)) {
    lo <- interval_breaks[j]
    hi <- interval_breaks[j + 1L]

    lo_cdf <- if (is.finite(lo) && lo >= grid_min) {
      (1 - area) + area * stats::approx(g$x, cdf_cond, xout = lo, rule = 2)$y
    } else {
      0
    }

    hi_cdf <- if (is.finite(hi) && hi <= grid_max) {
      (1 - area) + area * stats::approx(g$x, cdf_cond, xout = hi, rule = 2)$y
    } else {
      1
    }

    probs[j] <- max(0, hi_cdf - lo_cdf)
  }

  total <- sum(probs)
  if (total > 0) {
    probs <- probs / total
  }

  probs
}

.integration_interval_probs_from_entry <- function(entry, interval_breaks,
                                                   metric, prior, cached_state,
                                                   LSL, USL, target,
                                                   sigma_level = 3) {
  n_intervals <- length(interval_breaks) - 1L

  if (is.list(entry) && !is.null(entry$degenerate)) {
    probs <- vapply(seq_len(n_intervals), function(j) {
      .degenerate_metric_prob(entry$degenerate, interval_breaks[c(j, j + 1L)])
    }, numeric(1))
  } else if (is.list(entry) && !is.null(entry$samples)) {
    bin_counts <- table(cut(entry$samples,
                            breaks = interval_breaks, include.lowest = TRUE))
    probs <- as.numeric(bin_counts) / sum(bin_counts)
  } else {
    grid_probs <- .integration_interval_probs_from_grid(entry, interval_breaks)

    probs <- vapply(seq_len(n_intervals), function(j) {
      lo <- interval_breaks[j]
      hi <- interval_breaks[j + 1L]

      tryCatch(
        compute_cpk_prob_integration(
          numeric(0), LSL, USL, c(lo, hi), prior,
          metric = metric, target = target,
          cached_state = cached_state,
          sigma_level = sigma_level
        ),
        error = function(e) NA_real_
      )
    }, numeric(1))

    use_integration <- !anyNA(probs)

    # Placeholder grids for degenerate +/-Inf results carry no finite density
    # information, so only use the grid when it contains actual mass.
    if (use_integration && !is.null(grid_probs)) {
      suspect <- any(grid_probs > 0.05 & probs < 0.001)
      if (suspect) {
        use_integration <- FALSE
      }
    } else if (!use_integration && is.null(grid_probs)) {
      probs[is.na(probs)] <- 0
      use_integration <- TRUE
    }

    if (!use_integration) {
      probs <- grid_probs
    }
  }

  total <- sum(probs)
  if (total > 0) {
    probs <- probs / total
  }

  probs
}

.integration_interval_from_entry <- function(entry, ci, ci_level, x = NULL, density = NULL) {
  if (is.list(entry) && !is.null(entry$degenerate)) {
    return(.degenerate_metric_interval(entry$degenerate, ci, ci_level))
  }

  stats <- .integration_result_stats(entry)
  if (!is.null(stats) && isTRUE(all.equal(ci_level, 0.95))) {
    return(switch(ci,
      "central" = c(stats["Q2.5"], stats["Q97.5"]),
      "HPD" = c(stats["HDI_Lo"], stats["HDI_Hi"]),
      stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
    ))
  }

  if (is.list(entry) && !is.null(entry$samples)) {
    return(switch(ci,
      "central" = {
        h <- (1 - ci_level) / 2
        unname(stats::quantile(entry$samples, c(h, 1 - h), na.rm = TRUE))
      },
      "HPD" = .integration_hdi_from_samples(entry$samples, ci_level),
      stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
    ))
  }

  if (is.list(entry) && !is.null(entry$grid)) {
    x <- entry$grid$x
    density <- entry$grid$density
    area <- entry$area
  } else {
    area <- 1
    if (is.list(entry) && !is.null(entry$area)) {
      area <- entry$area
    }
  }

  if (!is.null(x) && !is.null(density)) {
    return(switch(ci,
      "central" = {
        h <- (1 - ci_level) / 2
        .integration_quantiles_from_density(x, density, c(h, 1 - h), area = area)
      },
      "HPD" = .integration_hdi_from_density(x, density, ci_level),
      stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
    ))
  }

  switch(ci,
    "central" = c(stats["Q2.5"], stats["Q97.5"]),
    "HPD" = c(stats["HDI_Lo"], stats["HDI_Hi"]),
    stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
  )
}
