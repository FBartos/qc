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

.integration_interval_from_entry <- function(entry, ci, ci_level, x = NULL, density = NULL) {
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
