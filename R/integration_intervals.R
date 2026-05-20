.integration_result_stats <- function(x) {
  if (is.list(x) && !is.null(x$stats)) {
    return(x$stats)
  }

  x
}

.integration_result_degenerate <- function(x) {
  if (is.list(x) && !is.null(x$degenerate)) {
    return(x$degenerate)
  }

  NULL
}

.integration_result_mode <- function(entry, x = NULL, density = NULL) {
  degenerate <- .integration_result_degenerate(entry)
  if (!is.null(degenerate)) {
    return(switch(degenerate$type,
      "point" = degenerate$value,
      "normal" = degenerate$mean,
      "pos_inf" = Inf,
      "neg_inf" = -Inf,
      stop("Unknown degenerate metric distribution type: ", degenerate$type)
    ))
  }

  stats <- .integration_result_stats(entry)
  if (!is.null(stats)) {
    mean_val <- unname(stats["Mean"])
    median_val <- unname(stats["Median"])
    sd_val <- unname(stats["SD"])

    if (length(sd_val) == 1L && is.finite(sd_val) && sd_val == 0 &&
        length(median_val) == 1L && !is.na(median_val)) {
      return(median_val)
    }
  }

  if (!is.null(x) && !is.null(density)) {
    x <- as.numeric(x)
    density <- as.numeric(density)
    finite_peak_idx <- which(is.finite(x) & is.finite(density) & density > 0)

    if (length(finite_peak_idx) > 0L) {
      return(x[finite_peak_idx[which.max(density[finite_peak_idx])]])
    }
  }

  if (!is.null(stats)) {
    if (length(median_val) == 1L && !is.na(median_val) && !is.finite(median_val)) {
      return(median_val)
    }

    if (length(mean_val) == 1L && !is.na(mean_val) && !is.finite(mean_val)) {
      return(mean_val)
    }
  }

  if (!is.null(stats)) {
    return(unname(stats["Median"]))
  }

  NA_real_
}

.integration_grid_edges <- function(x) {
  if (is.null(x)) {
    return(list(left = numeric(0), right = numeric(0)))
  }

  if (is.data.frame(x)) {
    if (all(c("x_left", "x_right") %in% names(x))) {
      return(list(left = as.numeric(x$x_left), right = as.numeric(x$x_right)))
    }

    if (all(c("x", "width") %in% names(x))) {
      half_width <- as.numeric(x$width) / 2
      return(list(
        left = as.numeric(x$x) - half_width,
        right = as.numeric(x$x) + half_width
      ))
    }

    x <- x$x %||% numeric(0)
  } else if (is.list(x) && !is.atomic(x)) {
    if (!is.null(x$x_left) && !is.null(x$x_right)) {
      return(list(left = as.numeric(x$x_left), right = as.numeric(x$x_right)))
    }

    if (!is.null(x$x) && !is.null(x$width)) {
      half_width <- as.numeric(x$width) / 2
      return(list(
        left = as.numeric(x$x) - half_width,
        right = as.numeric(x$x) + half_width
      ))
    }

    x <- x$x %||% numeric(0)
  }

  x <- as.numeric(x)
  n <- length(x)

  if (n == 0L) {
    return(list(left = numeric(0), right = numeric(0)))
  }

  if (n == 1L) {
    return(list(left = x - 0.5, right = x + 0.5))
  }

  edges <- c(
    x[1] - (x[2] - x[1]) / 2,
    (x[-n] + x[-1L]) / 2,
    x[n] + (x[n] - x[n - 1L]) / 2
  )

  list(left = edges[-length(edges)], right = edges[-1L])
}

.integration_grid_cell_widths <- function(x) {
  edges <- .integration_grid_edges(x)
  edges$right - edges$left
}

.integration_density_grid_info <- function(entry = NULL, x = NULL,
                                           density = NULL, area = 1) {
  grid <- NULL
  if (is.list(entry) && !is.null(entry$grid)) {
    grid <- entry$grid
  } else if (is.data.frame(x) || (is.list(x) && !is.atomic(x))) {
    grid <- x
  }

  if (is.null(density) && !is.null(grid) && !is.null(grid$density)) {
    density <- grid$density
  }

  edges <- .integration_grid_edges(grid %||% x)
  width <- edges$right - edges$left

  if (length(width) == 0L || is.null(density)) {
    return(NULL)
  }

  density <- as.numeric(density)
  if (length(density) != length(width)) {
    return(NULL)
  }

  valid <- is.finite(density) & density >= 0 & is.finite(width) & width > 0
  if (!any(valid)) {
    return(NULL)
  }

  density[!valid] <- 0
  width[!valid] <- 0

  cond_mass <- density * width
  total_cond_mass <- sum(cond_mass)
  if (!is.finite(total_cond_mass) || total_cond_mass <= 0) {
    return(NULL)
  }

  if (is.list(entry)) {
    grid_mass <- entry$grid_mass %||% entry$area %||% area
    mass_nonpositive <- entry$mass_nonpositive
    left_tail_mass <- entry$left_tail_mass
    right_tail_mass <- entry$right_tail_mass
    support_lower <- entry$support_lower
  } else {
    grid_mass <- area
    mass_nonpositive <- 0
    left_tail_mass <- 0
    right_tail_mass <- max(0, 1 - grid_mass)
    support_lower <- NULL
  }

  if (is.null(grid_mass) || !is.finite(grid_mass)) {
    grid_mass <- 1
  }
  grid_mass <- max(min(grid_mass, 1), 0)

  if (!is.finite(mass_nonpositive %||% 0)) {
    mass_nonpositive <- 0
  }
  mass_nonpositive <- max(min(mass_nonpositive %||% 0, 1), 0)
  left_tail_mass <- left_tail_mass %||% 0
  right_tail_mass <- right_tail_mass %||%
    max(0, 1 - mass_nonpositive - left_tail_mass - grid_mass)

  if (!is.finite(left_tail_mass)) {
    left_tail_mass <- 0
  }
  if (!is.finite(right_tail_mass)) {
    right_tail_mass <- 0
  }

  left_tail_mass <- max(left_tail_mass, 0)
  right_tail_mass <- max(right_tail_mass, 0)

  known_mass <- mass_nonpositive + left_tail_mass + grid_mass + right_tail_mass
  if (is.finite(known_mass) && abs(known_mass - 1) > 1e-8) {
    right_tail_mass <- max(0, right_tail_mass + (1 - known_mass))
  }

  support_lower <- support_lower %||% min(edges$left[width > 0], na.rm = TRUE)
  abs_density <- density * (grid_mass / total_cond_mass)
  abs_mass <- abs_density * width

  list(
    x = if (!is.null(grid) && !is.null(grid$x)) as.numeric(grid$x) else as.numeric(x),
    density = density,
    left = edges$left,
    right = edges$right,
    width = width,
    abs_density = abs_density,
    abs_mass = abs_mass,
    grid_mass = grid_mass,
    mass_nonpositive = mass_nonpositive,
    left_tail_mass = left_tail_mass,
    right_tail_mass = right_tail_mass,
    support_lower = support_lower
  )
}

.integration_validate_interval_breaks <- function(interval_breaks) {
  interval_breaks <- as.numeric(interval_breaks)

  if (length(interval_breaks) < 2L) {
    stop("`interval_breaks` must contain at least two break points.")
  }

  if (anyNA(interval_breaks)) {
    stop("`interval_breaks` must not contain missing values.")
  }

  interior_idx <- seq.int(2L, length(interval_breaks) - 1L)
  if (length(interior_idx) > 0L &&
      any(!is.finite(interval_breaks[interior_idx]))) {
    stop("`interval_breaks` must have finite interior break points.")
  }

  if (any(diff(interval_breaks) <= 0)) {
    stop("`interval_breaks` must be strictly increasing.")
  }

  interval_breaks
}

.integration_left_tail_endpoint <- function(grid_info,
                                            tail_mass = grid_info$left_tail_mass) {
  tail_lo <- grid_info$support_lower
  tail_hi <- grid_info$left[1]

  if (tail_mass <= 0 || grid_info$left_tail_mass <= 0) {
    return(tail_lo %||% tail_hi)
  }

  tail_mass <- max(min(tail_mass, grid_info$left_tail_mass), 0)
  tail_width <- tail_hi - tail_lo

  if (!is.finite(tail_lo) || !is.finite(tail_hi) ||
      !is.finite(tail_width) || tail_width <= 0) {
    return(tail_hi)
  }

  tail_lo + tail_width * tail_mass / grid_info$left_tail_mass
}

.integration_right_tail_endpoint <- function(grid_info,
                                             tail_mass = grid_info$right_tail_mass) {
  tail_start <- grid_info$right[length(grid_info$right)]

  if (tail_mass <= 0 || grid_info$right_tail_mass <= 0) {
    return(tail_start)
  }

  tail_mass <- max(min(tail_mass, grid_info$right_tail_mass), 0)
  tail_density <- tail(grid_info$abs_density[
    is.finite(grid_info$abs_density) & grid_info$abs_density > 0
  ], 1)

  if (!length(tail_density) || !is.finite(tail_density) || tail_density <= 0) {
    return(Inf)
  }

  tail_start + tail_mass / tail_density
}

.integration_interval_index <- function(value, interval_breaks) {
  interval_breaks <- .integration_validate_interval_breaks(interval_breaks)

  idx <- cut(
    value,
    breaks = interval_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )

  if (is.na(idx)) {
    if (is.finite(value)) {
      return(NA_integer_)
    }

    return(if (value > 0) length(interval_breaks) - 1L else 1L)
  }

  as.integer(idx)
}

.integration_interval_overlap <- function(left, right, lo, hi) {
  overlap_left <- max(left, lo)
  overlap_right <- min(right, hi)
  max(0, overlap_right - overlap_left)
}

.integration_quantiles_from_density <- function(x, density, probs, area = 1,
                                                entry = NULL) {
  grid_info <- .integration_density_grid_info(
    entry = entry,
    x = x,
    density = density,
    area = area
  )
  if (is.null(grid_info)) {
    return(rep(NA_real_, length(probs)))
  }

  cum_grid_mass <- cumsum(grid_info$abs_mass)

  vapply(probs, function(prob) {
    if (!is.finite(prob)) {
      return(NA_real_)
    }

    prob <- max(min(prob, 1), 0)

    if (prob <= grid_info$mass_nonpositive) {
      return(0)
    }

    prob_after_nonpositive <- prob - grid_info$mass_nonpositive

    if (grid_info$left_tail_mass > 0 && prob_after_nonpositive <= grid_info$left_tail_mass) {
      return(.integration_left_tail_endpoint(grid_info, prob_after_nonpositive))
    }

    prob_in_grid <- prob_after_nonpositive - grid_info$left_tail_mass

    if (prob_in_grid <= 0 || grid_info$grid_mass <= 0) {
      return(grid_info$left[1])
    }

    if (prob_in_grid >= grid_info$grid_mass) {
      return(.integration_right_tail_endpoint(
        grid_info,
        prob_in_grid - grid_info$grid_mass
      ))
    }

    idx <- which(cum_grid_mass >= prob_in_grid)[1]
    prev_mass <- if (idx > 1L) cum_grid_mass[idx - 1L] else 0

    if (!is.finite(grid_info$abs_density[idx]) || grid_info$abs_density[idx] <= 0) {
      return(grid_info$x[idx])
    }

    grid_info$left[idx] + (prob_in_grid - prev_mass) / grid_info$abs_density[idx]
  }, numeric(1))
}

.integration_hdi_from_density <- function(x, density, ci_level, area = 1,
                                          entry = NULL) {
  grid_info <- .integration_density_grid_info(
    entry = entry,
    x = x,
    density = density,
    area = area
  )
  if (is.null(grid_info)) {
    return(c(NA_real_, NA_real_))
  }

  ci_level <- max(min(ci_level, 1), 0)

  if (ci_level <= grid_info$mass_nonpositive) {
    return(c(0, 0))
  }

  remaining_mass <- ci_level - grid_info$mass_nonpositive
  if (remaining_mass <= grid_info$left_tail_mass) {
    return(c(
      grid_info$support_lower %||% grid_info$left[1],
      .integration_left_tail_endpoint(grid_info, remaining_mass)
    ))
  }

  target_mass <- remaining_mass - grid_info$left_tail_mass
  tail_mass_needed <- max(0, target_mass - grid_info$grid_mass)
  target_mass <- min(target_mass, grid_info$grid_mass)
  sorted_idx <- order(grid_info$abs_density, decreasing = TRUE)
  cum_mass <- cumsum(grid_info$abs_mass[sorted_idx])
  cutoff_idx <- which(cum_mass >= target_mass)[1]

  if (is.na(cutoff_idx)) {
    cutoff_idx <- length(sorted_idx)
  }

  hdi_idx <- sorted_idx[seq_len(cutoff_idx)]
  upper <- max(grid_info$right[hdi_idx])
  if (tail_mass_needed > 0) {
    upper <- max(upper, .integration_right_tail_endpoint(grid_info, tail_mass_needed))
  }

  c(min(grid_info$left[hdi_idx]), upper)
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

  interval_breaks <- .integration_validate_interval_breaks(interval_breaks)

  g <- entry$grid
  if (!is.data.frame(g) || !all(c("x", "density") %in% names(g)) || nrow(g) < 1L) {
    return(NULL)
  }

  grid_info <- .integration_density_grid_info(
    entry = entry,
    x = g,
    density = g$density,
    area = entry$area %||% 1
  )
  if (is.null(grid_info)) {
    return(NULL)
  }

  n_intervals <- length(interval_breaks) - 1L
  probs <- numeric(n_intervals)

  if (grid_info$mass_nonpositive > 0) {
    zero_idx <- .integration_interval_index(0, interval_breaks)
    if (!is.na(zero_idx)) {
      probs[zero_idx] <- probs[zero_idx] + grid_info$mass_nonpositive
    }
  }

  if (grid_info$left_tail_mass > 0) {
    tail_lo <- grid_info$support_lower
    tail_hi <- grid_info$left[1]
    tail_width <- tail_hi - tail_lo

    if (is.finite(tail_lo) && is.finite(tail_hi) && tail_width > 0) {
      for (j in seq_len(n_intervals)) {
        overlap <- .integration_interval_overlap(
          tail_lo,
          tail_hi,
          interval_breaks[j],
          interval_breaks[j + 1L]
        )
        if (overlap > 0) {
          probs[j] <- probs[j] + grid_info$left_tail_mass * overlap / tail_width
        }
      }
    } else {
      tail_idx <- .integration_interval_index(tail_hi, interval_breaks)
      if (!is.na(tail_idx)) {
        probs[tail_idx] <- probs[tail_idx] + grid_info$left_tail_mass
      }
    }
  }

  for (i in seq_along(grid_info$abs_mass)) {
    if (!is.finite(grid_info$abs_density[i]) || grid_info$abs_density[i] <= 0) {
      next
    }

    for (j in seq_len(n_intervals)) {
      overlap <- .integration_interval_overlap(
        grid_info$left[i],
        grid_info$right[i],
        interval_breaks[j],
        interval_breaks[j + 1L]
      )
      if (overlap > 0) {
        probs[j] <- probs[j] + grid_info$abs_density[i] * overlap
      }
    }
  }

  if (grid_info$right_tail_mass > 0) {
    tail_lo <- tail(grid_info$right, 1)
    tail_hi <- .integration_right_tail_endpoint(grid_info)
    tail_width <- tail_hi - tail_lo

    if (is.finite(tail_lo) && is.finite(tail_hi) && tail_width > 0) {
      for (j in seq_len(n_intervals)) {
        overlap <- .integration_interval_overlap(
          tail_lo,
          tail_hi,
          interval_breaks[j],
          interval_breaks[j + 1L]
        )
        if (overlap > 0) {
          probs[j] <- probs[j] + grid_info$right_tail_mass * overlap / tail_width
        }
      }
    } else if (is.infinite(interval_breaks[length(interval_breaks)])) {
      probs[n_intervals] <- probs[n_intervals] + grid_info$right_tail_mass
    }
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
                                                   sigma_level = 3,
                                                   request = NULL) {
  request <- .as_qc_integration_request(
    request = request,
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  metric <- request$metric
  prior <- request$prior
  cached_state <- request$cached_state
  LSL <- request$LSL
  USL <- request$USL
  target <- request$target
  sigma_level <- request$sigma_level

  interval_breaks <- .integration_validate_interval_breaks(interval_breaks)
  n_intervals <- length(interval_breaks) - 1L
  probs <- .integration_direct_interval_probs(entry, interval_breaks)

  if (is.null(probs)) {
    grid_probs <- .integration_interval_probs_from_grid(entry, interval_breaks)

    probs <- vapply(seq_len(n_intervals), function(j) {
      lo <- interval_breaks[j]
      hi <- interval_breaks[j + 1L]

      tryCatch(
        compute_cpk_prob_integration(
          bounds = c(lo, hi),
          request = request
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
        .integration_quantiles_from_density(
          x,
          density,
          c(h, 1 - h),
          area = area,
          entry = if (is.list(entry)) entry else NULL
        )
      },
      "HPD" = .integration_hdi_from_density(
        x,
        density,
        ci_level,
        area = area,
        entry = if (is.list(entry)) entry else NULL
      ),
      stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
    ))
  }

  stats <- .integration_result_stats(entry)
  if (!is.null(stats)) {
    if (!isTRUE(all.equal(ci_level, 0.95))) {
      stop(
        "Stats-only integration results only support ci_level = 0.95. ",
        "Provide samples or a density grid to request a different interval level."
      )
    }

    return(switch(ci,
      "central" = c(stats["Q2.5"], stats["Q97.5"]),
      "HPD" = c(stats["HDI_Lo"], stats["HDI_Hi"]),
      stop("Unknown ci for integration method. Only 'central' and 'HPD' are supported.")
    ))
  }

  stop("This integration metric distribution does not expose interval data.")
}
