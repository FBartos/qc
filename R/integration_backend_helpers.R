.integration_resolve_target <- function(metric, target, LSL, USL) {
  if (!is.null(target) || !metric %in% c("Cpm", "Cpc")) {
    return(target)
  }

  (LSL + USL) / 2
}

.integration_can_use_density_solver <- function(prior, metric, enabled = TRUE) {
  density_metrics <- c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")
  isTRUE(enabled) &&
    metric %in% density_metrics &&
    .integration_can_use_density(prior, metric)
}

.integration_density_support_lower <- function(metric) {
  if (.metric_can_be_negative(metric)) -Inf else 0
}

.integration_request_args <- function(request) {
  request <- .as_qc_integration_request(request = request)

  list(
    data = request$data,
    LSL = request$LSL,
    USL = request$USL,
    prior = request$prior,
    metric = request$metric,
    target = request$target,
    cached_state = request$cached_state,
    sigma_level = request$sigma_level
  )
}

.integration_make_solver <- function(data = NULL, LSL = NULL, USL = NULL,
                                     prior = NULL, metric = "Cpk",
                                     target = NULL, cached_state = NULL,
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

  do.call(make_solver, .integration_request_args(request))
}

.integration_make_density_solver <- function(data = NULL, LSL = NULL, USL = NULL,
                                             prior = NULL,
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

  do.call(make_density_solver, .integration_request_args(request))
}

.integration_compute_metric_moments <- function(use_analytic = TRUE,
                                                request) {
  request <- .as_qc_integration_request(request = request)

  do.call(
    compute_metric_moments,
    c(
      .integration_request_args(request),
      list(use_analytic = use_analytic)
    )
  )
}

.integration_interval_prob_from_survival <- function(S, bounds) {
  lower_bound <- min(bounds)
  upper_bound <- max(bounds)
  p_lower <- if (is.finite(lower_bound)) S(lower_bound) else 1
  p_upper <- if (is.finite(upper_bound)) S(upper_bound) else 0
  min(max(p_lower - p_upper, 0), 1)
}

.integration_interval_probs_from_samples <- function(samples, interval_breaks) {
  bin_counts <- table(cut(
    samples,
    breaks = interval_breaks,
    include.lowest = TRUE
  ))
  as.numeric(bin_counts) / sum(bin_counts)
}

.integration_interval_probs_from_degenerate <- function(dist, interval_breaks) {
  n_intervals <- length(interval_breaks) - 1L

  if (identical(dist$type, "point")) {
    probs <- numeric(n_intervals)
    point_idx <- .integration_interval_index(dist$value, interval_breaks)
    if (!is.na(point_idx)) {
      probs[point_idx] <- 1
    }
    return(probs)
  }

  probs <- vapply(seq_len(n_intervals), function(i) {
    .degenerate_metric_prob(dist, interval_breaks[c(i, i + 1L)])
  }, numeric(1))
  total <- sum(probs)
  if (total > 0) {
    probs <- probs / total
  }
  probs
}

.integration_basic_interval_probs_from_entry <- function(entry, interval_breaks) {
  if (!is.list(entry)) {
    return(NULL)
  }

  if (!is.null(entry$degenerate)) {
    return(.integration_interval_probs_from_degenerate(entry$degenerate, interval_breaks))
  }

  if (!is.null(entry$samples)) {
    return(.integration_interval_probs_from_samples(entry$samples, interval_breaks))
  }

  NULL
}

.integration_direct_interval_probs <- function(entry, interval_breaks) {
  interval_breaks <- .integration_validate_interval_breaks(interval_breaks)

  direct_probs <- .integration_basic_interval_probs_from_entry(entry, interval_breaks)
  if (!is.null(direct_probs)) {
    return(direct_probs)
  }

  .integration_interval_probs_from_grid(entry, interval_breaks)
}

.new_integration_grid_entry <- function(metric, grid, area, tail_info,
                                        stats = NULL,
                                        divergence_info = NULL) {
  structure(
    list(
      metric = metric,
      grid = grid,
      area = area,
      grid_mass = max(min(area, 1), 0),
      mass_nonpositive = tail_info$mass_nonpositive,
      left_tail_mass = tail_info$left_tail_mass,
      right_tail_mass = tail_info$right_tail_mass,
      support_lower = tail_info$support_lower,
      stats = stats,
      divergence_info = divergence_info
    ),
    class = c("qc_integration_entry", "list")
  )
}

.integration_metric_case <- function(data = NULL, LSL = NULL, USL = NULL,
                                     prior = NULL, metric = "Cpk",
                                     target = NULL, cached_state = NULL,
                                     sigma_level = 3,
                                     context = c("analysis", "probability"),
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

  context <- match.arg(context)
  metric_can_be_negative <- .metric_can_be_negative(request$metric)

  degenerate <- NULL
  if (inherits(request$prior, "PriorConjugate")) {
    posterior_info <- .compute_validated_conjugate_posterior(
      request$prior,
      request$data,
      request$cached_state,
      context = switch(
        context,
        "analysis" = "The conjugate posterior for density analysis",
        "probability" = "The conjugate posterior for interval probability computation"
      )
    )
    if (posterior_info$is_degenerate) {
      degenerate <- .degenerate_conjugate_metric_distribution(
        posterior_info$post$mu_n,
        posterior_info$post$k_n,
        request$LSL,
        request$USL,
        request$target,
        request$metric,
        sigma_level = request$sigma_level
      )
    }
  }

  list(
    target = request$target,
    metric_can_be_negative = metric_can_be_negative,
    degenerate = degenerate,
    request = request
  )
}

.integration_backend_resolver <- function(data = NULL, LSL = NULL, USL = NULL,
                                          prior = NULL, metric = "Cpk",
                                          target = NULL, cached_state = NULL,
                                          sigma_level = 3,
                                          prefer_density = TRUE,
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

  can_use_density <- .integration_can_use_density_solver(
    prior = request$prior,
    metric = request$metric,
    enabled = prefer_density
  )

  density_error <- NULL
  pdf_fn <- NULL
  if (can_use_density) {
    pdf_fn <- tryCatch(
      .integration_make_density_solver(request = request),
      error = function(e) {
        density_error <<- e
        NULL
      }
    )
    can_use_density <- !is.null(pdf_fn)
  }

  S <- NULL
  if (!can_use_density || !prefer_density) {
    S <- tryCatch(
      .integration_make_solver(request = request),
      error = function(e) NULL
    )
  }

  list(
    mode = if (can_use_density) "density" else "survival",
    pdf_fn = pdf_fn,
    S = S,
    can_use_density = can_use_density,
    density_error = density_error
  )
}

.integration_grid_entry <- function(metric, grid_df, area, backend,
                                    metric_can_be_negative,
                                    stats = NULL,
                                    divergence_info = NULL) {
  grid_x <- c(
    grid_df$x_left[1],
    grid_df$x_right
  )
  grid_mass <- max(min(area, 1), 0)

  tail_info <- if (identical(backend$mode, "density")) {
    .integration_density_tail_masses(
      pdf_fn = backend$pdf_fn,
      grid_x = grid_x,
      grid_mass = grid_mass,
      support_lower = 0
    )
  } else {
    .integration_survival_tail_masses(
      S = backend$S,
      grid_x = grid_x,
      grid_mass = grid_mass,
      support_lower = .integration_density_support_lower(metric),
      metric_can_be_negative = metric_can_be_negative
    )
  }

  .new_integration_grid_entry(
    metric = metric,
    grid = grid_df,
    area = grid_mass,
    tail_info = tail_info,
    stats = stats,
    divergence_info = divergence_info
  )
}
