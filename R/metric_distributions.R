# Internal adapter layer for backend-neutral capability-metric distributions.

.qc_metric_names <- function() {
  c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
}

.new_qc_metric_context <- function(prior = NULL,
                                   cached_state = NULL,
                                   LSL = NULL,
                                   USL = NULL,
                                   target = NULL,
                                   sigma_level = NULL,
                                   distribution = NULL,
                                   divergence = NULL) {
  list(
    prior = prior,
    cached_state = cached_state,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    distribution = distribution,
    divergence = divergence
  )
}

.new_qc_metric_distribution <- function(metric,
                                        samples = NULL,
                                        entry = NULL,
                                        stats = NULL,
                                        context = NULL) {
  if (is.null(samples) && is.null(entry) && is.null(stats)) {
    stop("A metric distribution requires samples, an integration entry, or stats.")
  }

  if (is.null(stats) && !is.null(entry) && is.list(entry) && !is.null(entry$stats)) {
    stats <- entry$stats
  }

  structure(
    list(
      metric = metric,
      samples = samples,
      entry = entry,
      stats = stats,
      context = context
    ),
    class = "qc_metric_distribution"
  )
}

.qc_match_metrics <- function(what) {
  match.arg(what, choices = .qc_metric_names(), several.ok = TRUE)
}

.integration_context_from_metrics <- function(obj) {
  .new_qc_metric_context(
    prior = attr(obj, "prior"),
    cached_state = attr(obj, "cached_state"),
    LSL = attr(obj, "LSL"),
    USL = attr(obj, "USL"),
    target = attr(obj, "target"),
    sigma_level = attr(obj, "sigma"),
    distribution = attr(obj, "distribution"),
    divergence = attr(obj, "divergence")
  )
}

.entry_distributions_from_results <- function(results, what, context = NULL) {
  setNames(lapply(what, function(metric_name) {
    .new_qc_metric_distribution(
      metric = metric_name,
      entry = results[[metric_name]],
      context = context
    )
  }), what)
}

.sample_distributions_from_metrics <- function(metrics, what) {
  setNames(lapply(what, function(metric_name) {
    .new_qc_metric_distribution(
      metric = metric_name,
      samples = metrics[[metric_name]]
    )
  }), what)
}

.as_qc_metric_distributions <- function(obj, what = .qc_metric_names(), ...) {
  UseMethod(".as_qc_metric_distributions")
}

.as_qc_metric_distributions.bpc <- function(obj, what = .qc_metric_names(), ...) {
  .as_qc_metric_distributions(obj$metrics, what = what, ...)
}

.as_qc_metric_distributions.bpc_summary <- function(obj, what = .qc_metric_names(), ...) {
  .as_qc_metric_distributions(obj$metrics, what = what, ...)
}

.pc_distribution_metrics <- function(obj,
                                     context = "Distribution queries for `pc` objects") {
  metrics <- obj$metrics_boot

  if (is.null(metrics)) {
    stop(
      context,
      " require bootstrap draws. Refit with `bootstrap = TRUE`."
    )
  }

  metrics
}

.as_qc_metric_distributions.pc <- function(obj, what = .qc_metric_names(), ...) {
  .as_qc_metric_distributions(.pc_distribution_metrics(obj), what = what, ...)
}

.as_qc_metric_distributions.pc_summary <- function(obj, what = .qc_metric_names(), ...) {
  .as_qc_metric_distributions(.pc_distribution_metrics(obj), what = what, ...)
}

.as_qc_metric_distributions.qc_integration_result <- function(obj, what = .qc_metric_names(), ...) {
  .as_qc_metric_distributions(obj$metrics, what = what, ...)
}

.as_qc_metric_distributions.capability_metrics <- function(obj, what = .qc_metric_names(), ...) {
  what <- .qc_match_metrics(what)
  distributions <- attr(obj, "distributions")

  if (!is.null(distributions)) {
    return(distributions[what])
  }

  if (identical(attr(obj, "method"), "integration")) {
    results <- attr(obj, "results")
    if (is.null(results)) {
      stop(
        "Integration-backed capability_metrics require attached integration results. ",
        "Use the parent bpc/bpc_summary object or refit the model."
      )
    }

    return(.entry_distributions_from_results(
      results,
      what = what,
      context = .integration_context_from_metrics(obj)
    ))
  }

  .sample_distributions_from_metrics(obj, what = what)
}

.coerce_qc_metric_distributions <- function(obj = NULL,
                                            what = .qc_metric_names(),
                                            stats_list = NULL,
                                            context = NULL) {
  what <- .qc_match_metrics(what)

  if (!is.null(stats_list)) {
    return(setNames(lapply(what, function(metric_name) {
      item <- stats_list[[metric_name]]
      if (inherits(item, "qc_metric_distribution")) {
        return(item)
      }

      if (is.null(item)) {
        stop("Missing distribution input for metric '", metric_name, "'.")
      }

      if (is.numeric(item) && !is.null(names(item))) {
        return(.new_qc_metric_distribution(
          metric = metric_name,
          stats = item,
          context = context
        ))
      }

      .new_qc_metric_distribution(
        metric = metric_name,
        entry = item,
        context = context
      )
    }), what))
  }

  if (is.null(obj)) {
    stop("Provide either `obj` or `stats_list`.")
  }

  .as_qc_metric_distributions(obj, what = what)
}

.density_frame_from_samples <- function(samples) {
  samples <- samples[is.finite(samples)]

  if (length(samples) < 2L) {
    stop("Cannot estimate a density: need at least two finite sample values.")
  }

  density_estimate <- stats::density(samples)
  tibble::tibble(
    x = density_estimate$x,
    density = density_estimate$y
  )
}

.density_frame_from_entry <- function(entry) {
  if (is.list(entry) && !is.null(entry$grid)) {
    grid <- entry$grid

    if (is.null(grid$x) || is.null(grid$density)) {
      stop("Integration grids must contain both 'x' and 'density' columns.")
    }

    return(tibble::tibble(
      x = grid$x,
      density = grid$density
    ))
  }

  if (is.list(entry) && !is.null(entry$samples)) {
    return(.density_frame_from_samples(entry$samples))
  }

  stop("This integration metric distribution does not expose density data.")
}

.qc_metric_distribution_density_data <- function(dist, levels = dist$metric) {
  df <- if (!is.null(dist$samples)) {
    .density_frame_from_samples(dist$samples)
  } else {
    .density_frame_from_entry(dist$entry)
  }

  tibble::tibble(
    x = df$x,
    density = df$density,
    metric = factor(dist$metric, levels = levels)
  )
}

.sample_metric_interval <- function(samples, ci, ci_level,
                                     ci_custom_left = NULL,
                                     ci_custom_right = NULL,
                                     bf_support = NULL) {
  h <- (1 - ci_level) / 2

  switch(ci,
    "central" = unname(stats::quantile(samples, c(h, 1 - h), na.rm = TRUE)),
    "HPD" = unname(HDInterval::hdi(samples, ci_level)),
    "custom" = c(ci_custom_left, ci_custom_right),
    "support" = c(bf_support$lower, bf_support$upper),
    stop("Unknown ci.")
  )
}

.sample_metric_stats <- function(samples, ci_level = 0.95) {
  central <- .sample_metric_interval(samples, "central", ci_level)
  hdi <- .sample_metric_interval(samples, "HPD", ci_level)

  c(
    Mean = mean(samples),
    Median = stats::median(samples),
    SD = stats::sd(samples),
    Q2.5 = central[1],
    Q97.5 = central[2],
    HDI_Lo = hdi[1],
    HDI_Hi = hdi[2]
  )
}

.qc_metric_distribution_stats <- function(dist, ci_level = 0.95) {
  if (!is.null(dist$stats)) {
    return(dist$stats)
  }

  if (!is.null(dist$samples)) {
    stats <- .sample_metric_stats(dist$samples, ci_level = ci_level)
    dist$stats <- stats
    return(stats)
  }

  .integration_result_stats(dist$entry)
}

.qc_metric_distribution_point_estimate <- function(dist,
                                                   point_estimate,
                                                   density_frame = NULL) {
  if (point_estimate == "none") {
    return(NULL)
  }

  if (is.null(dist$samples)) {
    stats <- .qc_metric_distribution_stats(dist)
    return(switch(point_estimate,
      "mean" = unname(stats["Mean"]),
      "median" = unname(stats["Median"]),
      "mode" = .integration_result_mode(
        dist$entry %||% dist$stats,
        x = if (is.null(density_frame)) NULL else density_frame$x,
        density = if (is.null(density_frame)) NULL else density_frame$density
      ),
      stop("Unknown point_estimate.")
    ))
  }

  switch(point_estimate,
    "mean" = mean(dist$samples),
    "median" = stats::median(dist$samples),
    "mode" = {
      if (is.null(density_frame)) {
        density_frame <- .density_frame_from_samples(dist$samples)
      }
      density_frame$x[which.max(density_frame$density)]
    },
    stop("Unknown point_estimate.")
  )
}

.qc_metric_distribution_interval <- function(dist,
                                             ci,
                                             ci_level,
                                             density_frame = NULL,
                                             ci_custom_left = NULL,
                                             ci_custom_right = NULL,
                                             bf_support = NULL) {
  if (!is.null(dist$samples)) {
    return(.sample_metric_interval(
      dist$samples,
      ci = ci,
      ci_level = ci_level,
      ci_custom_left = ci_custom_left,
      ci_custom_right = ci_custom_right,
      bf_support = bf_support
    ))
  }

  .integration_interval_from_entry(
    dist$entry %||% dist$stats,
    ci = ci,
    ci_level = ci_level,
    x = if (is.null(density_frame)) NULL else density_frame$x,
    density = if (is.null(density_frame)) NULL else density_frame$density
  )
}

.qc_metric_distribution_interval_probs <- function(dist, interval_breaks) {
  if (!is.null(dist$samples)) {
    return(.integration_interval_probs_from_samples(dist$samples, interval_breaks))
  }

  entry <- dist$entry
  direct_probs <- .integration_direct_interval_probs(entry, interval_breaks)
  if (!is.null(direct_probs)) {
    return(direct_probs)
  }

  context <- dist$context
  has_solver_context <- !is.null(context) &&
    !is.null(context$prior) &&
    !is.null(context$LSL) &&
    !is.null(context$USL)

  if (!has_solver_context) {
    grid_probs <- .integration_interval_probs_from_grid(entry, interval_breaks)
    if (!is.null(grid_probs)) {
      return(grid_probs)
    }

    stop("This integration metric distribution does not expose interval probabilities.")
  }

  request <- .new_qc_integration_request(
    data = numeric(0),
    LSL = context$LSL,
    USL = context$USL,
    prior = context$prior,
    metric = dist$metric,
    target = context$target,
    cached_state = context$cached_state,
    sigma_level = context$sigma_level
  )

  .integration_interval_probs_from_entry(
    entry = entry,
    interval_breaks = interval_breaks,
    metric = dist$metric,
    prior = context$prior,
    cached_state = context$cached_state,
    LSL = context$LSL,
    USL = context$USL,
    target = context$target,
    sigma_level = context$sigma_level,
    request = request
  )
}

.qc_metric_distribution_divergence <- function(dist) {
  if (is.list(dist$entry) && !is.null(dist$entry$divergence_info)) {
    return(dist$entry$divergence_info)
  }

  context <- dist$context
  if (!is.null(context) && !is.null(context$divergence)) {
    return(context$divergence[[dist$metric]])
  }

  NULL
}

.qc_metric_distributions_summary_table <- function(distributions, ci_level) {
  rows <- lapply(distributions, function(dist) {
    stats <- .qc_metric_distribution_stats(dist, ci_level = ci_level)
    interval <- .qc_metric_distribution_interval(
      dist,
      ci = "central",
      ci_level = ci_level
    )

    data.frame(
      metric = dist$metric,
      mean = unname(stats["Mean"]),
      median = unname(stats["Median"]),
      sd = unname(stats["SD"]),
      lower = interval[1],
      upper = interval[2],
      row.names = NULL
    )
  })

  tibble::as_tibble(do.call(rbind, rows))
}

.qc_metric_distributions_interval_summary <- function(distributions, interval_probability) {
  interval_breaks <- c(-Inf, interval_probability, Inf)

  rows <- lapply(distributions, function(dist) {
    probs <- .qc_metric_distribution_interval_probs(dist, interval_breaks)
    df <- as.data.frame(t(probs))
    names(df) <- levels(cut(0, breaks = interval_breaks, include.lowest = TRUE))
    df$metric <- dist$metric
    df[, c("metric", names(df)[names(df) != "metric"])]
  })

  tibble::as_tibble(do.call(rbind, rows))
}

.qc_metric_distributions_divergence_table <- function(distributions) {
  rows <- lapply(distributions, function(dist) {
    info <- .qc_metric_distribution_divergence(dist)
    if (is.null(info) || (!info$mean_divergent && !info$sd_divergent)) {
      return(NULL)
    }

    data.frame(
      metric = dist$metric,
      mean_divergent = info$mean_divergent,
      sd_divergent = info$sd_divergent,
      alpha = info$alpha,
      reason = if (is.null(info$reason)) "" else info$reason,
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(NULL)
  }

  tibble::as_tibble(do.call(rbind, rows))
}

.qc_metric_distributions_summary_bundle <- function(distributions,
                                                    ci_level,
                                                    interval_probability,
                                                    mean_override = NULL,
                                                    include_divergence = TRUE) {
  summary <- .qc_metric_distributions_summary_table(
    distributions,
    ci_level = ci_level
  )

  if (!is.null(mean_override)) {
    summary$mean <- unname(unlist(mean_override[summary$metric]))
  }

  list(
    summary = summary,
    interval_summary = .qc_metric_distributions_interval_summary(
      distributions,
      interval_probability = interval_probability
    ),
    divergence_diagnostics = if (isTRUE(include_divergence)) {
      .qc_metric_distributions_divergence_table(distributions)
    } else {
      NULL
    }
  )
}
