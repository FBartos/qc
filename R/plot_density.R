
#' Extract density data for capability metrics
#'
#' Extracts density data as a tibble for visualization or testing purposes.
#' Returns a standardized format regardless of whether the object was fitted
#' with MCMC or integration method.
#'
#' @param obj An object of class `bpc`, `bpc_summary`, or `capability_metrics`
#' @param what Character vector of metrics to extract (default: all six metrics)
#' @param ... Additional arguments (currently unused)
#' @return A tibble with columns:
#'   \describe{
#'     \item{x}{Numeric vector of evaluation points}
#'     \item{density}{Numeric vector of density values}
#'     \item{metric}{Factor indicating the metric name}
#'   }
#' @export
extract_density_data <- function(obj, ...) {
  UseMethod("extract_density_data")
}

#' @rdname extract_density_data
#' @export
extract_density_data.bpc <- function(obj, what = c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"), ...) {
  what <- match.arg(what, several.ok = TRUE)

  if (!is.null(obj$method) && obj$method == "integration") {
    # Integration method: extract from pre-computed grids
    .extract_density_integration(obj$integration_result, what = what)
  } else {
    # MCMC method: compute density from samples
    extract_density_data(obj$metrics, what = what, ...)
  }
}

#' @rdname extract_density_data
#' @export
extract_density_data.bpc_summary <- function(obj, what = c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"), ...) {
  what <- match.arg(what, several.ok = TRUE)

  if (!is.null(obj$integration_result)) {
    # Integration method: use stored results
    .extract_density_integration(obj$integration_result, what = what)
  } else {
    # MCMC method: compute density from samples
    extract_density_data(obj$metrics, what = what, ...)
  }
}

#' @rdname extract_density_data
#' @export
extract_density_data.capability_metrics <- function(obj, what = c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"), ...) {
  what <- match.arg(what, several.ok = TRUE)

  # Check method attribute - for integration metrics, we can't compute density from samples
 if (!is.null(attr(obj, "method")) && attr(obj, "method") == "integration") {
    stop("Cannot extract density from integration capability_metrics directly. ",
         "Use extract_density_data() on the parent bpc or bpc_summary object instead.")
  }

  # MCMC method: compute density from samples
  vctrs::vec_rbind(!!!lapply(what, function(name) {
    .extract_density_from_samples(obj[[name]], metric = name, levels = what)
  }))
}

#' Extract density data from raw samples
#' @param samples Numeric sample vector
#' @param metric Metric name
#' @param levels Factor levels for metric column
#' @return A tibble with columns: x, density, metric
#' @keywords internal
.extract_density_from_samples <- function(samples, metric, levels) {
  samples <- samples[is.finite(samples)]

  if (length(samples) < 2L) {
    stop("Cannot extract density for metric '", metric, "': ",
         "need at least two finite sample values.")
  }

  density_estimate <- stats::density(samples)

  tibble::tibble(
    x = density_estimate$x,
    density = density_estimate$y,
    metric = factor(metric, levels = levels)
  )
}

#' Extract density data from integration results
#' @param integration_result Integration result object containing results list
#' @param what Character vector of metrics to extract
#' @return A tibble with columns: x, density, metric
#' @keywords internal
.extract_density_integration <- function(integration_result, what) {
  results <- integration_result$results

  vctrs::vec_rbind(!!!lapply(what, function(name) {
    result <- results[[name]]
    grid <- result$grid

    if (!is.null(grid)) {
      if (is.null(grid$x) || is.null(grid$density)) {
        stop("Cannot extract integration density for metric '", name, "': ",
             "grid must contain both 'x' and 'density' columns.")
      }

      return(tibble::tibble(
        x = grid$x,
        density = grid$density,
        metric = factor(name, levels = what)
      ))
    }

    if (!is.null(result$samples)) {
      return(.extract_density_from_samples(result$samples, metric = name, levels = what))
    }

    stop("Cannot extract integration density for metric '", name, "': ",
         "result has neither a density grid nor samples.")
  }))
}

#' Extract point estimates for capability metrics
#'
#' @param obj A bpc or capability_metrics object (for MCMC) or NULL (for integration)
#' @param what Character vector of metrics
#' @param point_estimate Type of point estimate ("mean", "median", "mode", "none")
#' @param dfDensity Density data tibble from extract_density_data
#' @param stats_list Optional list of pre-computed stats or integration result entries
#'   (for integration method)
#' @return A tibble with columns: x, y, metric, or NULL if point_estimate is "none"
#' @export
extract_point_estimates <- function(obj, what, point_estimate, dfDensity, stats_list = NULL) {
  if (point_estimate == "none") return(NULL)

  # Build density functions for y-value lookup
  listOfFuns <- setNames(lapply(what, function(name) {
    subset_df <- dfDensity[dfDensity$metric == name, ]
    stats::approxfun(subset_df$x, subset_df$density, rule = 2, yleft = 0, yright = 0)
  }), what)

  vctrs::vec_rbind(!!!lapply(what, function(name) {
    subset_df <- dfDensity[dfDensity$metric == name, ]
    entry <- if (!is.null(stats_list)) stats_list[[name]] else NULL

    xValue <- if (!is.null(stats_list)) {
      # Integration method: use pre-computed stats
      stats <- .integration_result_stats(entry)
      switch(point_estimate,
             "mean"   = stats["Mean"],
             "median" = stats["Median"],
             "mode"   = .integration_result_mode(entry, x = subset_df$x, density = subset_df$density),
             stop("Unknown point_estimate."))
    } else {
      # MCMC method: compute from samples
      switch(point_estimate,
             "mean"   = mean(obj[[name]]),
             "median" = stats::median(obj[[name]]),
             "mode"   = subset_df$x[which.max(subset_df$density)],
             stop("Unknown point_estimate."))
    }

    # Skip infinite point estimates (analytically divergent moments)
    if (!isTRUE(is.finite(xValue)))
      return(tibble::tibble(x = numeric(0), y = numeric(0),
                            metric = factor(character(0), levels = what)))

    tibble::tibble(
      x = unname(xValue),
      y = listOfFuns[[name]](xValue),
      metric = factor(name, levels = what)
    )
  }))
}

#' Extract credible interval data for capability metrics
#'
#' @param obj A bpc or capability_metrics object (for MCMC) or NULL (for integration)
#' @param what Character vector of metrics
#' @param ci Type of CI ("central", "HPD", "custom", "none")
#' @param ci_level CI level (default 0.95)
#' @param dfDensity Density data tibble
#' @param stats_list Optional list of pre-computed stats or integration result entries
#'   (for integration method)
#' @param ci_custom_left Custom CI left bound
#' @param ci_custom_right Custom CI right bound
#' @return A list with dfCi (bounds tibble) and dfArea (filled area tibble), or NULL if ci is "none"
#' @export
extract_ci_data <- function(obj, what, ci, ci_level, dfDensity,
                            stats_list = NULL,
                            ci_custom_left = NULL, ci_custom_right = NULL) {
  if (ci == "none") return(NULL)

  # Build density functions for area computation
  listOfFuns <- setNames(lapply(what, function(name) {
    subset_df <- dfDensity[dfDensity$metric == name, ]
    stats::approxfun(subset_df$x, subset_df$density, rule = 2, yleft = 0, yright = 0)
  }), what)

  # Get CI bounds
  listOfCiEstimates <- setNames(lapply(what, function(name) {
    subset_df <- dfDensity[dfDensity$metric == name, ]

    xValue <- if (!is.null(stats_list)) {
      # Integration method: compute interval from stored result or density grid
      .integration_interval_from_entry(
        stats_list[[name]],
        ci = ci,
        ci_level = ci_level,
        x = subset_df$x,
        density = subset_df$density
      )
    } else {
      # MCMC method: compute from samples
      h <- (1 - ci_level) / 2
      switch(ci,
             "central" = stats::quantile(obj[[name]], c(h, 1 - h)),
             "HPD"     = HDInterval::hdi(obj[[name]], ci_level),
             "custom"  = c(ci_custom_left, ci_custom_right),
             "support" = stop("Support intervals are not implemented yet."),
             stop("Unknown ci."))
    }
    tibble::tibble(x = unname(xValue), metric = factor(name, levels = what))
  }), what)

  # Build dfCi (bounds for error bars)
  dfCi0 <- vctrs::vec_rbind(!!!listOfCiEstimates)
  dfCi <- tibble::tibble(
    xmin = dfCi0$x[seq(1, nrow(dfCi0), by = 2)],
    xmax = dfCi0$x[seq(2, nrow(dfCi0), by = 2)],
    metric = factor(as.character(dfCi0$metric[seq(1, nrow(dfCi0), by = 2)]), levels = what)
  )

  # Build dfArea (filled polygon under curve)
  dfArea <- vctrs::vec_rbind(!!!lapply(what, function(name) {
    est <- listOfCiEstimates[[name]]
    if (length(est$x) < 2L || any(!is.finite(est$x))) {
      return(tibble::tibble(
        x = numeric(0),
        y = numeric(0),
        metric = factor(character(0), levels = what)
      ))
    }
    xValues <- seq(min(est$x), max(est$x), length.out = 256)
    yValues <- listOfFuns[[name]](xValues)
    tibble::tibble(x = xValues, y = yValues, metric = factor(name, levels = what))
  }))

  list(dfCi = dfCi, dfArea = dfArea)
}

#' Default colors for capability region shading
#'
#' Returns the default color palette for the five capability regions.
#' @return A named character vector of colors
#' @export
default_region_colors <- function() {
  c(
    "Incapable"     = "#F87462",
    "Capable"       = "#FAA53D",
    "Satisfactory"  = "#F5CD47",
    "Excellent"     = "#579DFF",
    "Super"         = "#4BCE97"
  )
}

#' Default cutoffs for capability regions
#'
#' Returns the default cutoff values that separate the five capability regions.
#' @return A numeric vector of cutoff values
#' @export
default_region_cutoffs <- function() {
  c(1, 4/3, 1.5, 2)
}

#' Assign region labels based on cutoffs
#'
#' @param x Numeric vector of values
#' @param cutoffs Numeric vector of cutoff values (will be sorted)
#' @param region_names Character vector of region names (length = length(cutoffs) + 1)
#' @return Character vector of region labels
#' @keywords internal
.assign_regions <- function(x, cutoffs, region_names) {
  cutoffs <- sort(cutoffs)
  interval <- findInterval(x, cutoffs, left.open = FALSE) + 1L
  region_names[interval]
}

#' Process density data for region coloring
#'
#' Adds interpolated points at cutoffs to ensure clean color transitions.
#'
#' @param dfLines Data frame with x, y, metric columns
#' @param cutoffs Numeric vector of cutoff values
#' @param region_colors Named vector of colors for each region
#' @return Data frame with region column added and interpolated points at cutoffs
#' @keywords internal
.process_regions <- function(dfLines, cutoffs, region_colors) {
  cutoffs <- sort(cutoffs)
  region_names <- names(region_colors)

  split_df <- split(dfLines, dfLines$metric)

  processed <- lapply(split_df, function(d) {
    y_cut <- stats::approx(d$x, d$y, xout = cutoffs)$y

    cut_df <- d[rep(1L, 2L * length(cutoffs)), ]
    cut_df$x <- rep(cutoffs, each = 2L)
    cut_df$y <- rep(y_cut, each = 2L)

    region_idx <- rep(seq_along(cutoffs), each = 2L)
    region_idx <- region_idx + rep(c(0L, 1L), times = length(cutoffs))
    cut_df$region <- region_names[region_idx]

    d$region <- .assign_regions(d$x, cutoffs, region_names)

    d2 <- rbind(d, cut_df)
    d2[order(d2$x), ]
  })

  result <- do.call(rbind, processed)
  rownames(result) <- NULL
  result$region <- factor(result$region, levels = region_names)
  result
}

#' Build density plot from prepared data
#'
#' Single ggplot skeleton that creates the density plot from prepared dataframes.
#'
#' @param dfLines Density lines tibble with x, y, metric columns
#' @param dfPoints Point estimates tibble (or NULL)
#' @param ci_data List with dfCi and dfArea (or NULL)
#' @param what Character vector of metrics being plotted
#' @param point_estimate Type of point estimate used
#' @param ci Type of CI used
#' @param ci_level CI level
#' @param show_ci_text Whether to show CI text
#' @param show_ci_bar Whether to show CI bar
#' @param show_point_text Whether to show point estimate text
#' @param ci_fill Fill color for CI area
#' @param ci_fill_alpha Alpha for CI fill
#' @param linewidth Line width
#' @param single_panel Whether to use single panel
#' @param axes Axis scaling option
#' @param axes_custom Custom axis limits list
#' @return A ggplot object
#' @keywords internal
build_density_plot <- function(
    dfLines,
    dfPoints = NULL,
    ci_data = NULL,
    what,
    point_estimate = "none",
    ci = "none",
    ci_level = 0.95,
    show_ci_text = FALSE,
    show_ci_bar = FALSE,
    show_point_text = FALSE,
    ci_fill = "grey60",
    ci_fill_alpha = 0.2,
    linewidth = 1,
    single_panel = FALSE,
    axes = "automatic",
    axes_custom = list(),
    textsize = 18,
    colorScheme = NULL,
    stripTextFontsize = NULL,
    show_regions = FALSE,
    region_cutoffs = default_region_cutoffs(),
    region_colors = default_region_colors(),
    region_alpha = 0.55,
    show_cutoff_lines = TRUE
) {

  # Region coloring layers
  layer_region <- layer_cutoffs <- NULL
  if (show_regions) {
    dfRegions <- .process_regions(dfLines, region_cutoffs, region_colors)

    layer_region <- ggplot2::geom_area(
      data = dfRegions,
      mapping = ggplot2::aes(x = x, y = y, fill = region, group = interaction(metric, region)),
      alpha = region_alpha,
      color = NA,
      stat = "identity",
      position = "identity",
      inherit.aes = FALSE
    )

    if (show_cutoff_lines) {
      # If point estimates are shown, limit cutoff lines to mode height to avoid overlapping text
      if (!is.null(dfPoints) && point_estimate != "none") {
        # Use mode height (max density) for each metric
        dfCutoffSegments <- do.call(rbind, lapply(unique(dfLines$metric), function(m) {
          d <- dfLines[dfLines$metric == m, ]
          mode_height <- max(d$y)
          data.frame(
            x = region_cutoffs,
            yend = mode_height,
            metric = m
          )
        }))
        dfCutoffSegments$y <- 0

        layer_cutoffs <- ggplot2::geom_segment(
          data = dfCutoffSegments,
          mapping = ggplot2::aes(x = x, xend = x, y = y, yend = yend),
          linetype = "dashed",
          linewidth = 0.8,
          inherit.aes = FALSE
        )
      } else {
        layer_cutoffs <- ggplot2::geom_vline(
          xintercept = region_cutoffs,
          linetype = "dashed",
          linewidth = 0.8
        )
      }
    }
  }

  has_prior <- any(dfLines$type == "prior")

  # Density line layer
  layer_line <- if (has_prior) {
    ggplot2::geom_line(
      data = dfLines,
      mapping = ggplot2::aes(x = x, y = y, group = interaction(metric, type), color = metric, linetype = type),
      linewidth = linewidth
    )
  } else {
    ggplot2::geom_line(
      data = dfLines,
      mapping = ggplot2::aes(x = x, y = y, group = metric, color = metric),
      linewidth = linewidth
    )
  }

  # Set linetype scale: solid for posterior, dashed for prior
  scale_linetype <- if (has_prior) {
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c("posterior" = "solid", "prior" = "dashed"),
      labels = c("posterior" = "Posterior", "prior" = "Prior")
    )
  } else {
    NULL
  }

  # Point estimate layer
  layer_points <- NULL
  if (!is.null(dfPoints)) {
    layer_points <- ggplot2::geom_point(
      data = dfPoints,
      mapping = ggplot2::aes(x = x, y = y, group = metric, color = metric, fill = metric),
      size = 4,
      inherit.aes = FALSE
    )
  }

  # CI layers
  layer_area <- layer_cibar <- NULL
  dfCi <- NULL
  if (!is.null(ci_data)) {
    dfCi <- ci_data$dfCi
    dfArea <- ci_data$dfArea

    layer_area <- ggplot2::geom_area(
      data = dfArea,
      mapping = ggplot2::aes(x = x, y = y, group = metric, fill = metric),
      color = NA,
      alpha = ci_fill_alpha,
      inherit.aes = FALSE,
      stat = "identity",
      position = "identity"
    )

    if (show_ci_bar) {
      dfCi$y <- 1.15 * c(tapply(dfLines$y, dfLines$metric, max))
      layer_cibar <- ggplot2::geom_errorbar(
        data = dfCi,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y, group = metric, color = metric),
        linewidth = linewidth,
        inherit.aes = FALSE
      )
    }
  }

  # Text layer
  layer_text <- NULL
  if (show_ci_text || show_point_text) {
    ci_mult <- if (ci != "none" && show_ci_bar) 1.3 else 1.15
    x_center <- if (axes == "custom" && !is.null(axes_custom[["xmin"]]) && !is.null(axes_custom[["xmax"]])) {
      mean(c(axes_custom[["xmin"]], axes_custom[["xmax"]]))
    } else {
      NA_real_
    }
    metric_levels <- what[what %in% unique(as.character(dfLines$metric))]
    df_text <- tibble::tibble(
      metric = factor(metric_levels, levels = what),
      y = vapply(metric_levels, function(metric_name) {
        max(dfLines$y[dfLines$metric == metric_name]) * ci_mult
      }, numeric(1)),
      x = vapply(metric_levels, function(metric_name) {
        if (is.na(x_center)) {
          stats::median(dfLines$x[dfLines$metric == metric_name])
        } else {
          x_center
        }
      }, numeric(1))
    )

    point_labels <- rep("", nrow(df_text))
    ci_labels <- rep("", nrow(df_text))

    if (show_point_text && !is.null(dfPoints) && nrow(dfPoints) > 0L) {
      point_estimate_name <- switch(point_estimate,
                                    "mean"   = gettext("Mean"),
                                    "median" = gettext("Median"),
                                    "mode"   = gettext("Mode"),
                                    "")
      point_lookup <- sprintf("%s = %.3f", point_estimate_name, dfPoints$x)
      names(point_lookup) <- as.character(dfPoints$metric)
      point_labels <- unname(point_lookup[metric_levels])
      point_labels[is.na(point_labels)] <- ""
    }

    if (show_ci_text && !is.null(dfCi) && nrow(dfCi) > 0L) {
      ci_lookup <- switch(ci,
                          "central" = sprintf("%.1f%% CI [%.3f, %.3f]", 100 * ci_level, dfCi$xmin, dfCi$xmax),
                          "HPD"     = sprintf("%.1f%% CI<sub>HPD</sub> [%.3f, %.3f]", 100 * ci_level, dfCi$xmin, dfCi$xmax),
                          "custom"  = sprintf("Custom CI [%.3f, %.3f]", dfCi$xmin, dfCi$xmax),
                          "")
      names(ci_lookup) <- as.character(dfCi$metric)
      ci_labels <- unname(ci_lookup[metric_levels])
      ci_labels[is.na(ci_labels)] <- ""
    }

    labels <- point_labels
    has_point <- nzchar(labels)
    has_ci <- nzchar(ci_labels)
    labels[!has_point & has_ci] <- ci_labels[!has_point & has_ci]
    labels[has_point & has_ci] <- paste0(labels[has_point & has_ci], "; ", ci_labels[has_point & has_ci])

    df_text$labels <- labels
    df_text <- df_text[nzchar(df_text$labels), , drop = FALSE]

    if (nrow(df_text) > 0L) {
      layer_text <- ggtext::geom_richtext(
        data = df_text,
        mapping = ggplot2::aes(x = x, y = y, label = labels, group = metric),
        fill = NA, label.color = NA,
        size = textsize,
        nudge_y = 0.05 * max(dfLines$y)
      )
    }
  }

  # Scales and facets
  scale_x <- scale_y <- facet <- NULL
  if (length(what) == 1L || single_panel) {
    if (axes == "custom" && !is.null(axes_custom[["xmin"]]) && !is.null(axes_custom[["xmax"]])) {
      xBreaks <- getPrettyAxisBreaks(c(axes_custom[["xmin"]], axes_custom[["xmax"]]))
      scale_x <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = sort(c(axes_custom[["xmin"]], axes_custom[["xmax"]])))
    } else {
      xBreaks <- getPrettyAxisBreaks(dfLines$x)
      xLimits <- range(dfLines$x)
      scale_x <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = xLimits)
    }
    if (axes == "custom" && !is.null(axes_custom[["ymin"]]) && !is.null(axes_custom[["ymax"]])) {
      ybreaks <- getPrettyAxisBreaks(c(axes_custom[["ymin"]], axes_custom[["ymax"]]))
      scale_y <- ggplot2::scale_y_continuous(breaks = ybreaks, limits = sort(c(axes_custom[["ymin"]], axes_custom[["ymax"]])))
    }
  } else {
    scales <- switch(axes,
                     "automatic" = "free",
                     "fixed"     = "fixed",
                     "free"      = "free",
                     "custom"    = "fixed",
                     "free")
    if (axes == "custom") {
      if (!is.null(axes_custom[["xmin"]]) && !is.null(axes_custom[["xmax"]])) {
        scale_x <- ggplot2::scale_x_continuous(limits = sort(c(axes_custom[["xmin"]], axes_custom[["xmax"]])))
      }
      if (!is.null(axes_custom[["ymin"]]) && !is.null(axes_custom[["ymax"]])) {
        ybreaks <- getPrettyAxisBreaks(c(axes_custom[["ymin"]], axes_custom[["ymax"]]))
        scale_y <- ggplot2::scale_y_continuous(breaks = ybreaks, limits = sort(c(axes_custom[["ymin"]], axes_custom[["ymax"]])))
      }
    }
    facet <- ggplot2::facet_wrap(~metric, scales = scales)
  }

  # Build color palette
  # When colorScheme is NULL or "grey", use a single grey for all metrics
  # When colorScheme is a vector of colors, use those directly
  nColors <- length(unique(dfLines$metric))
  if (is.null(colorScheme) || identical(colorScheme, "grey")) {
    plotColors <- rep("grey50", nColors)
  } else if (is.character(colorScheme) && length(colorScheme) >= nColors) {
    # colorScheme is a vector of colors passed from caller
    plotColors <- colorScheme[seq_len(nColors)]
  } else {
    # Fallback: single grey
    plotColors <- rep("grey50", nColors)
  }

  # Strip text theme
  stripTheme <- NULL
  if (!is.null(stripTextFontsize)) {
    stripTheme <- ggplot2::theme(strip.text = ggplot2::element_text(size = stripTextFontsize))
  }

  # Build final plot
  plt <- ggplot2::ggplot() +
    layer_region +
    layer_cutoffs +
    layer_line +
    scale_linetype +
    layer_points +
    layer_area +
    layer_cibar +
    layer_text +
    ggplot2::scale_color_manual(values = plotColors)

  # Use region colors for fill if regions are shown, otherwise use metric colors
  if (show_regions) {
    plt <- plt + ggplot2::scale_fill_manual(
      values = region_colors,
      breaks = names(region_colors)
    )
  } else {
    plt <- plt + ggplot2::scale_fill_manual(values = plotColors)
  }

  plt <- plt +
    ggplot2::labs(
      group = "Capability Metric",
      color = "Capability Metric",
      fill = if (show_regions) NULL else "Capability Metric",
      x = "Value",
      y = "Density"
    ) +
    scale_x +
    scale_y +
    facet +
    stripTheme

  return(plt)
}

#' Plot density for the posterior distribution of one or more capability metrics
#'
#' @param obj An object of class `bpc`, `bpc_capability_metrics`, or `bpc_summary`.
#' @param LSL Lower Specification Limit
#' @param target Target value
#' @param USL Upper Specification Limit
#' @param ...
#'
#' @export
#'
plot_density <- function(obj, ...) {
  UseMethod("plot_density")
}

#' @export
plot_density.bpc <- function(obj, LSL = NULL, USL = NULL, target = NULL, ...) {

  # Validate that either all or none of LSL, USL, target are provided
  provided <- c(!is.null(LSL), !is.null(USL), !is.null(target))
  if (any(provided) && !all(provided)) {
    stop("If any of LSL, USL, or target are provided, all three must be specified.")
  }

  # For integration method
  if (!is.null(obj$method) && obj$method == "integration") {
    if (all(provided)) {
      new_result <- .bpc_fit_integration(
        data = numeric(0),
        LSL = LSL,
        USL = USL,
        target = target,
        prior_mu = obj$prior_mu %||% "Jeffreys_mu",
        prior_sigma = obj$prior_sigma %||% "Jeffreys_sigma",
        sigma = obj$sigma %||% obj$integration_result$sigma %||% 3,
        cached_state = obj$integration_result$cached_state
      )
      plot_density_integration_results(new_result, ...)
    } else {
      # Use pre-computed density grids
      plot_density_integration(obj, ...)
    }
  } else {
    # For MCMC method
    if (is.null(LSL) && is.null(USL) && is.null(target)) {
      # Use already computed metrics from the fit
      plot_density(obj$metrics, ...)
    } else {
      # Recompute metrics with new spec limits if provided
      plot_density(.compute_capability_metrics(fit = obj, LSL = LSL, USL = USL, target = target), ...)
    }
  }
}

#' @export
plot_density.capability_metrics <- function(
    obj,
    what = c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"),
    point_estimate  = c("none", "mean", "median", "mode"),
    ci              = c("none", "central", "HPD", "custom", "support"),
    ci_level        = 0.95,
    ci_custom_left  = NULL,
    ci_custom_right = NULL,
    bf_support      = NULL,
    show_ci_text    = ci != "none",
    show_ci_bar     = ci != "none",
    show_point_text = point_estimate != "none",
    ci_fill         = "grey60",
    ci_fill_alpha   = 0.7,
    linewidth       = 1,
    single_panel    = FALSE,
    axes            = c("automatic", "fixed", "free", "custom"),
    axes_custom     = list("xmin" = -10, "xmax" = 10, "ymin" = -10, "ymax" = 10),
    textsize        = 18,
    priorSummaryObject = NULL,
    colorScheme     = NULL,
    stripTextFontsize = NULL,
    show_regions    = FALSE,
    region_cutoffs  = default_region_cutoffs(),
    region_colors   = default_region_colors(),
    region_alpha    = 0.55,
    show_cutoff_lines = TRUE,
    ...
  ) {

  what <- match.arg(what, several.ok = TRUE)
  axes <- match.arg(axes)
  point_estimate <- match.arg(point_estimate)
  ci <- match.arg(ci)

  # Input validation
  if (ci == "central" || ci == "HPD") {
    BayesTools::check_real(ci_level, name = "ci_level", check_length = 1, lower = 0, upper = 1, allow_NA = FALSE)
  } else if (ci == "custom") {
    BayesTools::check_real(ci_custom_left,  name = "ci_custom_left",  check_length = 1, lower = -Inf,           upper = ci_custom_right, allow_NA = FALSE)
    BayesTools::check_real(ci_custom_right, name = "ci_custom_right", check_length = 1, lower = ci_custom_left, upper = Inf,             allow_NA = FALSE)
  } else if (ci == "support") {
    BayesTools::check_list(bf_support, name = "bf_support", check_names = c("lower", "upper"), allow_NULL = TRUE)
  }
  BayesTools::check_bool(show_ci_text,    name = "show_ci_text",    check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_ci_bar,     name = "show_ci_bar",     check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_point_text, name = "show_point_text", check_length = 1, allow_NA = FALSE)

  # Extract density data
  dfDensity <- extract_density_data(obj, what = what)
  dfDensity$type <- "posterior"

  # Extract prior density if available and not showing regions
  if (!is.null(priorSummaryObject) && !show_regions) {
    dfDensityPrior <- extract_density_data(priorSummaryObject, what = what)
    dfDensityPrior$type <- "prior"
    dfDensity <- vctrs::vec_rbind(dfDensity, dfDensityPrior)
  }

  # Extract point estimates (MCMC: pass obj for sample-based computation)
  dfPoints <- extract_point_estimates(
    obj = obj,
    what = what,
    point_estimate = point_estimate,
    dfDensity = dfDensity[dfDensity$type == "posterior", ],
    stats_list = NULL  # MCMC uses samples
  )

  # Extract CI data (MCMC: pass obj for sample-based computation)
  ci_data <- extract_ci_data(
    obj = obj,
    what = what,
    ci = ci,
    ci_level = ci_level,
    dfDensity = dfDensity[dfDensity$type == "posterior", ],
    stats_list = NULL,  # MCMC uses samples
    ci_custom_left = ci_custom_left,
    ci_custom_right = ci_custom_right
  )

  # Prepare dfLines for plotting (rename columns: density -> y, metric -> g)
  dfLines <- tibble::tibble(
    x = dfDensity$x,
    y = dfDensity$density,
    metric = dfDensity$metric,
    type = dfDensity$type
  )

  # Build plot using single skeleton
  build_density_plot(
    dfLines = dfLines,
    dfPoints = dfPoints,
    ci_data = ci_data,
    what = what,
    point_estimate = point_estimate,
    ci = ci,
    ci_level = ci_level,
    show_ci_text = show_ci_text,
    show_ci_bar = show_ci_bar,
    show_point_text = show_point_text,
    ci_fill = ci_fill,
    ci_fill_alpha = ci_fill_alpha,
    linewidth = linewidth,
    single_panel = single_panel,
    axes = axes,
    axes_custom = axes_custom,
    textsize = textsize,
    colorScheme = colorScheme,
    stripTextFontsize = stripTextFontsize,
    show_regions = show_regions,
    region_cutoffs = region_cutoffs,
    region_colors = region_colors,
    region_alpha = region_alpha,
    show_cutoff_lines = show_cutoff_lines
  )
}

#' @export
plot_density.bpc_summary <- function(obj, ..., priorSummaryObject = NULL) {
  # Check if this is from integration method
  if (!is.null(obj$integration_result)) {
    # Use the stored integration results for plotting
    plot_density_integration_results(obj$integration_result, ..., priorSummaryObject = priorSummaryObject)
  } else {
    # MCMC method: metrics contain samples
    plot_density(obj = obj$metrics, ..., priorSummaryObject = priorSummaryObject)
  }
}

#' Plot density for integration method using pre-computed density grids
#'
#' @param obj A bpc object fitted with method = "integration"
#' @param ... Additional arguments passed to plot_density_integration_results
#' @return A ggplot object
#' @keywords internal
plot_density_integration <- function(obj, ...) {
  # Get pre-computed results from integration and delegate to results plotter
  plot_density_integration_results(obj$integration_result, ...)
}

#' Plot density from integration results
#'
#' @param integration_result Integration result object containing results list with grid and stats
#' @param what Character vector of metrics to plot
#' @param point_estimate Type of point estimate to show
#' @param ci Type of credible interval to show
#' @param ci_level Credible interval level (for central/HPD)
#' @param show_ci_text Whether to show CI text
#' @param show_ci_bar Whether to show CI bar
#' @param show_point_text Whether to show point estimate text
#' @param ci_fill Fill color for CI area
#' @param ci_fill_alpha Alpha for CI fill
#' @param linewidth Line width
#' @param single_panel If TRUE, plot all metrics in a single panel
#' @param axes Axis scaling option
#' @param axes_custom Custom axis limits
#' @param ... Additional arguments (ignored)
#' @return A ggplot object
#' @keywords internal
plot_density_integration_results <- function(
    integration_result,
    what = c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"),
    point_estimate  = c("none", "mean", "median", "mode"),
    ci              = c("none", "central", "HPD"),
    ci_level        = 0.95,
    show_ci_text    = ci != "none",
    show_ci_bar     = ci != "none",
    show_point_text = point_estimate != "none",
    ci_fill         = "grey60",
    ci_fill_alpha   = 0.7,
    linewidth       = 1,
    single_panel    = FALSE,
    axes            = c("automatic", "fixed", "free", "custom"),
    axes_custom     = list("xmin" = -10, "xmax" = 10, "ymin" = -10, "ymax" = 10),
    textsize        = 18,
    priorSummaryObject = NULL,
    colorScheme     = NULL,
    stripTextFontsize = NULL,
    show_regions    = FALSE,
    region_cutoffs  = default_region_cutoffs(),
    region_colors   = default_region_colors(),
    region_alpha    = 0.55,
    show_cutoff_lines = TRUE,
    ...
  ) {

  what <- match.arg(what, several.ok = TRUE)
  axes <- match.arg(axes)
  point_estimate <- match.arg(point_estimate)
  ci <- match.arg(ci)

  BayesTools::check_bool(show_ci_text,    name = "show_ci_text",    check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_ci_bar,     name = "show_ci_bar",     check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_point_text, name = "show_point_text", check_length = 1, allow_NA = FALSE)

  # Extract density data using shared extractor
  dfDensity <- .extract_density_integration(integration_result, what = what)
  dfDensity$type <- "posterior"

  # Extract prior density if available and not showing regions
  if (!is.null(priorSummaryObject) && !show_regions) {
    dfDensityPrior <- extract_density_data(priorSummaryObject, what = what)
    dfDensityPrior$type <- "prior"
    dfDensity <- vctrs::vec_rbind(dfDensity, dfDensityPrior)
  }

  # Build stats list for integration method
  results <- integration_result$results
  stats_list <- setNames(lapply(what, function(name) {
    results[[name]]
  }), what)

  # Extract point estimates (integration: pass stats_list for pre-computed values)
  dfPoints <- extract_point_estimates(
    obj = NULL,  # Not needed for integration
    what = what,
    point_estimate = point_estimate,
    dfDensity = dfDensity[dfDensity$type == "posterior", ],
    stats_list = stats_list
  )

  # Extract CI data (integration: pass stats_list for pre-computed values)
  ci_data <- extract_ci_data(
    obj = NULL,  # Not needed for integration
    what = what,
    ci = ci,
    ci_level = ci_level,
    dfDensity = dfDensity[dfDensity$type == "posterior", ],
    stats_list = stats_list
  )

  # Prepare dfLines for plotting
  dfLines <- tibble::tibble(
    x = dfDensity$x,
    y = dfDensity$density,
    metric = dfDensity$metric,
    type = dfDensity$type
  )

  # Build plot using single skeleton
  build_density_plot(
    dfLines = dfLines,
    dfPoints = dfPoints,
    ci_data = ci_data,
    what = what,
    point_estimate = point_estimate,
    ci = ci,
    ci_level = ci_level,
    show_ci_text = show_ci_text,
    show_ci_bar = show_ci_bar,
    show_point_text = show_point_text,
    ci_fill = ci_fill,
    ci_fill_alpha = ci_fill_alpha,
    linewidth = linewidth,
    single_panel = single_panel,
    axes = axes,
    axes_custom = axes_custom,
    textsize = textsize,
    colorScheme = colorScheme,
    stripTextFontsize = stripTextFontsize,
    show_regions = show_regions,
    region_cutoffs = region_cutoffs,
    region_colors = region_colors,
    region_alpha = region_alpha,
    show_cutoff_lines = show_cutoff_lines
  )
}
