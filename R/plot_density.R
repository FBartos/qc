
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

  # For integration method, use pre-computed density grids directly
  if (!is.null(obj$method) && obj$method == "integration") {
    plot_density_integration(obj, ...)
  } else {
    # For MCMC method, use pre-computed metrics if spec limits not provided
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
    what = c("Cp", "CpU", "CpL", "Cpk", "Cpc", "Cpm"),
    point_estimate  = c("none", "mean", "median", "mode"),
    ci              = c("none", "central", "HPD", "custom", "support"),
    ci_level        = 0.95,
    ci_custom_left  = NULL,
    ci_custom_right = NULL,
    bf_support      = NULL,
    show_ci_text    = ci != "none",
    show_ci_bar     = ci != "none",
    show_point_text = point_estimate != "none",
    # note sure if these two make sense, we should overwrite them when we do not do
    # facetting, but they could also be vectors I guess
    ci_fill         = "grey60",
    ci_fill_alpha   = 0.8,
    linewidth       = 1,
    single_panel    = FALSE,
    axes            = c("automatic", "fixed", "free", "custom"),
    axes_custom     = list("xmin" = -10, "xmax" = 10, "ymin" = -10, "ymax" = 10),
    priorSummaryObject = NULL, #TODO: implement this!
    ...
  ) {

  what <- match.arg(what, several.ok = TRUE)
  axes <- match.arg(axes)
  point_estimate <- match.arg(point_estimate)
  ci <- match.arg(ci)

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

  listOfDensities <- setNames(lapply(what, function(name) {
    density(obj[[name]])
  }), what)

  dfLines <- vctrs::vec_rbind(!!!lapply(what, function(name) {
    density_estimate <- listOfDensities[[name]]
    tibble::tibble(
      x = density_estimate$x,
      y = density_estimate$y,
      g = factor(name),
    )
  }))

  layer_line <- ggplot2::geom_line(
    data = dfLines,
    mapping = ggplot2::aes(x = x, y = y, group = g, color = g),
    linewidth = linewidth,
  )

  # this is not stricly necessary, but it is pretty convenient.
  listOfFuns <- setNames(lapply(what, function(name) {
    density_estimate <- listOfDensities[[name]]
    stats::approxfun(density_estimate$x, density_estimate$y, rule = 2, yleft = 0, yright = 0)
  }), what)

  if (point_estimate == "none") {
    layer_points <- layer_point_text <- NULL
  } else {

    aesPoints <- ggplot2::aes(x = x, y = y, group = g, color = g, fill = g)

    listOfPointEstimates <- lapply(what, function(name) {

      xValue <- switch(point_estimate,
                       "mean"   = mean(obj[[name]]),
                       "median" = median(obj[[name]]),
                       "mode"   = listOfDensities[[name]]$x[which.max(listOfDensities[[name]]$y)],
                       "none"   = NULL,
                       stop("Unknown point_estimate.")
      )
      yValue <- listOfFuns[[name]](xValue)

      tibble::tibble(
        x = xValue,
        y = yValue,
        g = factor(name)
      )
    })

    dfPoints <- vctrs::vec_rbind(!!!listOfPointEstimates)

    layer_points <- ggplot2::geom_point(data = dfPoints, mapping = aesPoints, inherit.aes = FALSE)

    # TODO: there should only be one geom_text for both the text about the point estimate and the CI bar.
    if (show_point_text) {


    }

  }

  if (ci == "none") {
    layer_area <- layer_cibar <- NULL
  } else {


    listOfCiEstimates <- setNames(lapply(what, function(name) {

      xValue <- switch(ci,
                       "central" = {h <- (1 - ci_level); stats::quantile(obj[[name]], c(h, 1 - h))},
                       "HPD"     = HDInterval::hdi(obj[[name]], ci_level),
                       "custom"  = c(ci_custom_left, ci_custom_right),
                       "support" = stop("Support intervals are not implemented yet."),
                       stop("Unknown ci.")
      )

      xValues <- seq(min(xValue), max(xValue), length.out = 2^8)
      yValues <- c(listOfFuns[[name]](xValues))

      tibble::tibble(
        x = xValue,
        g = factor(name)
      )
    }), what)

    listOfCiEstimatesForArea <- lapply(what, function(name) {

      est <- listOfCiEstimates[[name]]
      xValues <- seq(min(est$x), max(est$x), length.out = 2^8)
      yValues <- c(listOfFuns[[name]](xValues))

      tibble::tibble(
        x = xValues,
        y = yValues,
        g = factor(name)
      )
    })

    dfArea <- vctrs::vec_rbind(!!!listOfCiEstimatesForArea)
    aesArea <- ggplot2::aes(x = x, y = y, group = g)

    layer_area <- ggplot2::geom_area(
      data = dfArea,
      mapping = aesArea,
      color = NA,
      fill = ci_fill,
      alpha = ci_fill_alpha,
      inherit.aes = FALSE,
      stat = "identity",
      position = "identity"
    )

    if (show_ci_bar) {
      dfCi0 <- vctrs::vec_rbind(!!!listOfCiEstimates)
      dfCi  <- tibble::tibble(
        xmin = dfCi0$x[seq(1, nrow(dfCi0), by = 2)],
        xmax = dfCi0$x[seq(2, nrow(dfCi0), by = 2)],
        g    = dfCi0$g[seq(1, nrow(dfCi0), by = 2)]
      )
      dfCi$y <- 1.15 * c(tapply(dfLines$y, dfLines$g, max))
      layer_cibar <- ggplot2::geom_errorbar(
        data = dfCi,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y, group = g, color = g),
        linewidth = linewidth,
        inherit.aes = FALSE
      )
    }

  }

  layer_text <- NULL
  if (show_ci_text || show_point_text) {

    ci_mult <- if (ci != "none" && show_ci_bar) 1.3 else 1.15
    df_text <- dfLines |>
      dplyr::group_by(g) |>
      dplyr::summarize(
        y = max(y) * ci_mult,
        x = median(x),
        .groups = "drop"
      )

    # could be done inside the dplyr::summarize above
    labels <- character(length(df_text$g))
    if (show_point_text) {
      point_estimate_name <- switch(point_estimate,
                                   "mean"   = gettext("Mean"),
                                   "median" = gettext("Median"),
                                   "mode"   = gettext("Mode"),
                                   "none"   = NULL,
                                   stop("Unknown point_estimate.")
      )
      point_estimate_txt <- sprintf("%s = %.3f", point_estimate_name, dfPoints$x)
      labels <- point_estimate_txt
    }

    if (show_ci_text) {

      if (ci == "custom") {
        ci_custom_mass <- listOfDensities <- lapply(what, function(name) {
          mean(obj[[name]] <= ci_custom_right & obj[[name]] >= ci_custom_left)
        })
      }

      ci_txt <- switch(ci,
                       "central" = sprintf("%.1f%% CI [%.3f, %.3f]", 100 * ci_level, dfCi$xmin, dfCi$xmax),
                       "HPD"     = sprintf("%.1f%% CI<sub>HPD</sub> [%.3f, %.3f]", 100 * ci_level, dfCi$xmin, dfCi$xmax),
                       "custom"  = sprintf("P(%.3f &le; &theta; &le; %.3f) =  %.3f", ci_custom_left, ci_custom_right, ci_custom_mass),
                       "support" = sprintf("Support<sub>BF = %.1f</sub>", bf_support),
                       stop("Unknown ci.")
      )

      if (all(labels == ""))
        labels <- ci_txt
      else
        labels <- paste0(labels, "; ", ci_txt)
    }

    df_text$labels <- labels

    layer_text <- ggtext::geom_richtext(
      data = df_text,
      mapping = ggplot2::aes(x = x, y = y, label = labels, group = g),
      fill = NA, label.color = NA, # remove border/ outline
      nudge_y = 0.05 * max(dfLines$y),
    )

  }

  # cannot set xbreaks with free scales!
  # A: yes we can, use
  # ggh4x::facetted_pos_scales
  scale_x <- scale_y <- facet <- NULL
  if (length(what) == 1L || single_panel) {
    xBreaks <- getPrettyAxisBreaks(dfLines$x)
    xLimits <- range(dfLines$x)
    scale_x <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = xLimits)
  } else {
    scales <- switch(axes,
                     "automatic" = "free",
                     "fixed"     = "fixed",
                     "free"      = "free",
                     "custom"    = "fixed",
                     stop("Unknown axes option.")
    )
    if (axes == "custom") {
      if (!is.null(axes_custom[["xmin"]]) && !is.null(axes_custom[["xmax"]])) {
        xbreaks <- getPrettyAxisBreaks(c(axes_custom[["xmin"]], axes_custom[["xmax"]]))
        scale_x <- ggplot2::scale_x_continuous(limits = sort(c(axes_custom[["xmin"]], axes_custom[["xmax"]])))
      }
      if (!is.null(axes_custom[["ymin"]]) && !is.null(axes_custom[["ymax"]])) {
        ybreaks <- getPrettyAxisBreaks(c(axes_custom[["ymin"]], axes_custom[["ymax"]]))
        scale_y <- ggplot2::scale_y_continuous(breaks = ybreaks, limits = sort(c(axes_custom[["ymin"]], axes_custom[["ymax"]])))
      }
    }
    facet <- ggplot2::facet_wrap(~g, scales = scales)
  }

  # TODO: maybe we shouldn't use color?
  plt <- ggplot2::ggplot() +
    layer_area +
    layer_line +
    layer_points +
    layer_cibar +
    layer_text +
    ggplot2::labs(group = "Capability Metric", color = "Capability Metric", fill = "Capability Metric", x = "Value", y = "Density") +
    scale_x +
    scale_y +
    facet

  return(plt)
}

#' @export
plot_density.bpc_summary <- function(obj, ...) {
  plot_density(obj = obj$metrics, ...)
}

#' Plot density for integration method using pre-computed density grids
#'
#' @param obj A bpc object fitted with method = "integration"
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
plot_density_integration <- function(
    obj,
    what = c("Cp", "CpU", "CpL", "Cpk", "Cpc", "Cpm"),
    point_estimate  = c("none", "mean", "median", "mode"),
    ci              = c("none", "central", "HPD"),
    ci_level        = 0.95,
    show_ci_text    = ci != "none",
    show_ci_bar     = ci != "none",
    show_point_text = point_estimate != "none",
    ci_fill         = "grey60",
    ci_fill_alpha   = 0.8,
    linewidth       = 1,
    single_panel    = FALSE,
    axes            = c("automatic", "fixed", "free", "custom"),
    axes_custom     = list("xmin" = -10, "xmax" = 10, "ymin" = -10, "ymax" = 10),
    ...
  ) {

  what <- match.arg(what, several.ok = TRUE)
  axes <- match.arg(axes)
  point_estimate <- match.arg(point_estimate)
  ci <- match.arg(ci)

  BayesTools::check_bool(show_ci_text,    name = "show_ci_text",    check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_ci_bar,     name = "show_ci_bar",     check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(show_point_text, name = "show_point_text", check_length = 1, allow_NA = FALSE)

  # Get pre-computed results from integration
  results <- obj$integration_result$results

  # Build dfLines from pre-computed density grids
  dfLines <- vctrs::vec_rbind(!!!lapply(what, function(name) {
    grid <- results[[name]]$grid
    tibble::tibble(
      x = grid$x,
      y = grid$density,
      g = factor(name),
    )
  }))

  layer_line <- ggplot2::geom_line(
    data = dfLines,
    mapping = ggplot2::aes(x = x, y = y, group = g, color = g),
    linewidth = linewidth,
  )

  # Create approxfun for density evaluation (for point estimates and CI areas)
  listOfFuns <- setNames(lapply(what, function(name) {
    grid <- results[[name]]$grid
    stats::approxfun(grid$x, grid$density, rule = 2, yleft = 0, yright = 0)
  }), what)

  # Get pre-computed stats for each metric
  listOfStats <- setNames(lapply(what, function(name) {
    results[[name]]$stats
  }), what)

  if (point_estimate == "none") {
    layer_points <- layer_point_text <- NULL
  } else {

    aesPoints <- ggplot2::aes(x = x, y = y, group = g, color = g, fill = g)

    listOfPointEstimates <- lapply(what, function(name) {
      stats <- listOfStats[[name]]
      grid <- results[[name]]$grid

      xValue <- switch(point_estimate,
                       "mean"   = stats["Mean"],
                       "median" = stats["Median"],
                       "mode"   = grid$x[which.max(grid$density)],
                       "none"   = NULL,
                       stop("Unknown point_estimate.")
      )
      yValue <- listOfFuns[[name]](xValue)

      tibble::tibble(
        x = unname(xValue),
        y = yValue,
        g = factor(name)
      )
    })

    dfPoints <- vctrs::vec_rbind(!!!listOfPointEstimates)

    layer_points <- ggplot2::geom_point(data = dfPoints, mapping = aesPoints, inherit.aes = FALSE)
  }

  if (ci == "none") {
    layer_area <- layer_cibar <- NULL
  } else {

    # Get CI bounds from pre-computed stats
    listOfCiEstimates <- setNames(lapply(what, function(name) {
      stats <- listOfStats[[name]]

      xValue <- switch(ci,
                       "central" = c(stats["Q2.5"], stats["Q97.5"]),
                       "HPD"     = c(stats["HDI_Lo"], stats["HDI_Hi"]),
                       stop("Unknown ci.")
      )

      # Note: For "central", the stored quantiles are at 2.5% and 97.5%
      # If ci_level != 0.95, we'd need to recompute, but that requires access
      # to the CDF which we don't have readily. For now, warn if ci_level != 0.95
      if (ci == "central" && ci_level != 0.95) {
        warning("For integration method, central CI uses pre-computed 95% quantiles. ",
                "ci_level argument is ignored.")
      }

      tibble::tibble(
        x = unname(xValue),
        g = factor(name)
      )
    }), what)

    listOfCiEstimatesForArea <- lapply(what, function(name) {

      est <- listOfCiEstimates[[name]]
      xValues <- seq(min(est$x), max(est$x), length.out = 2^8)
      yValues <- c(listOfFuns[[name]](xValues))

      tibble::tibble(
        x = xValues,
        y = yValues,
        g = factor(name)
      )
    })

    dfArea <- vctrs::vec_rbind(!!!listOfCiEstimatesForArea)
    aesArea <- ggplot2::aes(x = x, y = y, group = g)

    layer_area <- ggplot2::geom_area(
      data = dfArea,
      mapping = aesArea,
      color = NA,
      fill = ci_fill,
      alpha = ci_fill_alpha,
      inherit.aes = FALSE,
      stat = "identity",
      position = "identity"
    )

    if (show_ci_bar) {
      dfCi0 <- vctrs::vec_rbind(!!!listOfCiEstimates)
      dfCi  <- tibble::tibble(
        xmin = dfCi0$x[seq(1, nrow(dfCi0), by = 2)],
        xmax = dfCi0$x[seq(2, nrow(dfCi0), by = 2)],
        g    = dfCi0$g[seq(1, nrow(dfCi0), by = 2)]
      )
      dfCi$y <- 1.15 * c(tapply(dfLines$y, dfLines$g, max))
      layer_cibar <- ggplot2::geom_errorbar(
        data = dfCi,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y, group = g, color = g),
        linewidth = linewidth,
        inherit.aes = FALSE
      )
    }

  }

  layer_text <- NULL
  if (show_ci_text || show_point_text) {

    ci_mult <- if (ci != "none" && show_ci_bar) 1.3 else 1.15
    df_text <- dfLines |>
      dplyr::group_by(g) |>
      dplyr::summarize(
        y = max(y) * ci_mult,
        x = median(x),
        .groups = "drop"
      )

    labels <- character(length(df_text$g))
    if (show_point_text) {
      point_estimate_name <- switch(point_estimate,
                                   "mean"   = gettext("Mean"),
                                   "median" = gettext("Median"),
                                   "mode"   = gettext("Mode"),
                                   "none"   = NULL,
                                   stop("Unknown point_estimate.")
      )
      point_estimate_txt <- sprintf("%s = %.3f", point_estimate_name, dfPoints$x)
      labels <- point_estimate_txt
    }

    if (show_ci_text) {
      # Note: For integration, we use ci_level = 0.95 (the pre-computed level)
      display_level <- if (ci == "central") 95 else 95

      ci_txt <- switch(ci,
                       "central" = sprintf("%.1f%% CI [%.3f, %.3f]", display_level, dfCi$xmin, dfCi$xmax),
                       "HPD"     = sprintf("%.1f%% CI<sub>HPD</sub> [%.3f, %.3f]", display_level, dfCi$xmin, dfCi$xmax),
                       stop("Unknown ci.")
      )

      if (all(labels == ""))
        labels <- ci_txt
      else
        labels <- paste0(labels, "; ", ci_txt)
    }

    df_text$labels <- labels

    layer_text <- ggtext::geom_richtext(
      data = df_text,
      mapping = ggplot2::aes(x = x, y = y, label = labels, group = g),
      fill = NA, label.color = NA,
      nudge_y = 0.05 * max(dfLines$y),
    )

  }

  scale_x <- scale_y <- facet <- NULL
  if (length(what) == 1L || single_panel) {
    xBreaks <- getPrettyAxisBreaks(dfLines$x)
    xLimits <- range(dfLines$x)
    scale_x <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = xLimits)
  } else {
    scales <- switch(axes,
                     "automatic" = "free",
                     "fixed"     = "fixed",
                     "free"      = "free",
                     "custom"    = "fixed",
                     stop("Unknown axes option.")
    )
    if (axes == "custom") {
      if (!is.null(axes_custom[["xmin"]]) && !is.null(axes_custom[["xmax"]])) {
        xbreaks <- getPrettyAxisBreaks(c(axes_custom[["xmin"]], axes_custom[["xmax"]]))
        scale_x <- ggplot2::scale_x_continuous(breaks = xbreaks, limits = sort(c(axes_custom[["xmin"]], axes_custom[["xmax"]])))
      }
      if (!is.null(axes_custom[["ymin"]]) && !is.null(axes_custom[["ymax"]])) {
        ybreaks <- getPrettyAxisBreaks(c(axes_custom[["ymin"]], axes_custom[["ymax"]]))
        scale_y <- ggplot2::scale_y_continuous(breaks = ybreaks, limits = sort(c(axes_custom[["ymin"]], axes_custom[["ymax"]])))
      }
    }
    facet <- ggplot2::facet_wrap(~g, scales = scales)
  }

  plt <- ggplot2::ggplot() +
    layer_area +
    layer_line +
    layer_points +
    layer_cibar +
    layer_text +
    ggplot2::labs(group = "Capability Metric", color = "Capability Metric", fill = "Capability Metric", x = "Value", y = "Density") +
    scale_x +
    scale_y +
    facet

  return(plt)
}
