
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
plot_density.bpc <- function(obj, LSL = -1, USL = 1, target = 0, ...) {
  plot_density(.compute_capability_metrics(fit = obj, LSL = LSL, USL = USL, target = target), ...)
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
    ...
  ) {

  what <- match.arg(what, several.ok = TRUE)
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

  layer_line <- ggplot2::geom_line(data = dfLines, mapping = ggplot2::aes(x = x, y = y, group = g, color = g))

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

    layer_points <- jaspGraphs::geom_point(data = dfPoints, mapping = aesPoints, inherit.aes = FALSE)

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
      fill = "grey80",
      alpha = 0.5,
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
      dfCi$y <- 1.15 * max(dfLines$y)
      layer_cibar <- ggplot2::geom_errorbar(
        data = dfCi,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y, group = g, color = g),
        # width = 0.01,
        inherit.aes = FALSE
      )
    }

  }

  layer_text <- NULL
  if (show_ci_text || show_point_text) {

    df_text <- tibble::tibble(
      y = max(dfLines$y) * if (ci != "none" && show_ci_bar) 1.3 else 1.15,
      x = median(dfLines$x),
      g = factor(what)
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
  scale_x <- facet <- NULL
  if (length(what) == 1L) {
    xBreaks <- jaspGraphs::getPrettyAxisBreaks(dfLines$x)
    xLimits <- range(dfLines$x)
    scale_x <- ggplot2::scale_x_continuous(breaks = xBreaks, limits = xLimits)
  } else {
    facet <- ggplot2::facet_grid(cols = ggplot2::vars(g), scales = "free")
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
    facet

  return(plt)
}

#' @export
plot_density.bpc_summary <- function(obj, ...) {
  plot_density(obj = obj$metrics, ...)
}
