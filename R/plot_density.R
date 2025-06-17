
#' Plot density for the posterior distribution of one or more capability metrics
#'
#' @param obj
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
plot_density <- function(obj, ...) {
  UseMethod("plot_density")
}

#' @export
plot_density.bpc <- function(obj, LSL = -1, USL = 1, target = 0, ...) {
  plot_density(.bpc_compute_capability_metrics(fit = obj, LSL = LSL, USL = USL, target = target), ...)
}

#' @export
plot_density.bpc_capability_metrics <- function(obj, what = c("Cp", "CpU", "CpL", "Cpk", "Cpc", "Cpm"), ...) {

  what <- match.arg(what, several.ok = TRUE)

  nsamples <- length(obj[[1]])
  capability_metrics <- tibble::tibble(
    parameter = rep(what, each = nsamples),
    value     = unlist(obj[what], use.names = FALSE),
  )

  plt <- capability_metrics |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = ggplot2::after_stat(density))) +
    ggplot2::geom_density()

  if (length(what) > 1L)
    plt <- plt + ggplot2::facet_grid(cols = ggplot2::vars(parameter), scales = "free")

  return(plt)
}

#' @export
plot_density.bpc_summary <- function(obj, what = c("Cp", "CpU", "CpL", "Cpk", "Cpc", "Cpm"), ...) {
  plot_density(obj = obj$metrics, what = what, ...)
}
