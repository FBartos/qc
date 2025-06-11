
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

plot_density.bpc <- function(obj, ...) {
  plot_density(compute_capability_metrics(fit = obj), ...)
}

plot_density.capability_metrics <- function(obj, what = c("Cp", "CpU", "CpL", "Cpk", "Cpc", "Cpm"), ...) {

  what <- match.arg(what, several.ok = TRUE)

  nsamples <- length(obj[[1]])
  capability_metrics <- tibble::tibble(
    parameter = rep(what, each = nsamples),
    value     = unlist(obj[what], use.names = FALSE),
  )

  plt <- capability_metrics |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = ggplot2::after_stat(density))) +
    ggplot2::geom_density()

  if (length(what) > 0)
    plt <- plt + ggplot2::facet_grid(cols = ggplot2::vars(parameter), scales = "free")

  return(plt)
}
