#' Plot observations over time
#'
#' Creates a time series plot of raw observations against their observation
#' number. Specification limits and the target, when supplied, are drawn as
#' dashed horizontal lines.
#'
#' @param x Numeric vector of observations.
#' @param LSL Optional lower specification limit.
#' @param target Optional target value.
#' @param USL Optional upper specification limit.
#' @return A ggplot object.
#' @export
plot_time_series <- function(x, LSL = NULL, target = NULL, USL = NULL) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }

  specification_arguments <- list(LSL = LSL, target = target, USL = USL)
  valid_specifications <- vapply(specification_arguments, function(value) {
    is.null(value) || (is.numeric(value) && length(value) == 1L && is.finite(value))
  }, logical(1))
  if (!all(valid_specifications)) {
    stop("Specification limits and target must be finite numeric values.", call. = FALSE)
  }
  specification_values <- unlist(specification_arguments, use.names = FALSE)

  observations <- data.frame(observation = seq_along(x), value = x)
  specifications <- data.frame(value = unname(specification_values))

  yvalues <- c(x, specification_values)
  yvalues <- yvalues[is.finite(yvalues)]
  yscale <- NULL
  if (length(yvalues) > 0L) {
    yrange <- range(yvalues)
    if (yrange[1L] == yrange[2L]) {
      yrange <- yrange + c(-0.5, 0.5) * max(abs(yrange[1L]), 1)
    }
    ybreaks <- getPrettyAxisBreaks(yrange)
    yscale <- ggplot2::scale_y_continuous(limits = range(ybreaks), breaks = ybreaks)
  }

  plot <- ggplot2::ggplot(observations, ggplot2::aes(x = observation, y = value))
  if (nrow(specifications) > 0L) {
    plot <- plot + ggplot2::geom_hline(
      data = specifications,
      ggplot2::aes(yintercept = value),
      linetype = "dashed"
    )
  }

  plot +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3, colour = "black", fill = "grey", shape = 21) +
    yscale +
    ggplot2::labs(x = "Observation number", y = "Measurements")
}
