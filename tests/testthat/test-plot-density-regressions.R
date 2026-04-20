test_that("ci_fill controls the integration CI ribbon fill", {
  set.seed(1)
  fit <- bpc(rnorm(40, 0, 1), LSL = -3, USL = 3, target = 0, method = "integration")

  p <- plot_density(fit, what = "Cp", ci = "central", ci_fill = "red")
  area_layers <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomArea"), logical(1)))

  expect_length(area_layers, 1)

  built <- ggplot2::ggplot_build(p)
  expect_equal(unique(built$data[[area_layers]]$fill), "red")
})
