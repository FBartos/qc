test_that("ci_fill controls the integration CI ribbon fill", {
  set.seed(1)
  fit <- bpc(rnorm(40, 0, 1), LSL = -3, USL = 3, target = 0, method = "integration")

  p <- plot_density(fit, what = "Cp", ci = "central", ci_fill = "red")
  area_layers <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomArea"), logical(1)))

  expect_length(area_layers, 1)

  built <- ggplot2::ggplot_build(p)
  expect_equal(unique(built$data[[area_layers]]$fill), "red")
})

test_that("region cutoffs outside the density support keep finite boundary heights", {
  df_lines <- data.frame(
    x = c(3, 4, 5),
    y = c(0.2, 1, 0.3),
    metric = factor("Cp")
  )
  region_colors <- c("Needs Improvement" = "#000000", "Capable" = "#111111", "Excellent" = "#222222")

  processed <- qc:::.process_regions(df_lines, cutoffs = c(1, 7), region_colors = region_colors)

  expect_false(anyNA(processed$y))
  expect_equal(unique(processed$y[processed$x == 1]), 0.2)
  expect_equal(unique(processed$y[processed$x == 7]), 0.3)
})
