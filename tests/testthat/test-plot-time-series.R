test_that("plot_time_series plots observations above specification lines", {
  plot <- plot_time_series(c(1, 2, 3), LSL = 0, target = 2, USL = 4)

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[1L]]$geom, "GeomHline")
  expect_s3_class(plot$layers[[2L]]$geom, "GeomLine")
  expect_s3_class(plot$layers[[3L]]$geom, "GeomPoint")

  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1L]]$yintercept, c(0, 2, 4))
  expect_equal(unique(built$data[[1L]]$linetype), "dashed")
  expect_equal(built$data[[3L]]$x, 1:3)
  expect_equal(built$data[[3L]]$y, c(1, 2, 3))
})

test_that("plot_time_series only adds supplied reference lines", {
  plot <- plot_time_series(c(1, 2, 3), target = 2)

  expect_length(plot$layers, 3L)
  expect_equal(ggplot2::ggplot_build(plot)$data[[1L]]$yintercept, 2)
})
