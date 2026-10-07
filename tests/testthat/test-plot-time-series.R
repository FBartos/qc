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

test_that("plot_time_series uses pretty y breaks covering data and specifications", {
  plot <- plot_time_series(c(1.2, 2.7, 3.1), LSL = 0.4, USL = 4.6)
  scale <- plot$scales$get_scales("y")

  expect_equal(scale$breaks, pretty(c(0.4, 4.6)))
  expect_equal(scale$limits, range(pretty(c(0.4, 4.6))))
})

test_that("plot_time_series keeps specifications outside the data range visible", {
  plot <- plot_time_series(c(10, 11, 12), LSL = 2, USL = 30)
  limits <- plot$scales$get_scales("y")$limits

  expect_lte(limits[1L], 2)
  expect_gte(limits[2L], 30)
})

test_that("plot_time_series handles constant and empty data", {
  constant <- plot_time_series(c(5, 5, 5))
  limits <- constant$scales$get_scales("y")$limits
  expect_lt(limits[1L], 5)
  expect_gt(limits[2L], 5)
  expect_silent(ggplot2::ggplot_build(constant))

  empty <- plot_time_series(c(NA_real_, NA_real_))
  expect_null(empty$scales$get_scales("y"))
})

test_that("plot_time_series draws filled points and labels the axes", {
  plot <- plot_time_series(c(1, 2, 3))
  points <- ggplot2::ggplot_build(plot)$data[[2L]]

  expect_equal(unique(points$shape), 21)
  expect_equal(unique(points$fill), "grey")
  expect_equal(unique(points$colour), "black")
  expect_equal(unique(points$size), 3)
  expect_equal(plot$labels$x, "Observation number")
  expect_equal(plot$labels$y, "Measurements")
})
