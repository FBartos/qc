test_that("summary and plot_density require all specification limits together", {
  fit <- bpc(rnorm(30, 10, 2), LSL = 4, target = 10, USL = 16, method = "integration")

  expect_error(
    summary(fit, LSL = 2, target = 10),
    regexp = "all three must be specified"
  )

  expect_error(
    plot_density(fit, LSL = 2, target = 10),
    regexp = "all three must be specified"
  )
})

test_that("summary.pc requires all specification limits together", {
  fit <- pc(rnorm(30, 10, 2), LSL = 4, target = 10, USL = 16, samples = 20, seed = 1)

  expect_error(
    summary(fit, LSL = 2, target = 10),
    regexp = "all three must be specified"
  )
})
