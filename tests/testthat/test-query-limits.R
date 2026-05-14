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

test_that("summary recomputes cached metrics when sigma changes without new limits", {
  set.seed(123)
  x <- rnorm(30, 10, 1)

  fit_bpc <- bpc(
    x,
    LSL = 7,
    target = 10,
    USL = 13,
    method = "integration",
    prior = "Jeffreys",
    sigma = 3
  )
  bpc_sigma3 <- summary(fit_bpc, sigma = 3)
  bpc_sigma4 <- summary(fit_bpc, sigma = 4)
  bpc_sigma4_explicit <- summary(fit_bpc, LSL = 7, target = 10, USL = 13, sigma = 4)

  expect_equal(bpc_sigma4$summary$mean, bpc_sigma4_explicit$summary$mean, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(bpc_sigma3$summary$mean, bpc_sigma4$summary$mean)))

  fit_pc <- pc(
    x,
    LSL = 7,
    target = 10,
    USL = 13,
    bootstrap = FALSE,
    sigma = 3
  )
  pc_sigma3 <- summary(fit_pc, sigma = 3)
  pc_sigma4 <- summary(fit_pc, sigma = 4)
  pc_sigma4_explicit <- summary(fit_pc, LSL = 7, target = 10, USL = 13, sigma = 4)

  expect_equal(pc_sigma4$summary$mean, pc_sigma4_explicit$summary$mean)
  expect_false(isTRUE(all.equal(pc_sigma3$summary$mean, pc_sigma4$summary$mean)))
})
