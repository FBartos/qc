test_that("scaled t log density matches the shifted-scale parameterization", {
  x <- c(-1, 0, 1.5)
  df <- 7
  mu <- 1
  sigma <- 2

  expect_equal(
    qc:::lpdf_scaled_t(x, df = df, mu = mu, sigma = sigma),
    stats::dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)
  )
})

test_that("pc frequentist t fits return finite parameters and summaries", {
  set.seed(1)
  x <- 0.5 + 1.2 * stats::rt(200, df = 7)

  fit <- pc(
    x,
    LSL = -6,
    target = 0.5,
    USL = 7,
    distribution = "t",
    bootstrap = FALSE
  )

  expect_s3_class(fit, "pc")
  expect_equal(fit$distribution, "t")
  expect_equal(attr(fit$metrics, "method"), "pc")
  expect_true(all(is.finite(unlist(fit$fit$fit))))
  expect_gt(fit$fit$fit$scale, 0)
  expect_gt(fit$fit$fit$nu, 2)
  expect_lt(abs(fit$fit$fit$mu - mean(x)), 0.5)

  ss <- summary(fit)
  expect_s3_class(ss, "pc_summary")
  expect_true(all(is.finite(ss$summary$mean)))
  expect_true(all(is.na(ss$summary$median)))
})

test_that("pc t bootstrap supports shared density and support-interval contracts", {
  set.seed(2)
  x <- 0.25 + 1.1 * stats::rt(120, df = 9)

  fit <- pc(
    x,
    LSL = -6,
    target = 0.25,
    USL = 7,
    distribution = "t",
    samples = 20,
    seed = 2
  )

  expect_true(all(vapply(fit$metrics_boot, function(draws) all(is.finite(draws)), logical(1))))
  expect_equal(attr(fit$metrics, "method"), "pc")
  expect_equal(attr(fit$metrics_boot, "method"), "pc")

  ss <- summary(fit, interval_probability = c(0.8, 1.0, 1.33, 1.5))
  expect_true(all(is.finite(ss$summary$median)))
  expect_true(all(is.finite(ss$summary$sd)))

  df_density <- extract_density_data(fit, what = c("Cp", "Cpc"))
  expect_s3_class(df_density, "tbl_df")
  expect_equal(sort(unique(as.character(df_density$metric))), c("Cp", "Cpc"))

  support_ci <- qc:::extract_ci_data(
    fit,
    what = "Cpc",
    ci = "support",
    ci_level = 0.95,
    dfDensity = extract_density_data(fit, what = "Cpc"),
    bf_support = list(lower = 0.8, upper = 1.4)
  )
  expect_equal(support_ci$dfCi$xmin, 0.8)
  expect_equal(support_ci$dfCi$xmax, 1.4)

  plot_obj <- plot_density(
    fit,
    what = "Cpc",
    ci = "support",
    bf_support = list(lower = 0.8, upper = 1.4)
  )
  expect_s3_class(plot_obj, "ggplot")
  expect_false(inherits(try(ggplot2::ggplot_build(plot_obj), silent = TRUE), "try-error"))
})

test_that("pc bootstrap parallel path delegates to the parallel backend and closes the cluster", {
  set.seed(3)
  x <- rnorm(30, 10, 2)
  control <- list(samples = 10L, parallel = TRUE, cores = 2L, seed = 123L)
  calls <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    makeCluster = function(cores) {
      calls$cores <- cores
      structure(list(), class = "mock_cluster")
    },
    clusterEvalQ = function(cl, expr) {
      calls$cluster_evalq <- TRUE
      NULL
    },
    clusterExport = function(cl, varlist, envir) {
      calls$exported <- c(calls$exported, list(sort(varlist)))
      NULL
    },
    parLapplyLB = function(cl, X, fun) {
      calls$parallel_work <- length(X)
      lapply(X, fun)
    },
    stopCluster = function(cl) {
      calls$stopped <- TRUE
      NULL
    },
    .package = "parallel"
  )

  out <- qc:::.pc_bootstrap_fit("normal", list(x = x), control)

  expect_equal(nrow(out), 10L)
  expect_equal(calls$cores, 2L)
  expect_equal(calls$parallel_work, 10L)
  expect_true(any(vapply(
    calls$exported,
    identical,
    logical(1),
    c("bootstrap_seeds", "control", "distribution", "x")
  )))
  expect_true(any(vapply(
    calls$exported,
    function(x) ".pc_bootstrap_one" %in% x,
    logical(1)
  )))
  expect_false(isTRUE(calls$cluster_evalq))
  expect_true(isTRUE(calls$stopped))
})

test_that("pc t fitter returns NA fields when optimization fails", {
  testthat::local_mocked_bindings(
    optim = function(...) stop("boom"),
    .package = "stats"
  )

  fit <- qc:::pc_fit_distribution.t(
    structure("t", class = "t"),
    list(x = rnorm(20)),
    control = list()
  )

  expect_named(fit, c("mu", "scale", "nu"))
  expect_true(all(is.na(unlist(fit))))
})
