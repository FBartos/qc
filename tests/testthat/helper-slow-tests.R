qc_run_slow_tests <- function() {
  value <- tolower(trimws(Sys.getenv("QC_RUN_SLOW_TESTS", "false")))
  value %in% c("1", "true", "yes", "on")
}

skip_if_not_slow_tests <- function() {
  if (!qc_run_slow_tests()) {
    testthat::skip(
      "Slow regression test skipped. Set QC_RUN_SLOW_TESTS=true to run it."
    )
  }
}
