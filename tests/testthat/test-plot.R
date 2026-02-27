# Tests for plot methods

skip_if_not_installed("ggplot2")

# --- mr_result plots ---------------------------------------------------------

test_that("plot.mr_result returns NULL for non-success status", {
  res <- new_mr_result(
    status = "no_instruments",
    status_reason = "No instruments found"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

test_that("plot.mr_result returns NULL for empty results", {
  res <- new_mr_result(
    results = data.frame(
      exposure = character(),
      outcome = character(),
      method = character(),
      nsnp = integer(),
      b = numeric(),
      se = numeric(),
      pval = numeric()
    ),
    status = "success"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

# --- coloc_result plots ------------------------------------------------------

test_that("plot.coloc_result returns NULL for non-success status", {
  res <- new_coloc_result(
    status = "no_snps_in_region",
    status_reason = "No SNPs in the region"
  )

  expect_message(p <- plot(res), "Cannot plot")
  expect_null(p)
})

test_that("plot.coloc_result pp_bar returns ggplot with ABF results", {
  abf_summary <- c(
    nsnps = 50,
    PP.H0.abf = 0.01,
    PP.H1.abf = 0.05,
    PP.H2.abf = 0.02,
    PP.H3.abf = 0.12,
    PP.H4.abf = 0.80
  )

  res <- new_coloc_result(
    coloc_abf = list(summary = abf_summary),
    n_snps = 50L,
    status = "success"
  )

  p <- plot(res, type = "pp_bar")
  expect_s3_class(p, "ggplot")
})

test_that("plot.coloc_result pp_bar returns NULL without ABF", {
  res <- new_coloc_result(
    coloc_abf = NULL,
    n_snps = 50L,
    status = "success"
  )

  expect_message(p <- plot(res, type = "pp_bar"), "No ABF results")
  expect_null(p)
})

test_that("plot.coloc_result regional returns ggplot with harmonised data", {
  dat <- data.frame(
    pos.exposure = 1:10 * 1000,
    pval.exposure = runif(10, 1e-8, 0.1),
    pos.outcome = 1:10 * 1000,
    pval.outcome = runif(10, 1e-8, 0.1)
  )

  res <- new_coloc_result(
    harmonised_data = dat,
    n_snps = 10L,
    status = "success"
  )

  p <- plot(res, type = "regional")
  expect_s3_class(p, "ggplot")
})

test_that("plot.coloc_result regional returns NULL without data", {
  res <- new_coloc_result(
    harmonised_data = data.frame(),
    n_snps = 0L,
    status = "success"
  )

  expect_message(p <- plot(res, type = "regional"), "No harmonised data")
  expect_null(p)
})
