# run_mr() and run_coloc() accept a lazy dplyr table (here DuckDB via
# dbplyr) for `exposure` and `outcome` as well as a data frame (issue #42).
# Every test runs the same analysis on both input types and requires the
# results to be identical, not merely close.

skip_if_no_duckdb <- function() {
  testthat::skip_if_not_installed("DBI")
  testthat::skip_if_not_installed("dbplyr")
  testthat::skip_if_not_installed("duckdb")
}

# Write data frames into an in-memory DuckDB, one table per named argument,
# and return lazy tbl()s on that single connection plus the connection.
# The caller closes it: `on.exit(close_duckdb(db))`.
as_duckdb <- function(...) {
  con <- DBI::dbConnect(duckdb::duckdb())
  frames <- list(...)
  tables <- lapply(names(frames), function(nm) {
    DBI::dbWriteTable(con, nm, frames[[nm]])
    dplyr::tbl(con, nm)
  })
  c(stats::setNames(tables, names(frames)), list(con = con))
}

close_duckdb <- function(db) {
  DBI::dbDisconnect(db$con, shutdown = TRUE)
}

# Run `fun(exposure, outcome, ...)` with PLINK's console output captured,
# messages suppressed and warnings collected. The seed is fixed because
# TwoSampleMR::format_data() draws a random id.outcome; the previous
# allele-check record is cleared so a stale one can never pass for a match.
run_captured <- function(fun, exposure, outcome, ...) {
  .mrpipeline_env$last_allele_check <- NULL
  set.seed(42)
  warnings <- character()
  result <- NULL
  invisible(utils::capture.output(
    result <- withCallingHandlers(
      suppressMessages(fun(exposure = exposure, outcome = outcome, ...)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  ))
  list(
    result = result,
    warnings = warnings,
    allele_check = .mrpipeline_env$last_allele_check
  )
}

# Everything but the wall-clock timings must match exactly.
expect_same_run <- function(lazy, eager) {
  strip <- function(run) {
    run$result$timing <- NULL
    run
  }
  testthat::expect_identical(strip(lazy), strip(eager))
}

run_mr_both <- function(...) {
  db <- as_duckdb(exposure = cd40_exposure, outcome = sjogren_outcome)
  on.exit(close_duckdb(db))
  args <- list(
    exposure_id = "CD40",
    outcome_id = "SjD",
    verbose = FALSE,
    ...
  )
  eager <- do.call(
    run_captured,
    c(list(run_mr, cd40_exposure, sjogren_outcome), args)
  )
  lazy <- do.call(run_captured, c(list(run_mr, db$exposure, db$outcome), args))
  list(eager = eager, lazy = lazy)
}

# --- run_mr ------------------------------------------------------------------

test_that("cis-MR on a DuckDB table matches the data-frame run exactly", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  runs <- run_mr_both(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    ld_correct = TRUE,
    methods = c(
      "ivw_random",
      "ivw_fixed",
      "egger",
      "weighted_median",
      "heterogeneity",
      "pleiotropy",
      "loo"
    )
  )
  expect_equal(runs$eager$result$status, "success")
  expect_gt(nrow(runs$eager$result$instruments), 1)
  expect_false(is.null(runs$eager$allele_check))
  expect_same_run(runs$lazy, runs$eager)
})

test_that("manual-instrument MR on a DuckDB table matches the data-frame run", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")

  runs <- run_mr_both(
    instruments = cd40_three_snps(),
    methods = c("ivw_random", "ivw_fixed", "heterogeneity")
  )
  expect_equal(runs$eager$result$status, "success")
  expect_equal(nrow(runs$eager$result$instruments), 3L)
  expect_same_run(runs$lazy, runs$eager)
})

test_that("genome-wide MR on a DuckDB table matches the data-frame run", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  runs <- run_mr_both(
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    methods = c("ivw_random", "egger", "weighted_median")
  )
  expect_equal(runs$eager$result$status, "success")
  expect_same_run(runs$lazy, runs$eager)
})

test_that("early returns are the same for a DuckDB table", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")

  # No SNPs in the cis window: the filter collects zero rows.
  runs <- run_mr_both(
    instrument_region = list(chromosome = "1", start = 1e6, end = 2e6)
  )
  expect_equal(runs$eager$result$status, "no_instruments")
  expect_same_run(runs$lazy, runs$eager)

  # Manual instruments absent from the exposure.
  runs <- run_mr_both(instruments = c("rs_absent_1", "rs_absent_2"))
  expect_equal(runs$eager$result$status, "no_instruments")
  expect_same_run(runs$lazy, runs$eager)
})

# --- run_coloc ---------------------------------------------------------------

test_that("coloc on a DuckDB table matches the data-frame run exactly", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")
  skip_if_not_installed("coloc")
  bfile <- ld_bfile()

  methods <- if (rlang::is_installed("susieR")) {
    c("abf", "susie", "signals")
  } else {
    "abf"
  }
  db <- as_duckdb(exposure = cd40_exposure, outcome = sjogren_outcome)
  on.exit(close_duckdb(db))
  args <- list(
    exposure_id = "CD40",
    outcome_id = "SjD",
    gene_chr = 20,
    gene_start = cd40_region$start,
    gene_end = cd40_region$end,
    outcome_type = "cc",
    outcome_s = 6098 / 41420,
    bfile = bfile,
    methods = methods,
    verbose = FALSE
  )
  eager <- do.call(
    run_captured,
    c(list(run_coloc, cd40_exposure, sjogren_outcome), args)
  )
  lazy <- do.call(
    run_captured,
    c(list(run_coloc, db$exposure, db$outcome), args)
  )
  expect_equal(eager$result$status, "success")
  expect_gt(eager$result$n_snps, 0)
  expect_same_run(lazy, eager)

  # No SNPs in the window: both filters must collect before nrow().
  args$gene_chr <- 1
  eager <- do.call(
    run_captured,
    c(list(run_coloc, cd40_exposure, sjogren_outcome), args)
  )
  lazy <- do.call(
    run_captured,
    c(list(run_coloc, db$exposure, db$outcome), args)
  )
  expect_equal(eager$result$status, "no_snps_in_region")
  expect_same_run(lazy, eager)
})

# --- check_allele_orientation_gwas -------------------------------------------

test_that("the sampled allele check is identical for every input pairing", {
  skip_if_no_duckdb()
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_gwas_fixture()

  # n_sample = 10 of 37 shared non-instrument SNPs, so the radix-sorted
  # evenly spaced sample is drawn rather than every shared SNP taken.
  record <- function(exposure, outcome) {
    .mrpipeline_env$last_allele_check <- NULL
    set.seed(42)
    suppressWarnings(check_allele_orientation_gwas(
      exposure,
      outcome,
      f$instruments,
      allele_check = "none",
      n_sample = 10L
    ))
  }
  eager <- record(f$exposure, f$bug)
  expect_equal(eager$n_sampled, 10L)
  expect_equal(eager$status, "fail")

  # Same connection: the intersect runs in DuckDB.
  db <- as_duckdb(exposure = f$exposure, outcome = f$bug)
  on.exit(close_duckdb(db), add = TRUE)
  expect_true(same_backend(db$exposure, db$outcome))
  expect_identical(record(db$exposure, db$outcome), eager)

  # Mixed and cross-connection: the rsID columns are collected first.
  other <- as_duckdb(outcome = f$bug)
  on.exit(close_duckdb(other), add = TRUE)
  expect_false(same_backend(db$exposure, other$outcome))
  expect_identical(record(f$exposure, db$outcome), eager)
  expect_identical(record(db$exposure, f$bug), eager)
  expect_identical(record(db$exposure, other$outcome), eager)
})

test_that("shared_rsids intersects in the database on one connection", {
  skip_if_no_duckdb()
  exposure <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4", NA, "rs2"))
  outcome <- data.frame(rsids = c("rs4", "rs2", "rs3", NA, "rs9"))
  db <- as_duckdb(exposure = exposure, outcome = outcome)
  on.exit(close_duckdb(db))

  eager <- shared_rsids(exposure, outcome, instrument_snps = "rs3")
  expect_setequal(eager, c("rs2", "rs4"))
  expect_setequal(
    shared_rsids(db$exposure, db$outcome, instrument_snps = "rs3"),
    eager
  )
  # No instruments to exclude is not an error on either path.
  expect_setequal(
    shared_rsids(db$exposure, db$outcome, instrument_snps = character()),
    c("rs2", "rs3", "rs4")
  )

  # Two data frames are silent; a lazy input collected for want of a shared
  # connection says so, and only when verbose.
  expect_no_message(shared_rsids(exposure, outcome, "rs3", verbose = TRUE))
  expect_no_message(shared_rsids(
    db$exposure,
    db$outcome,
    "rs3",
    verbose = TRUE
  ))
  expect_message(
    shared_rsids(exposure, db$outcome, "rs3", verbose = TRUE),
    "one database"
  )
  expect_no_message(shared_rsids(exposure, db$outcome, "rs3", verbose = FALSE))
})
