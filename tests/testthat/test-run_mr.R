# Tests for run_mr()

# --- Argument validation ----------------------------------------------------

test_that("run_mr errors when ld_correct = TRUE and bfile = NULL", {
  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      ld_correct = TRUE,
      bfile = NULL
    ),
    "bfile"
  )
})

test_that("run_mr validates methods argument", {
  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      methods = "not_a_method"
    ),
    "Unknown method"
  )
})

test_that("run_mr validates exclude_regions argument", {
  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      exclude_regions = "not_a_df"
    ),
    "data frame"
  )

  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      exclude_regions = data.frame(chr = "6", start = 34e6, end = 26e6)
    ),
    "start.*<=.*end"
  )

  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      exclude_regions = data.frame(chr = "6")
    ),
    "missing"
  )
})

# --- Instrument selection: manual mode --------------------------------------

test_that("run_mr returns mr_result with no_instruments when manual instruments not found", {
  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-8,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  suppressMessages(
    expect_warning(
      expect_warning(
        result <- run_mr(
          exposure = exposure,
          exposure_id = "test_exp",
          outcome = data.frame(),
          outcome_id = "test_out",
          instruments = c("rs999"),
          instruments_strict = FALSE
        ),
        "not found"
      ),
      "No manual instruments"
    )
  )
  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "no_instruments")
  expect_true(nchar(result$status_reason) > 0)
})

test_that("run_mr errors with instruments_strict = TRUE for missing instruments", {
  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-8,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  expect_error(
    run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = data.frame(),
      outcome_id = "test_out",
      instruments = c("rs999"),
      instruments_strict = TRUE
    ),
    "not found"
  )
})

# --- Instrument selection: cis-MR mode --------------------------------------

test_that("run_mr returns mr_result with no_instruments when no instruments in cis region", {
  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 0.5,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = 100,
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = data.frame(),
      outcome_id = "test_out",
      instrument_region = list(chromosome = 1, start = 1, end = 200),
      pval_thresh = 5e-8
    ),
    "No significant instruments"
  )
  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "no_instruments")
  expect_true(nchar(result$status_reason) > 0)
})

# --- Instrument selection: genome-wide mode ---------------------------------

test_that("run_mr returns mr_result with no_instruments when no genome-wide significant instruments", {
  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 0.5,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = data.frame(),
      outcome_id = "test_out",
      pval_thresh = 5e-8
    ),
    "No genome-wide significant"
  )
  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "no_instruments")
  expect_true(nchar(result$status_reason) > 0)
})

# --- Region exclusion -------------------------------------------------------

test_that("run_mr removes instruments in excluded regions", {
  skip_if_not_installed("TwoSampleMR")

  # All instruments in MHC-like region
  exposure <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.1, 0.2),
    se.exposure = c(0.05, 0.05),
    effect_allele.exposure = c("A", "G"),
    other_allele.exposure = c("G", "T"),
    pval.exposure = c(1e-10, 1e-10),
    eaf.exposure = c(0.3, 0.4),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "6",
    pos.exposure = c(28e6, 30e6),
    stringsAsFactors = FALSE
  )

  mhc_region <- data.frame(chr = "6", start = 26e6, end = 34e6)

  expect_warning(
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = data.frame(),
      outcome_id = "test_out",
      instruments = c("rs1", "rs2"),
      exclude_regions = mhc_region
    ),
    "excluded regions"
  )
  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "no_instruments")
})

# --- mr_result S3 class -----------------------------------------------------

test_that("mr_result print works with empty results", {
  res <- new_mr_result()
  expect_message(print(res), "no results")
})

test_that("mr_result print works with results", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "exp1",
      outcome = "out1",
      method = "IVW",
      nsnp = 5,
      b = 0.1,
      se = 0.05,
      pval = 0.01,
      stringsAsFactors = FALSE
    ),
    f_stats = list(per_snp = rep(30, 5), mean = 30, min = 25)
  )
  expect_message(print(res), "exp1")
})

test_that("mr_result print shows status for failed results", {
  res <- new_mr_result(
    status = "no_instruments",
    status_reason = "No significant instruments in cis region for 'PCSK9'"
  )
  expect_message(print(res), "no_instruments")
  expect_message(print(res), "PCSK9")
})

test_that("mr_result summary works", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "exp1",
      outcome = "out1",
      method = c("IVW", "Egger"),
      nsnp = c(5, 5),
      b = c(0.1, 0.12),
      se = c(0.05, 0.06),
      pval = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    f_stats = list(per_snp = rep(30, 5), mean = 30, min = 25),
    methods_skipped = c(presso = "Requires >= 3 instruments"),
    params = list(exposure_id = "exp1", outcome_id = "out1")
  )
  expect_message(summary(res), "MR Results")
})

test_that("mr_result summary shows status for failed results", {
  res <- new_mr_result(
    status = "no_instruments",
    status_reason = "All instruments removed",
    params = list(exposure_id = "exp1", outcome_id = "out1")
  )
  expect_message(summary(res), "no_instruments")
})

# --- Single instrument Wald ratio (with mock) --------------------------------

test_that("run_mr returns Wald ratio for single instrument", {
  skip_if_not_installed("TwoSampleMR")

  # Build a minimal exposure with 1 SNP
  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.5,
    se.exposure = 0.1,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-10,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = 1000,
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  # Build outcome with matching SNP
  outcome <- data.frame(
    rsids = "rs1",
    beta = 0.1,
    se = 0.05,
    pval = 0.01,
    eaf = 0.3,
    effect_allele = "A",
    other_allele = "G",
    chr = "1",
    pos = 1000,
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = "rs1",
      methods = c("ivw", "egger", "weighted_median")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "success")
  expect_equal(nrow(result$results), 1)
  expect_equal(result$results$method, "Wald ratio")

  # IVW, egger, weighted_median should be skipped
  expect_true(all(
    c("ivw", "egger", "weighted_median") %in%
      names(result$methods_skipped)
  ))

  # F-stat should be computed
  expect_equal(length(result$f_stats$per_snp), 1)
  expect_true(result$f_stats$mean > 0)
})

# --- Method skipping for 2 SNPs -------------------------------------------

test_that("run_mr skips egger/weighted_median/presso with 2 instruments", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.5, 0.3),
    se.exposure = c(0.1, 0.1),
    effect_allele.exposure = c("A", "G"),
    other_allele.exposure = c("G", "T"),
    pval.exposure = c(1e-10, 1e-8),
    eaf.exposure = c(0.3, 0.4),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = c("1", "1"),
    pos.exposure = c(1000, 2000),
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = c("rs1", "rs2"),
    beta = c(0.1, 0.05),
    se = c(0.05, 0.03),
    pval = c(0.01, 0.1),
    eaf = c(0.3, 0.4),
    effect_allele = c("A", "G"),
    other_allele = c("G", "T"),
    chr = c("1", "1"),
    pos = c(1000, 2000),
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = c("rs1", "rs2"),
      methods = c("ivw", "egger", "weighted_median", "presso")
    )
  })

  expect_s3_class(result, "mr_result")

  # IVW should work
  expect_true(any(stringr::str_detect(
    result$results$method,
    "variance weighted"
  )))

  # Egger, weighted_median, presso should be skipped
  expect_true("egger" %in% names(result$methods_skipped))
  expect_true("weighted_median" %in% names(result$methods_skipped))
  expect_true("presso" %in% names(result$methods_skipped))
})

# --- Heterogeneity and leave-one-out ---------------------------------------

test_that("run_mr computes heterogeneity with 2 instruments but skips loo", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.5, 0.3),
    se.exposure = c(0.1, 0.1),
    effect_allele.exposure = c("A", "G"),
    other_allele.exposure = c("G", "T"),
    pval.exposure = c(1e-10, 1e-8),
    eaf.exposure = c(0.3, 0.4),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = c("1", "1"),
    pos.exposure = c(1000, 2000),
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = c("rs1", "rs2"),
    beta = c(0.1, 0.05),
    se = c(0.05, 0.03),
    pval = c(0.01, 0.1),
    eaf = c(0.3, 0.4),
    effect_allele = c("A", "G"),
    other_allele = c("G", "T"),
    chr = c("1", "1"),
    pos = c(1000, 2000),
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = c("rs1", "rs2"),
      methods = c("ivw", "heterogeneity", "loo")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_false(is.null(result$heterogeneity))
  expect_true(all(c("Q", "Q_df", "Q_pval") %in% names(result$heterogeneity)))

  expect_null(result$loo)
  expect_true("loo" %in% names(result$methods_skipped))
  expect_match(result$methods_skipped[["loo"]], "Requires >= 3")
})

test_that("run_mr computes leave-one-out with 3 instruments", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = c("rs1", "rs2", "rs3"),
    beta.exposure = c(0.5, 0.3, 0.4),
    se.exposure = c(0.1, 0.1, 0.1),
    effect_allele.exposure = c("A", "G", "C"),
    other_allele.exposure = c("G", "T", "A"),
    pval.exposure = c(1e-10, 1e-8, 1e-9),
    eaf.exposure = c(0.3, 0.4, 0.5),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = c("1", "1", "1"),
    pos.exposure = c(1000, 2000, 3000),
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = c("rs1", "rs2", "rs3"),
    beta = c(0.1, 0.05, 0.08),
    se = c(0.05, 0.03, 0.04),
    pval = c(0.01, 0.1, 0.05),
    eaf = c(0.3, 0.4, 0.5),
    effect_allele = c("A", "G", "C"),
    other_allele = c("G", "T", "A"),
    chr = c("1", "1", "1"),
    pos = c(1000, 2000, 3000),
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = c("rs1", "rs2", "rs3"),
      methods = c("ivw", "heterogeneity", "loo")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_false(is.null(result$loo))
  # Per-SNP rows plus the pooled "All" row
  expect_equal(nrow(result$loo), 4)
  expect_true("All" %in% result$loo$SNP)

  expect_false(is.null(result$heterogeneity))
})

test_that("run_mr skips heterogeneity and loo with 1 instrument", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.5,
    se.exposure = 0.1,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-10,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = 1000,
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = "rs1",
    beta = 0.1,
    se = 0.05,
    pval = 0.01,
    eaf = 0.3,
    effect_allele = "A",
    other_allele = "G",
    chr = "1",
    pos = 1000,
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = "rs1",
      methods = c("ivw", "heterogeneity", "loo")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_null(result$heterogeneity)
  expect_null(result$loo)
  expect_true("heterogeneity" %in% names(result$methods_skipped))
  expect_true("loo" %in% names(result$methods_skipped))
})

# --- flip_beta ---------------------------------------------------------------

test_that("run_mr flip_beta negates results, instruments, and loo consistently", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = c("rs1", "rs2", "rs3"),
    beta.exposure = c(0.5, 0.3, 0.4),
    se.exposure = c(0.1, 0.1, 0.1),
    effect_allele.exposure = c("A", "G", "C"),
    other_allele.exposure = c("G", "T", "A"),
    pval.exposure = c(1e-10, 1e-8, 1e-9),
    eaf.exposure = c(0.3, 0.4, 0.5),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = c("1", "1", "1"),
    pos.exposure = c(1000, 2000, 3000),
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = c("rs1", "rs2", "rs3"),
    beta = c(0.1, 0.05, 0.08),
    se = c(0.05, 0.03, 0.04),
    pval = c(0.01, 0.1, 0.05),
    eaf = c(0.3, 0.4, 0.5),
    effect_allele = c("A", "G", "C"),
    other_allele = c("G", "T", "A"),
    chr = c("1", "1", "1"),
    pos = c(1000, 2000, 3000),
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  run_args <- list(
    exposure = exposure,
    exposure_id = "test_exp",
    outcome = outcome,
    outcome_id = "test_out",
    instruments = c("rs1", "rs2", "rs3"),
    methods = c("ivw", "egger", "heterogeneity", "loo")
  )

  suppressMessages({
    unflipped <- do.call(run_mr, run_args)
    flipped <- do.call(run_mr, c(run_args, list(flip_beta = TRUE)))
  })

  expect_equal(flipped$results$b, -unflipped$results$b)
  expect_equal(flipped$results$or, 1 / unflipped$results$or)
  expect_equal(
    flipped$instruments$beta.exposure,
    -unflipped$instruments$beta.exposure
  )
  expect_equal(
    flipped$loo$b[flipped$loo$SNP != "All"],
    -unflipped$loo$b[unflipped$loo$SNP != "All"]
  )
  expect_true(flipped$params$flip_beta)
  expect_false(unflipped$params$flip_beta)

  # Sign-invariant outputs are unaffected by the flip
  expect_equal(flipped$f_stats, unflipped$f_stats)
  expect_equal(flipped$heterogeneity$Q, unflipped$heterogeneity$Q)

  # Audit trail: beta.exposure_original preserves the pre-transform value,
  # and effect_allele/other_allele/eaf are untouched by the flip (this is a
  # phenotype transformation, not an allele reorientation).
  expect_equal(
    flipped$instruments$beta.exposure_original,
    exposure$beta.exposure
  )
  expect_equal(
    unflipped$instruments$beta.exposure_original,
    exposure$beta.exposure
  )
  expect_equal(
    flipped$instruments$effect_allele.exposure,
    unflipped$instruments$effect_allele.exposure
  )
  expect_equal(
    flipped$instruments$other_allele.exposure,
    unflipped$instruments$other_allele.exposure
  )
  expect_equal(
    flipped$instruments$eaf.exposure,
    unflipped$instruments$eaf.exposure
  )

  # exposure_transformed/exposure_analysis_name on every output with an
  # "exposure" column
  expect_true(all(flipped$results$exposure_transformed))
  expect_true(all(flipped$instruments$exposure_transformed))
  expect_true(all(flipped$loo$exposure_transformed))
  expect_false(any(unflipped$results$exposure_transformed))

  expect_true(all(flipped$results$exposure_analysis_name == "-1 x test_exp"))
  expect_true(all(unflipped$results$exposure_analysis_name == "test_exp"))
})

test_that("run_mr exposure_label overrides the auto-generated exposure_analysis_name", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.5,
    se.exposure = 0.1,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-10,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = 1000,
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = "rs1",
    beta = 0.1,
    se = 0.05,
    pval = 0.01,
    eaf = 0.3,
    effect_allele = "A",
    other_allele = "G",
    chr = "1",
    pos = 1000,
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = "rs1",
      methods = "ivw",
      flip_beta = TRUE,
      exposure_label = "Genetically proxied test inhibition"
    )
  })

  expect_true(all(
    result$results$exposure_analysis_name ==
      "Genetically proxied test inhibition"
  ))
  expect_true(all(
    result$instruments$exposure_analysis_name ==
      "Genetically proxied test inhibition"
  ))
  expect_equal(
    result$params$exposure_label,
    "Genetically proxied test inhibition"
  )
})

test_that("run_mr flip_beta = FALSE (default) matches omitting the argument", {
  skip_if_not_installed("TwoSampleMR")

  exposure <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.5,
    se.exposure = 0.1,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-10,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = 1000,
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = "rs1",
    beta = 0.1,
    se = 0.05,
    pval = 0.01,
    eaf = 0.3,
    effect_allele = "A",
    other_allele = "G",
    chr = "1",
    pos = 1000,
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    default_call <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = "rs1",
      methods = "ivw"
    )
    explicit_false <- run_mr(
      exposure = exposure,
      exposure_id = "test_exp",
      outcome = outcome,
      outcome_id = "test_out",
      instruments = "rs1",
      methods = "ivw",
      flip_beta = FALSE
    )
  })

  expect_equal(default_call$results$b, explicit_false$results$b)
})
