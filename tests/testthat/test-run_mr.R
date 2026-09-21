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

test_that("run_mr validates allele_check argument", {
  expect_error(
    run_mr(
      exposure = data.frame(),
      exposure_id = "test",
      outcome = data.frame(),
      outcome_id = "test",
      allele_check = "bogus"
    ),
    "allele_check"
  )
})

test_that("run_mr validates harmonise_action argument", {
  for (bad in list(0, 4, "2", c(2, 3), NA)) {
    expect_error(
      run_mr(
        exposure = data.frame(),
        exposure_id = "test",
        outcome = data.frame(),
        outcome_id = "test",
        harmonise_action = bad
      ),
      "action"
    )
  }
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
      ld_corrected = FALSE,
      model = "random",
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
      ld_corrected = FALSE,
      model = c("random", "random"),
      stringsAsFactors = FALSE
    ),
    f_stats = list(per_snp = rep(30, 5), mean = 30, min = 25),
    methods_skipped = c(presso = "Requires >= 3 instruments"),
    params = list(exposure_id = "exp1", outcome_id = "out1")
  )
  expect_message(summary(res), "MR Results")
})

test_that("mr_result summary reports the harmonisation breakdown", {
  raw <- data.frame(
    SNP = paste0("rs", 1:5),
    mr_keep = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    palindromic = c(FALSE, FALSE, TRUE, TRUE, FALSE),
    ambiguous = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    remove = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  res <- new_mr_result(
    results = data.frame(
      exposure = "exp1",
      outcome = "out1",
      method = "IVW",
      nsnp = 3,
      b = 0.1,
      se = 0.05,
      pval = 0.01,
      ld_corrected = FALSE,
      model = "random",
      stringsAsFactors = FALSE
    ),
    f_stats = list(per_snp = rep(30, 3), mean = 30, min = 25),
    harmonisation = raw,
    params = list(exposure_id = "exp1", outcome_id = "out1")
  )
  expect_message(summary(res), "Harmonisation")
  expect_message(summary(res), "5 candidate SNPs -> 3 kept, 2 dropped")
  expect_message(summary(res), "ambiguous")
})

test_that("mr_result summary omits harmonisation when there is none", {
  res <- new_mr_result(
    results = data.frame(
      exposure = "exp1",
      outcome = "out1",
      method = "IVW",
      nsnp = 3,
      b = 0.1,
      se = 0.05,
      pval = 0.01,
      ld_corrected = FALSE,
      model = "random",
      stringsAsFactors = FALSE
    ),
    f_stats = list(per_snp = rep(30, 3), mean = 30, min = 25),
    params = list(exposure_id = "exp1", outcome_id = "out1")
  )
  # Older result objects carry no harmonisation field at all
  expect_message(summary(res), "MR Results")
  expect_false(any(grepl(
    "Harmonisation",
    capture_messages(summary(res))
  )))
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
      methods = c("ivw_random", "egger", "weighted_median")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "success")
  expect_equal(nrow(result$results), 1)
  expect_equal(result$results$method, "Wald ratio")

  # run_mr() carries the unfiltered harmonisation through (issue #17):
  # `instruments` stays the kept set, `harmonisation` explains the rest
  expect_true(all(
    c("mr_keep", "palindromic", "ambiguous", "remove") %in%
      names(result$harmonisation)
  ))
  expect_gte(nrow(result$harmonisation), nrow(result$instruments))
  expect_true(all(result$instruments$mr_keep))

  # IVW, egger, weighted_median should be skipped
  expect_true(all(
    c("ivw_random", "egger", "weighted_median") %in%
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
      methods = c("ivw_random", "egger", "weighted_median", "presso")
    )
  })

  expect_s3_class(result, "mr_result")

  # IVW should work, and the row says which estimator ran
  expect_true("IVW (random effects)" %in% result$results$method)
  ivw_row <- result$results[result$results$method == "IVW (random effects)", ]
  expect_equal(ivw_row$model, "random")
  expect_false(ivw_row$ld_corrected)

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
      methods = c("ivw_random", "heterogeneity", "loo")
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
      methods = c("ivw_random", "heterogeneity", "loo")
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
      methods = c("ivw_random", "heterogeneity", "loo")
    )
  })

  expect_s3_class(result, "mr_result")
  expect_null(result$heterogeneity)
  expect_null(result$loo)
  expect_true("heterogeneity" %in% names(result$methods_skipped))
  expect_true("loo" %in% names(result$methods_skipped))
})

# --- Allele orientation check -----------------------------------------------

test_that("run_mr detects a swapped-allele outcome even with 3 instruments", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_gwas_fixture()
  run <- function(outcome, ...) {
    suppressMessages(run_mr(
      exposure = f$exposure,
      exposure_id = "exp",
      outcome = outcome,
      outcome_id = "out",
      instruments = f$instruments,
      methods = "ivw_random",
      verbose = FALSE,
      ...
    ))
  }

  expect_error(run(f$bug), class = "mrpipeline_allele_check_error")
  rec <- last_allele_check()
  expect_equal(rec$status, "fail")
  expect_equal(rec$n_sampled, 37L) # the sampled path, not just the 3 instruments

  # Correctly labelled outcome runs and records a pass
  result <- run(f$ok)
  expect_s3_class(result, "mr_result")
  expect_equal(result$status, "success")
  expect_equal(nrow(result$instruments), 3L)
  expect_equal(last_allele_check()$status, "pass")
  expect_true("allele_check" %in% names(result$timing))
  expect_equal(result$params$allele_check, "error")

  # "warn" carries on and records the mode
  expect_warning(
    result <- run(f$bug, allele_check = "warn"),
    class = "mrpipeline_allele_check_warning"
  )
  expect_equal(result$status, "success")
  expect_equal(result$params$allele_check, "warn")
})

# --- LD correction (integration: requires bfile + plink) --------------------

# Resolve the bundled reference panel, skipping cleanly where it is absent.
ld_bfile <- function() {
  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
  testthat::skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )
  bfile
}

cd40_region <- list(chromosome = "20", start = 44746911, end = 44758502)

# The bundled panel is real 1000 Genomes EUR LD, so the default
# rsq_thresh = 0.001 clumps the CD40 cis region down to one instrument. At
# 0.3 seven survive, with pairwise |r| up to 0.53 -- enough LD that the
# corrected and uncorrected fits differ materially.
cd40_rsq <- 0.3

# Three of the CD40 cis instruments: none palindromic, pairwise r 0.45,
# -0.27 and -0.16 in the panel (rcond of R 0.36), and ConMix converges on
# them (it does not on every triple -- "wrong sign in by argument" from
# mr_conmix()). Avoid rs4810485: it is in perfect LD with rs1883832 (r = 1),
# which makes the GLS weight matrix singular.
cd40_three_snps <- function() {
  c("rs1883832", "rs4810486", "rs77048809")
}

# run_mr() on the bundled CD40 -> SjD data with PLINK's console output
# captured; warnings are returned alongside the result so tests can assert
# on them without the nested expect_warning() dance.
run_cd40 <- function(...) {
  warnings <- character()
  result <- NULL
  invisible(utils::capture.output(
    result <- withCallingHandlers(
      suppressMessages(run_mr(
        exposure = cd40_exposure,
        exposure_id = "CD40",
        outcome = sjogren_outcome,
        outcome_id = "SjD",
        verbose = FALSE,
        ...
      )),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  ))
  list(result = result, warnings = warnings)
}

test_that("ivw_fixed with ld_correct = TRUE is the hand-computed GLS estimate", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
  )
  res <- out$result
  expect_equal(res$status, "success")
  expect_true(all(c("ld_corrected", "model") %in% names(res$results)))

  h <- res$instruments
  ld <- res$ld_matrix
  expect_equal(nrow(h), nrow(ld))
  expect_gt(nrow(h), 2)

  # GLS with weight matrix diag(se_y) R diag(se_y):
  #   b  = (bx' O^-1 bx)^-1 bx' O^-1 by,  se = sqrt((bx' O^-1 bx)^-1)
  omega <- diag(h$se.outcome) %*% ld %*% diag(h$se.outcome)
  omega_inv <- solve(omega)
  bx <- h$beta.exposure
  by <- h$beta.outcome
  b_gls <- as.numeric(
    solve(t(bx) %*% omega_inv %*% bx) %*% (t(bx) %*% omega_inv %*% by)
  )
  se_gls <- sqrt(as.numeric(solve(t(bx) %*% omega_inv %*% bx)))

  fixed <- res$results[res$results$method == "IVW (fixed effects)", ]
  expect_equal(fixed$b, b_gls, tolerance = 1e-8)
  expect_equal(fixed$se, se_gls, tolerance = 1e-8)
  expect_true(fixed$ld_corrected)
  expect_equal(fixed$model, "fixed")

  # Random effects: same estimate, standard error never below fixed
  random <- res$results[res$results$method == "IVW (random effects)", ]
  expect_equal(random$b, fixed$b)
  expect_gte(random$se, fixed$se)
  expect_true(random$ld_corrected)
  expect_equal(random$model, "random")

  # Egger has a correlated form; weighted median does not
  egger <- res$results[res$results$method == "MR Egger", ]
  expect_true(egger$ld_corrected)
  expect_equal(egger$model, "random")
  wm <- res$results[res$results$method == "Weighted median", ]
  expect_false(wm$ld_corrected)
  expect_true(is.na(wm$model))

  # ... and said so by name
  expect_true(any(grepl("weighted_median", out$warnings)))
  expect_false(any(grepl("ivw_random|ivw_fixed|egger", out$warnings)))
})

test_that("ld_correct = FALSE records every row as uncorrected and warns about nothing", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
  )
  expect_equal(out$result$status, "success")
  expect_false(any(out$result$results$ld_corrected))
  expect_equal(
    out$result$results$model,
    c("random", "fixed", "random", NA_character_)
  )
  expect_length(out$warnings, 0)
})

test_that("ivw_random stays random effects below 4 instruments when LD-corrected", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instruments = cd40_three_snps(),
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "ivw_fixed")
  )
  res <- out$result
  expect_equal(res$status, "success")
  expect_equal(nrow(res$instruments), 3)
  # MendelianRandomization's "default" model would be fixed effects here;
  # the shortcut pins random (issue #27)
  expect_equal(res$results$model, c("random", "fixed"))
  expect_true(all(res$results$ld_corrected))
  expect_gte(res$results$se[1], res$results$se[2])
  expect_equal(res$results$b[1], res$results$b[2])
})

test_that("methods with no LD-corrected form warn by name and are recorded as uncorrected", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instruments = cd40_three_snps(),
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "weighted_median", "conmix")
  )
  res <- out$result
  expect_true(any(grepl(
    '"weighted_median" has no LD-corrected form',
    out$warnings
  )))
  expect_true(any(grepl('"conmix" has no LD-corrected form', out$warnings)))

  by_method <- stats::setNames(res$results$ld_corrected, res$results$method)
  expect_true(by_method[["IVW (random effects)"]])
  expect_false(by_method[["Weighted median"]])
  expect_false(by_method[["ConMix"]])

  # summary() lists both sides; print() tags the LD-corrected primary row
  expect_message(summary(res), "Applied to: IVW \\(random effects\\)")
  expect_message(summary(res), "Not applied to: Weighted median, ConMix")
  expect_message(print(res), "LD-corrected")
})

test_that("a method skipped for instrument count does not also warn about LD", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instruments = cd40_three_snps()[1:2],
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "weighted_median")
  )
  res <- out$result
  expect_equal(nrow(res$instruments), 2)
  expect_match(res$methods_skipped[["weighted_median"]], "Requires >= 3")
  expect_false(any(grepl("weighted_median", out$warnings)))
  expect_equal(res$results$method, "IVW (random effects)")
  expect_true(res$results$ld_corrected)
})

test_that("a single instrument with ld_correct = TRUE gives the Wald ratio and says LD is not applicable", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  # rs1883832's panel A1 is "T": on a single-SNP .bim ieugwasr reads that as
  # logical TRUE, which compute_ld_matrix() must map back to "T" or the SNP
  # is dropped at alignment and nothing is left to analyse.
  out <- run_cd40(
    instruments = "rs1883832",
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "weighted_median")
  )
  res <- out$result
  expect_equal(res$status, "success")
  expect_equal(res$results$method, "Wald ratio")
  expect_false(res$results$ld_corrected)
  expect_true(is.na(res$results$model))
  expect_match(res$methods_skipped[["ivw_random"]], "Only 1 instrument")
  expect_false(any(grepl("LD-corrected form", out$warnings)))
  expect_message(summary(res), "Not applicable: 1 instrument")
})

test_that("on the bundled panel the diagnostics match the correlated fits exactly", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  # Checks the plumbing against hand-computed correlated quantities; the
  # next test checks that the correction is material on this panel.
  out <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "egger", "heterogeneity", "loo")
  )
  res <- out$result
  expect_equal(res$status, "success")

  h <- res$instruments
  ld <- res$ld_matrix
  n <- nrow(h)
  omega <- diag(h$se.outcome) %*% ld %*% diag(h$se.outcome)
  w <- solve(omega)
  bx <- h$beta.exposure
  by <- h$beta.outcome
  b_gls <- as.numeric(solve(t(bx) %*% w %*% bx) %*% (t(bx) %*% w %*% by))
  r <- by - b_gls * bx
  q_hand <- as.numeric(t(r) %*% w %*% r)

  ht <- res$heterogeneity
  expect_equal(ht$method, c("MR Egger", "Inverse variance weighted"))
  expect_true(all(ht$ld_corrected))
  ivw_row <- ht[ht$method == "Inverse variance weighted", ]
  expect_equal(ivw_row$Q, q_hand, tolerance = 1e-8)
  expect_equal(ivw_row$Q_df, n - 1)
  expect_equal(
    ivw_row$Q_pval,
    stats::pchisq(q_hand, n - 1, lower.tail = FALSE),
    tolerance = 1e-8
  )
  expect_equal(ht$Q_df[ht$method == "MR Egger"], n - 2)

  input <- MendelianRandomization::mr_input(
    bx = bx,
    bxse = h$se.exposure,
    by = by,
    byse = h$se.outcome,
    correlation = ld
  )
  egger <- MendelianRandomization::mr_egger(input, correl = TRUE)
  expect_equal(
    res$pleiotropy$egger_intercept,
    egger@Intercept,
    tolerance = 1e-8
  )
  expect_equal(res$pleiotropy$pval, egger@Pvalue.Int, tolerance = 1e-8)
  expect_true(res$pleiotropy$ld_corrected)

  expect_equal(nrow(res$loo), n + 1)
  expect_true(all(res$loo$ld_corrected))
  all_row <- res$loo[res$loo$SNP == "All", ]
  ivw_random <- res$results[res$results$method == "IVW (random effects)", ]
  expect_equal(all_row$b, ivw_random$b, tolerance = 1e-8)
  expect_equal(all_row$se, ivw_random$se, tolerance = 1e-8)

  # summary() says which frames came from the correlated fits
  expect_message(
    summary(res),
    "Heterogeneity test \\(Cochran's Q\\) \\[LD-corrected\\]"
  )
  expect_message(summary(res), "Leave-one-out analysis \\[LD-corrected\\]")
  expect_message(summary(res), "Diagnostics from the correlated fits")

  # ... and the leave-one-out frame keeps TwoSampleMR's shape for plotting
  skip_if_not_installed("ggplot2")
  # TwoSampleMR::mr_leaveoneout_plot() returns one ggplot per exposure/outcome
  p <- plot(res, type = "loo")
  expect_type(p, "list")
  expect_s3_class(p[[1]], "ggplot")
})

test_that("perfectly correlated bundled instruments give a singular_ld_matrix status", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  # rs1883832 and rs4810485 are in perfect LD in 1000 Genomes EUR, so the
  # GLS weight matrix is exactly singular and MendelianRandomization's
  # solve() used to abort the whole call (issue #36).
  out <- run_cd40(
    instruments = c("rs1883832", "rs4810485", "rs4810486"),
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "egger")
  )
  res <- out$result
  expect_equal(res$status, "singular_ld_matrix")
  expect_match(res$status_reason, "rs1883832/rs4810485|rs4810485/rs1883832")
  expect_equal(nrow(res$results), 0)

  # Dropping one of the pair leaves a matrix that solves
  ok <- run_cd40(
    instruments = c("rs1883832", "rs4810486"),
    bfile = bfile,
    ld_correct = TRUE,
    methods = "ivw_random"
  )
  expect_equal(ok$result$status, "success")
  expect_true(ok$result$results$ld_corrected)
})

test_that("a singular matrix no method solves leaves the run alone", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  # weighted_median has no correlated form, so nothing solves the matrix
  out <- run_cd40(
    instruments = c("rs1883832", "rs4810485", "rs4810486"),
    bfile = bfile,
    ld_correct = TRUE,
    methods = "weighted_median"
  )
  expect_equal(out$result$status, "success")
  expect_false(out$result$results$ld_corrected)
  expect_true(any(grepl("near-singular", out$warnings)))
})

test_that("on the bundled panel LD correction changes the estimates and diagnostics", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  methods <- c("ivw_fixed", "heterogeneity")
  unc <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    methods = methods
  )$result
  corr <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    ld_correct = TRUE,
    methods = methods
  )$result
  expect_equal(corr$instruments$SNP, unc$instruments$SNP)

  # the clumped instruments are in real LD ...
  off_diag <- corr$ld_matrix[upper.tri(corr$ld_matrix)]
  expect_gt(max(abs(off_diag)), 0.5)

  # ... so the correlated fits move the estimate, its standard error and Q
  q_unc <- unc$heterogeneity$Q[
    unc$heterogeneity$method == "Inverse variance weighted"
  ]
  q_corr <- corr$heterogeneity$Q[
    corr$heterogeneity$method == "Inverse variance weighted"
  ]
  expect_gt(abs(q_corr - q_unc) / q_unc, 0.05)
  expect_gt(abs(corr$results$se - unc$results$se) / unc$results$se, 0.05)
  expect_gt(abs(corr$results$b - unc$results$b) / abs(unc$results$b), 0.05)
})

test_that("ld_correct = FALSE leaves the diagnostics as TwoSampleMR returns them, marked uncorrected", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instrument_region = cd40_region,
    rsq_thresh = cd40_rsq,
    bfile = bfile,
    methods = c("ivw_random", "egger", "heterogeneity", "loo")
  )
  res <- out$result
  expect_false(any(res$heterogeneity$ld_corrected))
  expect_false(res$pleiotropy$ld_corrected)
  expect_false(any(res$loo$ld_corrected))

  h <- res$instruments
  expect_equal(
    res$heterogeneity[, setdiff(names(res$heterogeneity), "ld_corrected")],
    TwoSampleMR::mr_heterogeneity(h),
    ignore_attr = TRUE
  )
  expect_equal(
    res$loo[, setdiff(names(res$loo), "ld_corrected")],
    TwoSampleMR::mr_leaveoneout(h),
    ignore_attr = TRUE
  )
  expect_false(any(grepl(
    "Cochran's Q\\) \\[LD-corrected\\]",
    capture_messages(summary(res))
  )))
})

test_that("two LD-corrected instruments give an IVW-only heterogeneity row", {
  skip_if_not_installed("TwoSampleMR")
  bfile <- ld_bfile()

  out <- run_cd40(
    instruments = cd40_three_snps()[1:2],
    bfile = bfile,
    ld_correct = TRUE,
    methods = c("ivw_random", "heterogeneity")
  )
  ht <- out$result$heterogeneity
  expect_equal(ht$method, "Inverse variance weighted")
  expect_equal(ht$Q_df, 1)
  expect_true(ht$ld_corrected)
})

# --- summary(): Steiger ------------------------------------------------------

test_that("summary() reports the per-SNP Steiger result it is given", {
  # $steiger holds TwoSampleMR::steiger_filtering() output -- one row per SNP
  # with a logical steiger_dir -- not directionality_test()'s single
  # correct_causal_direction, which summary() once read and printed as blank.
  res <- new_mr_result(
    results = data.frame(
      method = "IVW (random effects)",
      nsnp = 3L,
      b = 0.1,
      se = 0.05,
      pval = 0.05,
      model = "random",
      ld_corrected = FALSE
    ),
    steiger = data.frame(
      SNP = c("rs1", "rs2", "rs3"),
      steiger_dir = c(TRUE, TRUE, FALSE),
      steiger_pval = c(1e-10, 1e-5, 0.2)
    )
  )
  msgs <- paste(capture_messages(summary(res)), collapse = "")
  expect_match(msgs, "2/3 SNPs explain more variance in the exposure")
  expect_match(msgs, "Largest Steiger p-value: 0.2")
  expect_match(msgs, "Not in the expected direction: \"rs3\"")
})
