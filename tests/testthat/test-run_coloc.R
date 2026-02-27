# Tests for run_coloc()

# --- Argument validation ----------------------------------------------------

test_that("run_coloc errors when bfile is missing", {
  expect_error(
    run_coloc(
      exposure = data.frame(),
      outcome = data.frame(),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000
    ),
    "bfile"
  )
})

test_that("run_coloc validates methods argument", {
  expect_error(
    run_coloc(
      exposure = data.frame(),
      outcome = data.frame(),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000,
      bfile = "/fake/path",
      methods = "not_a_method"
    ),
    "should be one of"
  )
})

test_that("run_coloc errors when outcome_type = 'cc' without outcome_s", {
  expect_error(
    run_coloc(
      exposure = data.frame(
        SNP = "rs1",
        chr.exposure = "1",
        pos.exposure = 1500,
        beta.exposure = 0.1,
        se.exposure = 0.05,
        stringsAsFactors = FALSE
      ),
      outcome = data.frame(),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000,
      bfile = "/fake/path",
      outcome_type = "cc"
    ),
    "outcome_s"
  )
})

test_that("run_coloc errors when exposure_type = 'cc' without exposure_s", {
  expect_error(
    run_coloc(
      exposure = data.frame(
        SNP = "rs1",
        chr.exposure = "1",
        pos.exposure = 1500,
        beta.exposure = 0.1,
        se.exposure = 0.05,
        stringsAsFactors = FALSE
      ),
      outcome = data.frame(),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000,
      bfile = "/fake/path",
      exposure_type = "cc"
    ),
    "exposure_s"
  )
})

# --- Status returns (no external data needed) --------------------------------

test_that("run_coloc returns no_snps_in_region when exposure has no SNPs in window", {
  exposure <- data.frame(
    SNP = "rs1",
    chr.exposure = "2",
    pos.exposure = 500000,
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-5,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- run_coloc(
      exposure = exposure,
      outcome = data.frame(
        rsids = "rs1",
        chr = "1",
        pos = 1500,
        beta = 0.1,
        se = 0.05,
        eaf = 0.3,
        pval = 0.01,
        n = 5000,
        effect_allele = "A",
        other_allele = "G",
        stringsAsFactors = FALSE
      ),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000,
      bfile = "/fake/path",
      methods = "abf"
    ),
    "No exposure SNPs"
  )
  expect_s3_class(result, "coloc_result")
  expect_equal(result$status, "no_snps_in_region")
})

test_that("run_coloc returns no_snps_in_region when outcome has no SNPs in window", {
  exposure <- data.frame(
    SNP = "rs1",
    chr.exposure = "1",
    pos.exposure = 1500,
    beta.exposure = 0.1,
    se.exposure = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    pval.exposure = 1e-5,
    eaf.exposure = 0.3,
    exposure = "test_exp",
    id.exposure = "exp1",
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- run_coloc(
      exposure = exposure,
      outcome = data.frame(
        rsids = "rs1",
        chr = "2",
        pos = 500000,
        beta = 0.1,
        se = 0.05,
        eaf = 0.3,
        pval = 0.01,
        n = 5000,
        effect_allele = "A",
        other_allele = "G",
        stringsAsFactors = FALSE
      ),
      gene_chr = 1,
      gene_start = 1000,
      gene_end = 2000,
      bfile = "/fake/path",
      methods = "abf"
    ),
    "No outcome SNPs"
  )
  expect_s3_class(result, "coloc_result")
  expect_equal(result$status, "no_snps_in_region")
})

# --- coloc_result S3 class ---------------------------------------------------

test_that("coloc_result print works with non-success status", {
  res <- new_coloc_result(
    status = "no_snps_in_region",
    status_reason = "No exposure SNPs in region chr1:990-2010"
  )
  expect_message(print(res), "no_snps_in_region")
  expect_message(print(res), "No exposure SNPs")
})

test_that("coloc_result print works with ABF results", {
  res <- new_coloc_result(
    coloc_abf = list(
      summary = c(
        nsnps = 50,
        PP.H0.abf = 0.01,
        PP.H1.abf = 0.02,
        PP.H2.abf = 0.03,
        PP.H3.abf = 0.04,
        PP.H4.abf = 0.90
      )
    ),
    n_snps = 50L
  )
  expect_message(print(res), "PP.H4")
  expect_message(print(res), "50 SNPs")
})

test_that("coloc_result print works with no method results", {
  res <- new_coloc_result(n_snps = 10L)
  expect_message(print(res), "no method results")
})

test_that("coloc_result summary works with ABF results", {
  res <- new_coloc_result(
    coloc_abf = list(
      summary = c(
        nsnps = 50,
        PP.H0.abf = 0.01,
        PP.H1.abf = 0.02,
        PP.H2.abf = 0.03,
        PP.H3.abf = 0.04,
        PP.H4.abf = 0.90
      )
    ),
    n_snps = 50L,
    params = list()
  )
  expect_message(summary(res), "Colocalization Results")
  expect_message(summary(res), "PP.H4")
})

test_that("coloc_result summary shows skipped methods", {
  res <- new_coloc_result(
    n_snps = 10L,
    methods_skipped = c(signals = "SuSiE failed"),
    params = list()
  )
  expect_message(summary(res), "Skipped methods")
  expect_message(summary(res), "SuSiE failed")
})

test_that("coloc_result summary shows non-success status", {
  res <- new_coloc_result(
    status = "too_few_snps",
    status_reason = "Only 2 SNPs after harmonisation (need >= 3)",
    params = list()
  )
  expect_message(summary(res), "too_few_snps")
})

# --- Integration tests (require bfile) ---------------------------------------

test_that("run_coloc ABF-only returns correct coloc_result structure", {
  skip_if_not_installed("TwoSampleMR")
  skip_if_not_installed("coloc")

  bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )

  # Read SNPs from bim file to build realistic test data

  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  test_snps <- head(bim$V2, 20)
  test_chr <- as.character(bim$V1[1])
  test_positions <- bim$V4[seq_along(test_snps)]
  min_pos <- min(test_positions)
  max_pos <- max(test_positions)

  exposure <- data.frame(
    SNP = test_snps,
    chr.exposure = test_chr,
    pos.exposure = test_positions,
    beta.exposure = rnorm(length(test_snps), 0, 0.1),
    se.exposure = runif(length(test_snps), 0.01, 0.05),
    effect_allele.exposure = bim$V5[seq_along(test_snps)],
    other_allele.exposure = bim$V6[seq_along(test_snps)],
    pval.exposure = runif(length(test_snps), 1e-8, 1e-3),
    eaf.exposure = runif(length(test_snps), 0.1, 0.9),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = test_snps,
    chr = test_chr,
    pos = test_positions,
    beta = rnorm(length(test_snps), 0, 0.1),
    se = runif(length(test_snps), 0.01, 0.05),
    eaf = runif(length(test_snps), 0.1, 0.9),
    pval = runif(length(test_snps), 1e-5, 0.5),
    n = 20000,
    effect_allele = bim$V5[seq_along(test_snps)],
    other_allele = bim$V6[seq_along(test_snps)],
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages({
    result <- run_coloc(
      exposure = exposure,
      outcome = outcome,
      gene_chr = test_chr,
      gene_start = min_pos,
      gene_end = max_pos,
      coloc_window = 10000L,
      outcome_n = 20000,
      bfile = bfile,
      methods = "abf"
    )
  })

  expect_s3_class(result, "coloc_result")
  expect_equal(result$status, "success")
  expect_true(!is.null(result$coloc_abf))
  expect_true("PP.H4.abf" %in% names(result$coloc_abf$summary))
  expect_true(result$n_snps > 0)
  expect_true(nrow(result$harmonised_data) > 0)
})

test_that("run_coloc skips prop_test when susie not requested", {
  skip_if_not_installed("TwoSampleMR")
  skip_if_not_installed("coloc")

  bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )

  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  test_snps <- head(bim$V2, 20)
  test_chr <- as.character(bim$V1[1])
  test_positions <- bim$V4[seq_along(test_snps)]

  exposure <- data.frame(
    SNP = test_snps,
    chr.exposure = test_chr,
    pos.exposure = test_positions,
    beta.exposure = rnorm(length(test_snps), 0, 0.1),
    se.exposure = runif(length(test_snps), 0.01, 0.05),
    effect_allele.exposure = bim$V5[seq_along(test_snps)],
    other_allele.exposure = bim$V6[seq_along(test_snps)],
    pval.exposure = runif(length(test_snps), 1e-8, 1e-3),
    eaf.exposure = runif(length(test_snps), 0.1, 0.9),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    rsids = test_snps,
    chr = test_chr,
    pos = test_positions,
    beta = rnorm(length(test_snps), 0, 0.1),
    se = runif(length(test_snps), 0.01, 0.05),
    eaf = runif(length(test_snps), 0.1, 0.9),
    pval = runif(length(test_snps), 1e-5, 0.5),
    n = 20000,
    effect_allele = bim$V5[seq_along(test_snps)],
    other_allele = bim$V6[seq_along(test_snps)],
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )

  suppressMessages(suppressWarnings({
    result <- run_coloc(
      exposure = exposure,
      outcome = outcome,
      gene_chr = test_chr,
      gene_start = min(test_positions),
      gene_end = max(test_positions),
      coloc_window = 10000L,
      outcome_n = 20000,
      bfile = bfile,
      methods = c("abf", "prop_test")
    )
  }))

  expect_s3_class(result, "coloc_result")
  expect_true("prop_test" %in% names(result$methods_skipped))
})
