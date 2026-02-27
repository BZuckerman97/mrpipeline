# Tests for internal helper functions in R/helpers.R

# --- eaf_to_maf --------------------------------------------------------------

test_that("eaf_to_maf converts frequencies correctly", {
  expect_equal(eaf_to_maf(0.1), 0.1)
  expect_equal(eaf_to_maf(0.9), 0.1)
  expect_equal(eaf_to_maf(0.5), 0.5)
  expect_equal(eaf_to_maf(0.01), 0.01)
  expect_equal(eaf_to_maf(0.99), 0.01)
})

test_that("eaf_to_maf is vectorised", {
  result <- eaf_to_maf(c(0.1, 0.5, 0.9))
  expect_equal(result, c(0.1, 0.5, 0.1))
})

test_that("eaf_to_maf handles edge cases", {
  expect_equal(eaf_to_maf(0), 0)
  expect_equal(eaf_to_maf(1), 0)
  expect_true(is.na(eaf_to_maf(NA)))
})

# --- resolve_sample_size -----------------------------------------------------

test_that("resolve_sample_size prefers explicit_n", {
  result <- resolve_sample_size(
    explicit_n = 5000,
    data_column = c(1000, 2000, 3000)
  )
  expect_identical(result, 5000L)
})

test_that("resolve_sample_size uses median of data column when no explicit_n", {
  result <- resolve_sample_size(
    explicit_n = NULL,
    data_column = c(1000, 2000, 3000),
    label = "test"
  )
  expect_identical(result, 2000L)
})

test_that("resolve_sample_size returns NULL when nothing available", {
  result <- resolve_sample_size(
    explicit_n = NULL,
    data_column = NULL
  )
  expect_null(result)
})

test_that("resolve_sample_size handles NA in data column", {
  result <- resolve_sample_size(
    explicit_n = NULL,
    data_column = c(NA, 2000, NA),
    label = "test"
  )
  expect_identical(result, 2000L)

  # All NAs
  result_all_na <- resolve_sample_size(
    explicit_n = NULL,
    data_column = c(NA, NA)
  )
  expect_null(result_all_na)
})

test_that("resolve_sample_size coerces to integer", {
  result <- resolve_sample_size(explicit_n = 5000.7)
  expect_identical(result, 5000L)
})

# --- harmonise_and_filter -----------------------------------------------------

test_that("harmonise_and_filter returns filtered deduplicated data", {
  skip_if_not_installed("TwoSampleMR")

  # Minimal exposure data
  exposure <- data.frame(
    SNP = c("rs1", "rs2", "rs1"),
    beta.exposure = c(0.1, 0.2, 0.15),
    se.exposure = c(0.05, 0.05, 0.05),
    effect_allele.exposure = c("A", "G", "A"),
    other_allele.exposure = c("G", "T", "G"),
    pval.exposure = c(1e-5, 1e-3, 1e-4),
    eaf.exposure = c(0.3, 0.4, 0.3),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  outcome <- data.frame(
    SNP = c("rs1", "rs2", "rs1"),
    beta.outcome = c(0.05, 0.1, 0.06),
    se.outcome = c(0.02, 0.03, 0.02),
    effect_allele.outcome = c("A", "G", "A"),
    other_allele.outcome = c("G", "T", "G"),
    pval.outcome = c(0.01, 0.05, 0.02),
    eaf.outcome = c(0.3, 0.4, 0.3),
    outcome = "test_out",
    id.outcome = "out1",
    mr_keep.outcome = TRUE,
    pval_origin.outcome = "reported",
    stringsAsFactors = FALSE
  )

  result <- harmonise_and_filter(exposure, outcome)

  # Should have no duplicate SNPs

  expect_false(any(duplicated(result$SNP)))
  # Should only contain mr_keep == TRUE rows
  expect_true(all(result$mr_keep))
})

# --- align_to_ld_matrix ------------------------------------------------------

test_that("align_to_ld_matrix subsets and reorders correctly", {
  data <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs4"),
    beta = c(0.1, 0.2, 0.3, 0.4),
    stringsAsFactors = FALSE
  )

  ld <- matrix(1, nrow = 3, ncol = 3)
  rownames(ld) <- colnames(ld) <- c("rs2", "rs3", "rs1")

  result <- align_to_ld_matrix(data, ld)

  # intersect preserves order from data: rs1, rs2, rs3
  expect_equal(result$data$SNP, c("rs1", "rs2", "rs3"))
  expect_equal(nrow(result$ld_matrix), 3)
  expect_equal(rownames(result$ld_matrix), c("rs1", "rs2", "rs3"))
  # Both are in the same order
  expect_equal(result$data$SNP, rownames(result$ld_matrix))
})

test_that("align_to_ld_matrix errors when no shared SNPs", {
  data <- data.frame(SNP = c("rs1", "rs2"), stringsAsFactors = FALSE)
  ld <- matrix(1, nrow = 1, ncol = 1)
  rownames(ld) <- colnames(ld) <- "rs99"

  expect_error(align_to_ld_matrix(data, ld), "No SNPs in common")
})

# --- compute_ld_matrix --------------------------------------------------------
# Integration test: requires local bfile + plink

test_that("compute_ld_matrix works with local reference panel", {
  bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")
  skip_if_not(
    file.exists(paste0(bfile, ".bed")),
    "LD reference panel not available"
  )

  # Get some SNPs from the .bim file
  bim <- utils::read.table(
    paste0(bfile, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  test_snps <- head(bim$V2, 3)

  result <- compute_ld_matrix(test_snps, bfile)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(test_snps))
  # Row/col names should be clean rsIDs (no allele suffixes)
  expect_false(any(grepl("_", rownames(result))))
})

# --- clump_instruments --------------------------------------------------------
# Integration test: requires local bfile + plink

test_that("clump_instruments works with local reference panel", {
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
  test_dat <- data.frame(
    rsid = head(bim$V2, 10),
    pval = runif(min(10, nrow(bim)), min = 1e-10, max = 1e-5),
    id = "test",
    stringsAsFactors = FALSE
  )

  result <- clump_instruments(
    dat = test_dat,
    rsq_thresh = 0.1,
    bfile = bfile
  )

  expect_true(is.data.frame(result))
  expect_true("rsid" %in% colnames(result))
  expect_true(nrow(result) <= nrow(test_dat))
})
