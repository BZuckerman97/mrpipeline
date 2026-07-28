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
# ld_matrix() is the list returned by compute_ld_matrix(): list(ld, alleles).
# alleles$ld_a1 is the allele the LD matrix's sign is anchored to.

make_ld_list <- function(snps, ld_a1, ld_a2, ld = NULL) {
  n <- length(snps)
  if (is.null(ld)) {
    ld <- diag(n)
  }
  rownames(ld) <- colnames(ld) <- snps
  list(
    ld = ld,
    alleles = data.frame(
      SNP = snps,
      ld_a1 = ld_a1,
      ld_a2 = ld_a2,
      stringsAsFactors = FALSE
    )
  )
}

test_that("align_to_ld_matrix subsets and reorders correctly", {
  data <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs4"),
    beta.exposure = c(0.1, 0.2, 0.3, 0.4),
    beta.outcome = c(0.05, 0.1, 0.15, 0.2),
    effect_allele.exposure = c("A", "A", "A", "A"),
    other_allele.exposure = c("G", "G", "G", "G"),
    eaf.exposure = c(0.3, 0.3, 0.3, 0.3),
    eaf.outcome = c(0.3, 0.3, 0.3, 0.3),
    stringsAsFactors = FALSE
  )

  ld_matrix <- make_ld_list(
    snps = c("rs2", "rs3", "rs1"),
    ld_a1 = c("A", "A", "A"),
    ld_a2 = c("G", "G", "G")
  )

  result <- suppressMessages(align_to_ld_matrix(data, ld_matrix))

  # intersect preserves order from data: rs1, rs2, rs3
  expect_equal(result$data$SNP, c("rs1", "rs2", "rs3"))
  expect_equal(nrow(result$ld_matrix), 3)
  expect_equal(rownames(result$ld_matrix), c("rs1", "rs2", "rs3"))
  # Both are in the same order
  expect_equal(result$data$SNP, rownames(result$ld_matrix))
  # All alleles matched the LD panel -- no flips
  expect_equal(result$ld_sign, c(1, 1, 1))
})

test_that("align_to_ld_matrix errors when no shared SNPs", {
  data <- data.frame(
    SNP = c("rs1", "rs2"),
    effect_allele.exposure = c("A", "A"),
    other_allele.exposure = c("G", "G"),
    stringsAsFactors = FALSE
  )
  ld_matrix <- make_ld_list(snps = "rs99", ld_a1 = "A", ld_a2 = "G")

  expect_error(align_to_ld_matrix(data, ld_matrix), "No SNPs in common")
})

test_that("align_to_ld_matrix flips ld_sign when LD panel allele order is swapped", {
  data <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.1, 0.2),
    beta.outcome = c(0.05, 0.1),
    effect_allele.exposure = c("A", "C"),
    other_allele.exposure = c("G", "T"),
    eaf.exposure = c(0.3, 0.2),
    eaf.outcome = c(0.3, 0.2),
    stringsAsFactors = FALSE
  )

  # rs1: LD panel's A1/A2 match the harmonised effect/other allele -> "match"
  # rs2: LD panel's A1/A2 are swapped relative to harmonised -> "flip"
  ld_matrix <- make_ld_list(
    snps = c("rs1", "rs2"),
    ld_a1 = c("A", "T"),
    ld_a2 = c("G", "C")
  )

  result <- suppressMessages(align_to_ld_matrix(data, ld_matrix))

  expect_equal(result$data$SNP, c("rs1", "rs2"))
  expect_equal(result$ld_sign, c(1, -1))
  # beta.exposure/beta.outcome themselves are untouched -- ld_sign is applied
  # separately when building coloc/SuSiE dataset objects, not by mutating data
  expect_equal(result$data$beta.exposure, c(0.1, 0.2))
})

test_that("align_to_ld_matrix drops SNPs unresolvable against the LD panel", {
  data <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.1, 0.2),
    beta.outcome = c(0.05, 0.1),
    effect_allele.exposure = c("A", "C"),
    other_allele.exposure = c("G", "T"),
    eaf.exposure = c(0.3, 0.2),
    eaf.outcome = c(0.3, 0.2),
    stringsAsFactors = FALSE
  )

  # rs2's LD panel alleles (A/G) share no letters with the harmonised C/T --
  # an unresolvable allele mismatch (e.g. a mislabelled/multi-allelic site).
  ld_matrix <- make_ld_list(
    snps = c("rs1", "rs2"),
    ld_a1 = c("A", "A"),
    ld_a2 = c("G", "G")
  )

  result <- suppressMessages(align_to_ld_matrix(data, ld_matrix))

  expect_equal(result$data$SNP, "rs1")
  expect_equal(nrow(result$ld_matrix), 1)
  expect_equal(length(result$ld_sign), 1)
})

test_that("align_to_ld_matrix drops ambiguous palindromic SNPs", {
  data <- data.frame(
    SNP = c("rs1", "rs2"),
    beta.exposure = c(0.1, 0.2),
    beta.outcome = c(0.05, 0.1),
    effect_allele.exposure = c("A", "A"),
    other_allele.exposure = c("G", "T"),
    # rs2 is a palindromic A/T SNP with EAF in the ambiguous 0.42-0.58 zone
    eaf.exposure = c(0.3, 0.5),
    eaf.outcome = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )

  ld_matrix <- make_ld_list(
    snps = c("rs1", "rs2"),
    ld_a1 = c("A", "A"),
    ld_a2 = c("G", "T")
  )

  result <- suppressMessages(align_to_ld_matrix(data, ld_matrix))

  expect_equal(result$data$SNP, "rs1")
})

test_that("align_to_ld_matrix keeps non-ambiguous palindromic SNPs", {
  data <- data.frame(
    SNP = "rs1",
    beta.exposure = 0.1,
    beta.outcome = 0.05,
    effect_allele.exposure = "A",
    other_allele.exposure = "T",
    # Far from 0.5 -- not ambiguous, should be resolved by allele matching
    eaf.exposure = 0.1,
    eaf.outcome = 0.1,
    stringsAsFactors = FALSE
  )

  ld_matrix <- make_ld_list(snps = "rs1", ld_a1 = "A", ld_a2 = "T")

  result <- suppressMessages(align_to_ld_matrix(data, ld_matrix))

  expect_equal(result$data$SNP, "rs1")
  expect_equal(result$ld_sign, 1)
})

# --- compute_ld_matrix --------------------------------------------------------
# Integration test: requires local bfile + plink

test_that("compute_ld_matrix works with local reference panel", {
  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
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

  expect_true(is.list(result))
  expect_true(is.matrix(result$ld))
  expect_equal(nrow(result$ld), length(test_snps))
  # Row/col names should be clean rsIDs (no allele suffixes)
  expect_false(any(grepl("_", rownames(result$ld))))
  # alleles carries the parsed-out A1/A2 for every SNP in ld
  expect_true(is.data.frame(result$alleles))
  expect_setequal(result$alleles$SNP, rownames(result$ld))
  expect_true(all(c("ld_a1", "ld_a2") %in% names(result$alleles)))
})

# --- clump_instruments --------------------------------------------------------
# Integration test: requires local bfile + plink

test_that("clump_instruments works with local reference panel", {
  bfile <- sub(
    "\\.bed$",
    "",
    system.file("extdata", "ld_ref.bed", package = "mrpipeline")
  )
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
