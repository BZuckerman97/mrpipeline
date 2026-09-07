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

# --- check_allele_orientation / last_allele_check ----------------------------
# Fixtures: make_allele_fixture() / make_allele_gwas_fixture() in
# helper-allele-fixture.R. Under the bug, every non-palindromic outcome beta
# is sign-inverted and eaf.outcome becomes the complement of eaf.exposure,
# while palindromic SNPs look untouched -- see the "palindromic" test below.

test_that("allele check errors by default on a swapped-allele outcome", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()

  expect_error(
    suppressMessages(harmonise_and_filter(f$exposure, f$bug)),
    class = "mrpipeline_allele_check_error"
  )
  # Message is actionable: counts, an offender, and the fix
  expect_error(
    suppressMessages(harmonise_and_filter(f$exposure, f$bug)),
    "between \"exp\" and \"out\": 10/10"
  )
  expect_error(
    suppressMessages(harmonise_and_filter(f$exposure, f$bug)),
    "rs[0-9]+ \\([ACGT]\\): eaf.exposure = "
  )
  expect_error(
    suppressMessages(harmonise_and_filter(f$exposure, f$bug)),
    "col_map"
  )

  rec <- last_allele_check()
  expect_equal(rec$status, "fail")
  expect_equal(rec$n, 10L)
  expect_equal(rec$n_complementary, 10L)
  expect_equal(rec$prop, 1)
  expect_equal(rec$exposure, "exp")
  expect_equal(rec$outcome, "out")
  expect_equal(rec$id.exposure, "exp1")
  expect_equal(rec$id.outcome, "out1")
  expect_true(is.na(rec$n_sampled))
  expect_equal(nrow(rec$variants), 12L)
  expect_equal(sum(rec$variants$palindromic), 2L)
  expect_true(all(rec$variants$informative == !rec$variants$palindromic))
  expect_true(all(is.na(rec$variants$complementary[rec$variants$palindromic])))
})

test_that("palindromic betas are unchanged by the bug, so the check must ignore them", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()

  ok <- suppressMessages(harmonise_and_filter(
    f$exposure,
    f$ok,
    allele_check = "none"
  ))
  bug <- suppressMessages(harmonise_and_filter(
    f$exposure,
    f$bug,
    allele_check = "none"
  ))
  bug <- bug[match(ok$SNP, bug$SNP), ]
  pal <- ok$palindromic

  expect_equal(sum(pal), 2L)
  # Palindromic: strand resolved from the mis-assigned frequency, errors cancel
  expect_equal(bug$beta.outcome[pal], ok$beta.outcome[pal])
  expect_equal(bug$eaf.outcome[pal], ok$eaf.outcome[pal])
  # Non-palindromic: every beta sign-inverted, eaf.outcome the complement
  expect_equal(bug$beta.outcome[!pal], -ok$beta.outcome[!pal])
  expect_equal(bug$eaf.outcome[!pal], 1 - ok$eaf.outcome[!pal])
})

test_that("allele check passes on genuine EAF scatter and records a pass", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture(noise_sd = 0.08)

  expect_no_condition(
    res <- suppressMessages(harmonise_and_filter(f$exposure, f$ok))
  )
  expect_equal(nrow(res), 12L)

  rec <- last_allele_check()
  expect_equal(rec$status, "pass")
  expect_equal(rec$n, 10L)
  expect_lt(rec$prop, 0.7)
})

test_that("allele_check = 'none' reproduces plain harmonisation exactly", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()

  expect_no_condition(
    res <- suppressMessages(
      harmonise_and_filter(f$exposure, f$bug, allele_check = "none")
    )
  )
  plain <- suppressMessages(TwoSampleMR::harmonise_data(f$exposure, f$bug)) |>
    dplyr::filter(.data$mr_keep == TRUE) |>
    dplyr::filter(!duplicated(.data$SNP))
  expect_identical(res, plain)

  # The record is still stored so the verdict can be inspected after opting out
  rec <- last_allele_check()
  expect_equal(rec$status, "fail")
  expect_equal(rec$allele_check, "none")
})

test_that("allele_check = 'warn' warns with its class and still returns data", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()

  expect_warning(
    res <- suppressMessages(
      harmonise_and_filter(f$exposure, f$bug, allele_check = "warn")
    ),
    class = "mrpipeline_allele_check_warning"
  )
  expect_equal(nrow(res), 12L)
  expect_equal(last_allele_check()$status, "fail")
})

test_that("allele check is skipped below 10 informative SNPs", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()
  # Drop one non-palindromic SNP: 9 informative + 2 palindromic
  keep <- f$exposure$SNP != "rs10"

  expect_no_condition(
    res <- suppressMessages(
      harmonise_and_filter(f$exposure[keep, ], f$bug[keep, ])
    )
  )
  expect_equal(nrow(res), 11L)

  rec <- last_allele_check()
  expect_equal(rec$status, "skipped")
  expect_equal(rec$n, 9L)
  expect_equal(rec$n_complementary, 9L)
})

test_that("verbose = TRUE reports pass and skipped verdicts", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()

  expect_message(
    suppressMessages(
      harmonise_and_filter(f$exposure, f$ok, verbose = TRUE),
      classes = "simpleMessage"
    ),
    "check passed"
  )
  keep <- f$exposure$SNP != "rs10"
  expect_message(
    suppressMessages(
      harmonise_and_filter(f$exposure[keep, ], f$ok[keep, ], verbose = TRUE),
      classes = "simpleMessage"
    ),
    "check skipped"
  )
})

test_that("last_allele_check returns NULL before any check and the record after", {
  skip_if_not_installed("TwoSampleMR")
  old <- the$last_allele_check
  on.exit(the$last_allele_check <- old)
  the$last_allele_check <- NULL

  expect_message(
    expect_null(last_allele_check()),
    "No allele orientation check"
  )

  f <- make_allele_fixture()
  suppressMessages(harmonise_and_filter(f$exposure, f$ok))
  expect_equal(last_allele_check()$status, "pass")

  suppressMessages(harmonise_and_filter(
    f$exposure,
    f$bug,
    allele_check = "none"
  ))
  expect_equal(last_allele_check()$status, "fail")
})

test_that("allele check is skipped when allele frequencies are all missing", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_fixture()
  # The shape TwoSampleMR::format_data() produces when a file has no EAF:
  # columns present, all NA
  f$exposure$eaf.exposure <- NA_real_
  f$bug$eaf.outcome <- NA_real_

  expect_no_condition(
    res <- suppressMessages(harmonise_and_filter(f$exposure, f$bug))
  )
  expect_equal(nrow(res), 10L) # palindromic SNPs dropped without EAF

  rec <- last_allele_check()
  expect_equal(rec$status, "skipped")
  expect_equal(rec$n, 0L)
  expect_equal(nrow(rec$variants), 0L)
})

test_that("check_allele_orientation handles hand-built frames directly", {
  # No TwoSampleMR needed: operates on any frame with the harmonised columns
  harmonised <- data.frame(
    SNP = c(paste0("rs", 1:12), "rs1"),
    eaf.exposure = c(rep(0.2, 12), 0.2),
    eaf.outcome = c(rep(0.8, 10), 0.5, 0.8, 0.8),
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    palindromic = FALSE,
    remove = c(rep(FALSE, 11), TRUE, FALSE),
    id.exposure = "e",
    id.outcome = "o",
    stringsAsFactors = FALSE
  )

  expect_error(
    check_allele_orientation(harmonised),
    class = "mrpipeline_allele_check_error"
  )
  rec <- last_allele_check()
  # rs12 excluded (remove == TRUE), duplicate rs1 counted once, rs11 is a tie
  expect_equal(nrow(rec$variants), 12L)
  expect_equal(rec$n, 11L)
  expect_equal(rec$n_complementary, 10L)
  expect_false(rec$variants$complementary[rec$variants$SNP == "rs11"])
  expect_true(is.na(rec$variants$complementary[rec$variants$SNP == "rs12"]))

  # Missing columns -> skipped, never an error
  expect_no_condition(check_allele_orientation(data.frame(SNP = "rs1")))
  expect_equal(last_allele_check()$status, "skipped")
  expect_no_condition(check_allele_orientation(data.frame()))
  expect_equal(last_allele_check()$status, "skipped")

  # Invalid mode is rejected
  expect_error(
    check_allele_orientation(harmonised, allele_check = "bogus"),
    "allele_check"
  )
})

test_that("check_allele_orientation_gwas checks instruments plus shared SNPs", {
  skip_if_not_installed("TwoSampleMR")
  f <- make_allele_gwas_fixture()

  expect_error(
    check_allele_orientation_gwas(f$exposure, f$bug, f$instruments),
    class = "mrpipeline_allele_check_error"
  )
  rec <- last_allele_check()
  expect_equal(rec$status, "fail")
  expect_equal(rec$n_sampled, 37L)
  expect_equal(rec$n, 40L)
  expect_equal(rec$n_complementary, 40L)

  expect_no_condition(
    check_allele_orientation_gwas(f$exposure, f$ok, f$instruments)
  )
  expect_equal(last_allele_check()$status, "pass")
  expect_equal(last_allele_check()$n_complementary, 0L)

  # n_sample caps the number of non-instrument SNPs added
  expect_no_condition(
    check_allele_orientation_gwas(
      f$exposure,
      f$bug,
      f$instruments,
      n_sample = 5L
    )
  )
  rec <- last_allele_check()
  expect_equal(rec$n_sampled, 5L)
  expect_equal(rec$n, 8L)
  expect_equal(rec$status, "skipped")

  # Outcome without EAF, or with nothing in common -> skipped
  out_no_eaf <- f$bug
  out_no_eaf$eaf <- NULL
  expect_no_condition(
    check_allele_orientation_gwas(f$exposure, out_no_eaf, f$instruments)
  )
  expect_equal(last_allele_check()$status, "skipped")

  out_disjoint <- f$bug
  out_disjoint$rsids <- paste0("rs", 100 + seq_len(nrow(out_disjoint)))
  expect_no_condition(
    check_allele_orientation_gwas(f$exposure, out_disjoint, f$instruments)
  )
  expect_equal(last_allele_check()$status, "skipped")
  expect_no_condition(
    check_allele_orientation_gwas(f$exposure, data.frame(), f$instruments)
  )
  expect_equal(last_allele_check()$status, "skipped")
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

# --- compute_ld_to_index ------------------------------------------------------
# Integration test: requires local bfile + plink

test_that("compute_ld_to_index matches compute_ld_matrix's squared, signed r", {
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
  test_snps <- head(bim$V2, 5)
  index_snp <- test_snps[1]

  result <- compute_ld_to_index(test_snps, index_snp, bfile)

  expect_true(is.data.frame(result))
  expect_setequal(names(result), c("SNP", "r2"))
  expect_setequal(result$SNP, test_snps)

  m <- compute_ld_matrix(test_snps, bfile)
  r2_direct <- as.numeric(m$ld[, index_snp])^2
  expect_equal(
    result$r2[match(rownames(m$ld), result$SNP)],
    r2_direct
  )

  # LD of the index SNP against itself is always 1.
  expect_equal(result$r2[result$SNP == index_snp], 1)
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
