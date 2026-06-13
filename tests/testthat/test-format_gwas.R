# Tests for format_gwas()

# -- Helper -------------------------------------------------------------------

make_gwas <- function(...) {
  defaults <- list(
    rsids        = paste0("rs", 1:5),
    chr          = as.character(1:5),
    pos          = as.integer(seq(1e6, 5e6, by = 1e6)),
    beta         = c(0.10, -0.20, 0.30, -0.10, 0.05),
    se           = c(0.02,  0.03, 0.02,  0.04, 0.01),
    eaf          = c(0.30,  0.40, 0.50,  0.20, 0.45),
    pval         = c(1e-8,  1e-5, 1e-6,  1e-4, 1e-7),
    effect_allele = c("A", "G", "C", "T", "A"),
    other_allele  = c("G", "A", "T", "C", "G"),
    n            = rep(10000L, 5)
  )
  overrides <- list(...)
  merged <- c(defaults[setdiff(names(defaults), names(overrides))], overrides)
  do.call(data.frame, c(merged, stringsAsFactors = FALSE))
}

# -- Data frame input ---------------------------------------------------------

test_that("format_gwas accepts a data frame as path", {
  result <- format_gwas(make_gwas(), phenotype_id = "TEST")
  expect_s3_class(result, "data.frame")
  expect_equal(result$phenotype[1], "TEST")
})

test_that("format_gwas errors when path is not a file or data frame", {
  expect_error(format_gwas(123L, phenotype_id = "TEST"), "file path or data frame")
})

test_that("format_gwas errors when file path does not exist", {
  expect_error(
    format_gwas("/nonexistent/path/file.tsv", phenotype_id = "TEST"),
    "not found"
  )
})

# -- Column alias resolution --------------------------------------------------

test_that("format_gwas renames SNP to rsids", {
  df <- make_gwas()
  names(df)[names(df) == "rsids"] <- "SNP"
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("rsids" %in% names(result))
})

test_that("format_gwas renames P to pval", {
  df <- make_gwas()
  names(df)[names(df) == "pval"] <- "P"
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("pval" %in% names(result))
})

test_that("format_gwas renames BETA to beta", {
  df <- make_gwas()
  names(df)[names(df) == "beta"] <- "BETA"
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("beta" %in% names(result))
})

test_that("format_gwas renames SE to se", {
  df <- make_gwas()
  names(df)[names(df) == "se"] <- "SE"
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("se" %in% names(result))
})

test_that("format_gwas renames ALLELE1/ALLELE0 to effect/other_allele", {
  df <- make_gwas()
  names(df)[names(df) == "effect_allele"] <- "ALLELE1"
  names(df)[names(df) == "other_allele"]  <- "ALLELE0"
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("effect_allele" %in% names(result))
  expect_true("other_allele"  %in% names(result))
})

test_that("format_gwas col_map takes precedence over built-in aliases", {
  df <- make_gwas()
  names(df)[names(df) == "pval"] <- "PVALUE"
  result <- format_gwas(df, phenotype_id = "TEST", col_map = list(pval = "PVALUE"))
  expect_true("pval" %in% names(result))
})

# -- Allele uppercasing -------------------------------------------------------

test_that("format_gwas uppercases allele columns", {
  df <- make_gwas(
    effect_allele = c("a", "g", "c", "t", "a"),
    other_allele  = c("g", "a", "t", "c", "g")
  )
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_equal(result$effect_allele, c("A", "G", "C", "T", "A"))
  expect_equal(result$other_allele,  c("G", "A", "T", "C", "G"))
})

# -- Chromosome prefix stripping ----------------------------------------------

test_that("format_gwas strips chr prefix from chromosome column", {
  df <- make_gwas(chr = paste0("chr", 1:5))
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_false(any(grepl("^chr", result$chr)))
  expect_equal(result$chr, as.character(1:5))
})

# -- log10_pval back-transform ------------------------------------------------

test_that("format_gwas back-transforms log10_pval = TRUE", {
  df <- make_gwas(pval = -log10(c(1e-8, 1e-5, 1e-6, 1e-4, 1e-7)))
  result <- format_gwas(df, phenotype_id = "TEST", log10_pval = TRUE)
  expect_equal(result$pval, c(1e-8, 1e-5, 1e-6, 1e-4, 1e-7), tolerance = 1e-12)
})

test_that("format_gwas auto-detects and back-transforms neg_log_10_p_value column", {
  df <- make_gwas()
  expected <- df$pval
  names(df)[names(df) == "pval"] <- "neg_log_10_p_value"
  df$neg_log_10_p_value <- -log10(expected)
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_equal(result$pval, expected, tolerance = 1e-12)
})

# -- flip_beta ----------------------------------------------------------------

test_that("format_gwas negates beta when flip_beta = TRUE", {
  df    <- make_gwas()
  orig  <- df$beta
  result <- format_gwas(df, phenotype_id = "TEST", flip_beta = TRUE)
  expect_equal(result$beta, -orig)
})

# -- OR -> beta/se derivation -------------------------------------------------

test_that("format_gwas derives beta = log(OR) when beta column is absent", {
  df     <- make_gwas()
  df$or  <- exp(df$beta)
  df     <- df[, setdiff(names(df), "beta")]
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("beta" %in% names(result))
  expect_equal(result$beta, log(df$or), tolerance = 1e-10)
})

test_that("format_gwas derives se via Z-score when OR present and se absent", {
  df     <- make_gwas()
  df$or  <- exp(df$beta)
  df     <- df[, setdiff(names(df), c("beta", "se"))]
  result <- format_gwas(df, phenotype_id = "TEST")
  expect_true("se" %in% names(result))
  expect_true(all(result$se > 0))
})

test_that("format_gwas errors when OR present but pval absent for SE derivation", {
  df     <- make_gwas()
  df$or  <- exp(df$beta)
  df     <- df[, setdiff(names(df), c("beta", "pval"))]
  expect_error(format_gwas(df, phenotype_id = "TEST"), "cannot derive")
})

# -- marker_col parsing -------------------------------------------------------

test_that("format_gwas parses chr and pos from marker_col", {
  df <- make_gwas()
  df$MarkerName <- paste(df$chr, df$pos, df$effect_allele, sep = ":")
  df <- df[, setdiff(names(df), c("chr", "pos"))]
  result <- format_gwas(df, phenotype_id = "TEST", marker_col = "MarkerName")
  expect_equal(result$chr, as.character(1:5))
  expect_equal(result$pos, as.integer(seq(1e6, 5e6, by = 1e6)))
})

test_that("format_gwas errors when marker_col is not a column in the data", {
  expect_error(
    format_gwas(make_gwas(), phenotype_id = "TEST", marker_col = "NoSuchCol"),
    "not found"
  )
})

# -- Explicit n argument ------------------------------------------------------

test_that("format_gwas adds n column when absent and n argument supplied", {
  df <- make_gwas()[, setdiff(names(make_gwas()), "n")]
  result <- format_gwas(df, phenotype_id = "TEST", n = 50000L)
  expect_true("n" %in% names(result))
  expect_true(all(result$n == 50000L))
})

test_that("format_gwas does not overwrite an existing n column", {
  df     <- make_gwas()
  result <- format_gwas(df, phenotype_id = "TEST", n = 99999L)
  expect_equal(unique(result$n), 10000L)
})

# -- Validation errors --------------------------------------------------------

test_that("format_gwas errors on missing required columns", {
  df <- make_gwas()[, setdiff(names(make_gwas()), c("beta", "se"))]
  expect_error(format_gwas(df, phenotype_id = "TEST"), "required column")
})

test_that("format_gwas errors when rsids absent and no bim_path supplied", {
  df <- make_gwas()[, setdiff(names(make_gwas()), "rsids")]
  expect_error(format_gwas(df, phenotype_id = "TEST"), "rsID")
})

# -- type = "exposure" --------------------------------------------------------

test_that("format_gwas type = exposure returns TwoSampleMR-formatted data frame", {
  skip_if_not_installed("TwoSampleMR")
  result <- format_gwas(make_gwas(), phenotype_id = "TEST", type = "exposure")
  expect_s3_class(result, "data.frame")
  expect_true("SNP" %in% names(result))
  expect_true("beta.exposure" %in% names(result))
  expect_true(all(result$exposure == "TEST"))
})
