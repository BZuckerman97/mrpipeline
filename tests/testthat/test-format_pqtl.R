# Tests for format_pqtl_decode() and format_pqtl_ukbppp()

skip_if_not_installed("TwoSampleMR")

# -- Helpers ------------------------------------------------------------------

make_decode_gwas <- function(n = 5, chrom = paste0("chr", seq_len(n))) {
  data.frame(
    Name = paste0("rs", seq_len(n)),
    rsids = paste0("rs", seq_len(n)),
    Beta = c(0.10, -0.20, 0.05, -0.15, 0.08)[seq_len(n)],
    SE = c(0.02, 0.03, 0.02, 0.04, 0.01)[seq_len(n)],
    Pval = c(1e-8, 1e-5, 1e-7, 1e-4, 1e-6)[seq_len(n)],
    effectAllele = c("A", "G", "C", "T", "A")[seq_len(n)],
    otherAllele = c("G", "A", "T", "C", "G")[seq_len(n)],
    Pos = as.integer(seq(1e6, by = 1e5, length.out = n)),
    Chrom = chrom,
    N = rep(35559L, n),
    stringsAsFactors = FALSE
  )
}

make_decode_included <- function(gwas) {
  data.frame(
    Name = gwas$Name,
    effectAlleleFreq = seq(0.1, 0.5, length.out = nrow(gwas)),
    stringsAsFactors = FALSE
  )
}

make_ukbppp_gwas <- function(n = 5, chrom = as.character(seq_len(n))) {
  ids <- paste(chrom, seq(1e6, by = 1e5, length.out = n), "A", "G", sep = ":")
  data.frame(
    ID = ids,
    BETA = c(0.10, -0.20, 0.05, -0.15, 0.08)[seq_len(n)],
    SE = c(0.02, 0.03, 0.02, 0.04, 0.01)[seq_len(n)],
    A1FREQ = seq(0.1, 0.5, length.out = n),
    ALLELE1 = rep("A", n),
    ALLELE0 = rep("G", n),
    LOG10P = c(8, 5, 7, 4, 6)[seq_len(n)],
    CHROM = chrom,
    N = rep(50000L, n),
    GENPOS = as.integer(seq(1e6, by = 1e5, length.out = n)),
    stringsAsFactors = FALSE
  )
}

make_ukbppp_rsid <- function(gwas) {
  data.frame(
    ID = gwas$ID,
    rsid = paste0("rs", seq_len(nrow(gwas))),
    ALT = rep("A", nrow(gwas)),
    REF = rep("G", nrow(gwas)),
    POS38 = gwas$GENPOS,
    stringsAsFactors = FALSE
  )
}

# -- format_pqtl_decode -------------------------------------------------------

test_that("format_pqtl_decode returns a list with an exposure element", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40")
  expect_type(result, "list")
  expect_named(result, "exposure")
  expect_s3_class(result$exposure, "data.frame")
})

test_that("format_pqtl_decode sets exposure column to pqtl_assay", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "IL18")
  expect_true(all(result$exposure$exposure == "IL18"))
})

test_that("format_pqtl_decode strips chr prefix from Chrom", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40")
  chr_col <- result$exposure$chr.exposure
  expect_false(any(grepl("^chr", chr_col), na.rm = TRUE))
})

test_that("format_pqtl_decode converts chromosome 23 to X", {
  gwas <- make_decode_gwas(n = 3, chrom = c("chr1", "chr23", "chr22"))
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "TEST")
  chr_col <- result$exposure$chr.exposure
  expect_false("23" %in% chr_col)
  expect_true("X" %in% chr_col)
})

test_that("format_pqtl_decode retains all input SNPs via inner join on Name", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40")
  expect_equal(nrow(result$exposure), nrow(gwas))
})

test_that("format_pqtl_decode drops SNPs absent from included_variants", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)[1:3, ]
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40")
  expect_equal(nrow(result$exposure), 3L)
})

# -- format_pqtl_ukbppp -------------------------------------------------------

test_that("format_pqtl_ukbppp returns a data frame", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "IL6")
  expect_s3_class(result, "data.frame")
})

test_that("format_pqtl_ukbppp sets exposure column to pqtl_assay", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "CD40")
  expect_true(all(result$exposure == "CD40"))
})

test_that("format_pqtl_ukbppp converts LOG10P to p-values in (0, 1]", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "IL6")
  pvals <- result$pval.exposure
  expect_true(all(pvals > 0, na.rm = TRUE))
  expect_true(all(pvals <= 1, na.rm = TRUE))
})

test_that("format_pqtl_ukbppp converts chromosome 23 to X", {
  gwas <- make_ukbppp_gwas(n = 3, chrom = c("1", "23", "2"))
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "TEST")
  expect_false("23" %in% result$chr.exposure)
  expect_true("X" %in% result$chr.exposure)
})

test_that("format_pqtl_ukbppp joins on ID column when present in both inputs", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "IL6")
  expect_equal(nrow(result), nrow(gwas))
  expect_true("SNP" %in% names(result))
})

test_that("format_pqtl_ukbppp drops SNPs absent from rsid file", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)[1:3, ]
  result <- format_pqtl_ukbppp(gwas, rsids, pqtl_assay = "IL6")
  expect_equal(nrow(result), 3L)
})
