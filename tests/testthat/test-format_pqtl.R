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

test_that("format_pqtl_decode outcome type returns a plain data frame with format_gwas(outcome) schema", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40", type = "outcome")
  expect_s3_class(result, "data.frame")
  # Not wrapped in list(exposure = ...) like the type = "exposure" return.
  expect_false("exposure" %in% names(result))
  expect_true(all(
    c(
      "rsids",
      "chr",
      "pos",
      "beta",
      "se",
      "eaf",
      "pval",
      "n",
      "effect_allele",
      "other_allele",
      "phenotype"
    ) %in%
      names(result)
  ))
  expect_equal(nrow(result), nrow(gwas))
})

test_that("format_pqtl_decode outcome type carries N through as n when present", {
  gwas <- make_decode_gwas()
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40", type = "outcome")
  expect_true(all(result$n == 35559L))
})

test_that("format_pqtl_decode outcome type defaults n to NA when N is absent", {
  gwas <- make_decode_gwas()
  gwas$N <- NULL
  inc <- make_decode_included(gwas)
  result <- format_pqtl_decode(gwas, inc, pqtl_assay = "CD40", type = "outcome")
  expect_true(all(is.na(result$n)))
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

test_that("format_pqtl_ukbppp outcome type returns a data frame with format_gwas(outcome) schema", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(
    gwas,
    rsids,
    pqtl_assay = "IL6",
    type = "outcome"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(
    c(
      "rsids",
      "chr",
      "pos",
      "beta",
      "se",
      "eaf",
      "pval",
      "n",
      "effect_allele",
      "other_allele",
      "phenotype"
    ) %in%
      names(result)
  ))
  expect_equal(nrow(result), nrow(gwas))
})

test_that("format_pqtl_ukbppp outcome type carries N through as n when present", {
  gwas <- make_ukbppp_gwas()
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(
    gwas,
    rsids,
    pqtl_assay = "IL6",
    type = "outcome"
  )
  expect_true(all(result$n == 50000L))
})

test_that("format_pqtl_ukbppp outcome type defaults n to NA when N is absent, rather than erroring", {
  # Regression test: N was previously renamed unconditionally via
  # dplyr::rename(n = dplyr::all_of("N")) in preprocessing shared by both
  # exposure and outcome, which hard-errored on a missing N column for
  # *both* types -- silently defeating the outcome branch's own
  # `if (!"n" %in% names(ukbppp)) ukbppp$n <- NA_integer_` fallback, which
  # could never be reached.
  gwas <- make_ukbppp_gwas()
  gwas$N <- NULL
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(
    gwas,
    rsids,
    pqtl_assay = "IL6",
    type = "outcome"
  )
  expect_true(all(is.na(result$n)))
})

test_that("format_pqtl_ukbppp exposure type still works when N is absent", {
  gwas <- make_ukbppp_gwas()
  gwas$N <- NULL
  rsids <- make_ukbppp_rsid(gwas)
  result <- format_pqtl_ukbppp(
    gwas,
    rsids,
    pqtl_assay = "IL6",
    type = "exposure"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(gwas))
})

make_decode_linker <- function() {
  data.frame(
    seqID = c("SEQ1", "SEQ2"),
    identifier = c("file1.txt.gz", "file2.txt.gz"),
    stringsAsFactors = FALSE
  )
}

test_that("decode_pqtl_file_name returns the decode path and id for a matching seqID", {
  linker <- make_decode_linker()
  result <- decode_pqtl_file_name("SEQ1", linker, "some_dir")
  expect_equal(result$decode, file.path("some_dir", "file1.txt.gz"))
  expect_equal(result$id, "SEQ1")
})

test_that("decode_pqtl_file_name errors when unique_id matches zero or multiple rows", {
  linker <- make_decode_linker()
  expect_error(decode_pqtl_file_name("MISSING", linker, "some_dir"))

  linker_dup <- rbind(linker, linker[1, ])
  expect_error(decode_pqtl_file_name("SEQ1", linker_dup, "some_dir"))
})

make_ukbppp_linker_row <- function(chr = "7") {
  data.frame(
    Code = "syn51468683",
    Assay = "IL6",
    chr = chr,
    UKBPPP_ProteinID = "IL6:P05231:OID20101:v1",
    Panel = "Cardiometabolic",
    Docname = "IL6_P05231_OID20101_v1_Cardiometabolic.tar",
    stringsAsFactors = FALSE
  )
}

test_that("ukbppp_pqtl_file_name uses the protein's own chromosome by default", {
  linker <- make_ukbppp_linker_row(chr = "7")
  paths <- ukbppp_pqtl_file_name(
    synapse_id = "syn51468683",
    olink_linker_file = linker,
    olink_dir = "olink",
    olink_rsid_dir = "olink_rsid"
  )
  expect_match(paths$ukbppp, "discovery_chr7_")
  expect_match(paths$ukbppp_rsid, "_chr7_patched_v2")
})

test_that("ukbppp_pqtl_file_name overrides the chromosome when region_chr is supplied", {
  linker <- make_ukbppp_linker_row(chr = "7")
  paths <- ukbppp_pqtl_file_name(
    synapse_id = "syn51468683",
    olink_linker_file = linker,
    olink_dir = "olink",
    olink_rsid_dir = "olink_rsid",
    region_chr = "1"
  )
  expect_match(paths$ukbppp, "discovery_chr1_")
  expect_match(paths$ukbppp_rsid, "_chr1_patched_v2")
})
