# Tests for R/format_single_cell_eqtl.R:
# read_vcf_data() (internal), format_single_cell_onek1k(),
# format_sceqtl_1m_scbloodnl(), format_sceqtl_dice(),
# format_sceqtl_dynamic_cseqtl()

skip_if_not_installed("TwoSampleMR")

# -- Helpers ------------------------------------------------------------------

make_onek1k_df <- function(n = 5, multi_cell = FALSE) {
  df <- data.frame(
    CHR = rep(1L, n),
    POS = as.integer(seq(1e6, by = 1e4, length.out = n)),
    RSID = paste0("rs", seq_len(n)),
    A1 = rep("A", n),
    A2 = rep("G", n),
    A2_FREQ_ONEK1K = seq(0.1, 0.5, length.out = n),
    SPEARMANS_RHO = c(0.15, -0.20, 0.10, -0.25, 0.30)[seq_len(n)],
    P_VALUE = c(1e-6, 1e-7, 1e-8, 1e-6, 1e-9)[seq_len(n)],
    GENE = rep("CD40", n),
    stringsAsFactors = FALSE
  )
  if (multi_cell) {
    df$CELL_ID <- rep(c("cd4nc", "cd8nc"), length.out = n)
  }
  df
}

write_tsv <- function(df, path) {
  # compress = "none": some paths below use a .gz extension (to exercise
  # filename-based cell-type inference), but data.table::fread() only
  # requires the optional R.utils package for genuinely gzip-compressed
  # content, not for the extension alone -- writing plain text keeps these
  # tests independent of that unlisted dependency.
  data.table::fwrite(df, path, sep = "\t", compress = "none")
  path
}

make_scbloodnl_df <- function(n = 5) {
  data.frame(
    SNPName = paste0("rs", seq_len(n) + 100L),
    SNPChr = rep(1L, n),
    SNPChrPos = as.integer(seq(2e6, by = 1e4, length.out = n)),
    CisTrans = rep("cis", n),
    PValue = c(1e-6, 1e-7, 1e-8, 1e-6, 1e-9)[seq_len(n)],
    OverallZScore = c(5.2, -6.1, 4.8, -5.5, 7.0)[seq_len(n)],
    DatasetsNrSamples = rep("100;200;150", n),
    SNPType = rep("A/G", n),
    AlleleAssessed = rep("A", n),
    HGNCName = rep("IL6", n),
    stringsAsFactors = FALSE
  )
}

make_dice_vcf_lines <- function(
  n = 5,
  ref = c("A", "C", "G", "A", "T")[seq_len(n)],
  alt = c("G", "T", "C", "C", "G")[seq_len(n)],
  pval = c(1e-6, 1e-7, 1e-8, 1e-6, 1e-9)[seq_len(n)],
  beta = c(0.1, -0.2, 0.15, -0.1, 0.2)[seq_len(n)],
  gene_symbol = rep("CD40", n)
) {
  info <- sprintf(
    "Gene=ENSG%d;GeneSymbol=%s;Pvalue=%s;Beta=%s",
    seq_len(n),
    gene_symbol,
    pval,
    beta
  )
  info[is.na(gene_symbol)] <- sprintf(
    "Gene=ENSG%d;Pvalue=%s;Beta=%s",
    seq_len(n)[is.na(gene_symbol)],
    pval[is.na(gene_symbol)],
    beta[is.na(gene_symbol)]
  )
  data_lines <- paste(
    "6",
    seq(1e5, by = 1000L, length.out = n),
    paste0("rs", seq_len(n) + 1000L),
    ref,
    alt,
    ".",
    "PASS",
    info,
    sep = "\t"
  )
  c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=6>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    data_lines
  )
}

make_dynamic_df <- function(n = 5) {
  data.frame(
    gene = rep("IL6", n),
    SNP = paste0("rs", seq_len(n) + 2000L),
    beta = c(0.1, -0.2, 0.15, -0.1, 0.2)[seq_len(n)],
    se = c(0.02, 0.03, 0.02, 0.025, 0.03)[seq_len(n)],
    effect_allele = rep("A", n),
    other_allele = rep("G", n),
    pval = c(1e-6, 1e-7, 1e-8, 1e-6, 1e-9)[seq_len(n)],
    chr = rep(1L, n),
    pos = as.integer(seq(3e6, by = 1e4, length.out = n)),
    stringsAsFactors = FALSE
  )
}

# -- read_vcf_data (internal) --------------------------------------------------

test_that("read_vcf_data skips ## lines and parses the #CHROM header", {
  path <- tempfile(fileext = ".vcf")
  writeLines(c("##meta1", "##meta2", "#CHROM\tPOS\tID", "1\t100\trs1"), path)
  df <- read_vcf_data(path)
  expect_equal(names(df), c("CHROM", "POS", "ID"))
  expect_equal(nrow(df), 1L)
  expect_equal(as.character(df$ID[1]), "rs1")
})

test_that("read_vcf_data aborts when no #CHROM header line is found", {
  path <- tempfile(fileext = ".vcf")
  writeLines(c("##meta1", "##meta2"), path)
  expect_error(read_vcf_data(path), "No VCF column header")
})

# -- format_single_cell_onek1k -------------------------------------------------

test_that("format_single_cell_onek1k returns a TwoSampleMR exposure data frame", {
  path <- write_tsv(make_onek1k_df(), tempfile(fileext = ".tsv"))
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(mapping, "cd4nc")
  expect_s3_class(result, "data.frame")
  expect_true(all(
    c("chr.exposure", "pos.exposure", "eaf.exposure") %in% names(result)
  ))
})

test_that("format_single_cell_onek1k builds a GENE___cell_type phenotype label", {
  path <- write_tsv(make_onek1k_df(), tempfile(fileext = ".tsv"))
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(mapping, "cd4nc")
  expect_true(all(result$exposure == "CD40___cd4nc"))
})

test_that("format_single_cell_onek1k filters by CELL_ID on a combined file", {
  path <- write_tsv(
    make_onek1k_df(n = 6, multi_cell = TRUE),
    tempfile(fileext = ".tsv")
  )
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(mapping, "cd4nc", pval_thresh = NULL)
  expect_equal(nrow(result), 3L)
})

test_that("format_single_cell_onek1k aborts when onek1k_mapping is missing required columns", {
  mapping <- data.frame(cell_type = "cd4nc", stringsAsFactors = FALSE)
  expect_error(format_single_cell_onek1k(mapping, "cd4nc"), "missing column")
})

test_that("format_single_cell_onek1k applies pval_thresh filtering", {
  df <- make_onek1k_df()
  df$P_VALUE <- c(1e-8, 1e-2, 1e-8, 1e-2, 1e-8)
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(mapping, "cd4nc", pval_thresh = 1e-5)
  expect_equal(nrow(result), 3L)
})

test_that("format_single_cell_onek1k excludes variants in mhc_region", {
  df <- make_onek1k_df()
  df$CHR <- 6L
  df$POS <- c(20e6, 27e6, 30e6, 33e6, 40e6)
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(
    mapping,
    "cd4nc",
    pval_thresh = NULL,
    mhc_region = data.frame(chr = "6", start = 25e6, end = 34e6)
  )
  expect_equal(nrow(result), 2L)
})

test_that("format_single_cell_onek1k adds samplesize.exposure when sample_size is supplied", {
  path <- write_tsv(make_onek1k_df(), tempfile(fileext = ".tsv"))
  mapping <- data.frame(
    cell_type = "cd4nc",
    path_to_eqtl_file = path,
    stringsAsFactors = FALSE
  )
  result <- format_single_cell_onek1k(mapping, "cd4nc", sample_size = 463528L)
  expect_true(all(result$samplesize.exposure == 463528L))
})

# -- format_sceqtl_1m_scbloodnl ------------------------------------------------

test_that("format_sceqtl_1m_scbloodnl returns a TwoSampleMR exposure data frame", {
  path <- write_tsv(make_scbloodnl_df(), tempfile(fileext = ".tsv"))
  result <- format_sceqtl_1m_scbloodnl(path, cell_type = "cd4nc")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("chr.exposure", "pos.exposure") %in% names(result)))
})

test_that("format_sceqtl_1m_scbloodnl derives beta/se from z-score and total N", {
  path <- write_tsv(make_scbloodnl_df(n = 1), tempfile(fileext = ".tsv"))
  result <- format_sceqtl_1m_scbloodnl(path, cell_type = "cd4nc")
  expected_se <- 1 / sqrt(450)
  expect_equal(result$se.exposure[1], expected_se, tolerance = 1e-10)
  expect_equal(result$beta.exposure[1], 5.2 * expected_se, tolerance = 1e-10)
})

test_that("format_sceqtl_1m_scbloodnl drops rows with ambiguous other allele", {
  df <- make_scbloodnl_df(n = 3)
  df$SNPType <- rep("A/G", 3)
  df$AlleleAssessed <- c("A", "T", "A")
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  result <- format_sceqtl_1m_scbloodnl(
    path,
    cell_type = "cd4nc",
    pval_thresh = NULL
  )
  expect_equal(nrow(result), 2L)
})

test_that("format_sceqtl_1m_scbloodnl cis_only filters to CisTrans == 'cis'", {
  df <- make_scbloodnl_df(n = 4)
  df$CisTrans <- c("cis", "trans", "cis", "trans")
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  result <- format_sceqtl_1m_scbloodnl(
    path,
    cell_type = "cd4nc",
    pval_thresh = NULL
  )
  expect_equal(nrow(result), 2L)
})

test_that("format_sceqtl_1m_scbloodnl infers cell_type from filename", {
  # A .gz-suffixed path requires the (Suggests-only) R.utils package for
  # data.table::fread() to open it at all, regardless of actual content.
  skip_if_not_installed("R.utils")
  path <- file.path(tempdir(), "CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz")
  write_tsv(make_scbloodnl_df(), path)
  result <- format_sceqtl_1m_scbloodnl(path)
  expect_true(all(grepl("___CD4T$", result$exposure)))
})

# -- format_sceqtl_dice ---------------------------------------------------------

test_that("format_sceqtl_dice parses the VCF INFO field and derives se", {
  path <- tempfile(fileext = ".vcf")
  writeLines(make_dice_vcf_lines(), path)
  result <- format_sceqtl_dice(path, cell_type = "cd4nc")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("chr.exposure", "pos.exposure") %in% names(result)))
})

test_that("format_sceqtl_dice drops palindromic SNPs", {
  lines <- make_dice_vcf_lines(
    n = 2,
    ref = c("A", "A"),
    alt = c("T", "G"),
    pval = c(1e-6, 1e-6),
    beta = c(0.1, 0.1)
  )
  path <- tempfile(fileext = ".vcf")
  writeLines(lines, path)
  result <- format_sceqtl_dice(path, cell_type = "cd4nc")
  expect_equal(nrow(result), 1L)
})

test_that("format_sceqtl_dice drops non-SNP (indel) variants", {
  lines <- make_dice_vcf_lines(
    n = 2,
    ref = c("A", "AT"),
    alt = c("G", "A"),
    pval = c(1e-6, 1e-6),
    beta = c(0.1, 0.1)
  )
  path <- tempfile(fileext = ".vcf")
  writeLines(lines, path)
  result <- format_sceqtl_dice(path, cell_type = "cd4nc")
  expect_equal(nrow(result), 1L)
})

test_that("format_sceqtl_dice falls back to Gene when GeneSymbol is absent", {
  lines <- make_dice_vcf_lines(n = 1, gene_symbol = NA_character_)
  path <- tempfile(fileext = ".vcf")
  writeLines(lines, path)
  result <- format_sceqtl_dice(path, cell_type = "cd4nc")
  expect_true(all(result$exposure == "ENSG1___cd4nc"))
})

test_that("format_sceqtl_dice infers cell_type from filename", {
  path <- file.path(tempdir(), "t_cell_cd4_naive.vcf")
  writeLines(make_dice_vcf_lines(), path)
  result <- format_sceqtl_dice(path)
  expect_true(all(grepl("___t_cell_cd4_naive$", result$exposure)))
})

# -- format_sceqtl_dynamic_cseqtl ----------------------------------------------

test_that("format_sceqtl_dynamic_cseqtl passes through pre-formatted columns", {
  path <- write_tsv(make_dynamic_df(), tempfile(fileext = ".tsv"))
  result <- format_sceqtl_dynamic_cseqtl(path, cell_type = "cd4nc")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("chr.exposure", "pos.exposure") %in% names(result)))
  expect_true(all(result$exposure == "IL6___cd4nc"))
})

test_that("format_sceqtl_dynamic_cseqtl applies pval_thresh filtering", {
  df <- make_dynamic_df()
  df$pval <- c(1e-8, 1e-2, 1e-8, 1e-2, 1e-8)
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  result <- format_sceqtl_dynamic_cseqtl(
    path,
    cell_type = "cd4nc",
    pval_thresh = 1e-5
  )
  expect_equal(nrow(result), 3L)
})

test_that("format_sceqtl_dynamic_cseqtl excludes variants in mhc_region", {
  df <- make_dynamic_df()
  df$chr <- 6L
  df$pos <- c(20e6, 27e6, 30e6, 33e6, 40e6)
  path <- write_tsv(df, tempfile(fileext = ".tsv"))
  result <- format_sceqtl_dynamic_cseqtl(
    path,
    cell_type = "cd4nc",
    pval_thresh = NULL,
    mhc_region = data.frame(chr = "6", start = 25e6, end = 34e6)
  )
  expect_equal(nrow(result), 2L)
})

test_that("format_sceqtl_dynamic_cseqtl infers cell_type from filename", {
  skip_if_not_installed("R.utils")
  path <- file.path(tempdir(), "CD4T_500kb_combined.MR.tsv.gz")
  write_tsv(make_dynamic_df(), path)
  result <- format_sceqtl_dynamic_cseqtl(path)
  expect_true(all(grepl("___CD4T$", result$exposure)))
})
