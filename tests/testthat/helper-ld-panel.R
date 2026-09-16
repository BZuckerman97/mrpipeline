# Fixtures for LD-correction tests that need instruments in genuine LD.
#
# The bundled inst/extdata/ld_ref panel has random genotypes -- the largest
# |r| between the CD40 instruments is 0.023 -- so on it corrected and
# uncorrected results all but coincide. That is why nothing noticed the
# diagnostics were uncorrected (issue #31). make_ld_panel() writes a small
# PLINK fileset with a chosen pair of SNPs in strong LD, so an end-to-end
# run_mr(ld_correct = TRUE) has something to correct for.

# Five synthetic, non-palindromic instruments on chr 1, as an exposure /
# outcome pair in the shapes run_mr() takes.
synthetic_five_snps <- function() {
  exposure <- data.frame(
    SNP = paste0("rs", 1:5),
    beta.exposure = c(0.50, 0.30, 0.40, 0.35, 0.45),
    se.exposure = c(0.10, 0.10, 0.10, 0.10, 0.10),
    effect_allele.exposure = c("A", "G", "C", "T", "G"),
    other_allele.exposure = c("G", "T", "A", "C", "A"),
    pval.exposure = c(1e-10, 1e-8, 1e-9, 1e-8, 1e-9),
    eaf.exposure = c(0.30, 0.40, 0.50, 0.35, 0.45),
    exposure = "test_exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    chr.exposure = "1",
    pos.exposure = c(1000, 2000, 3000, 4000, 5000),
    samplesize.exposure = 10000,
    stringsAsFactors = FALSE
  )
  outcome <- data.frame(
    rsids = paste0("rs", 1:5),
    beta = c(0.10, 0.05, 0.08, 0.02, 0.09),
    se = c(0.05, 0.03, 0.04, 0.03, 0.04),
    pval = c(0.01, 0.10, 0.05, 0.50, 0.03),
    eaf = c(0.30, 0.40, 0.50, 0.35, 0.45),
    effect_allele = c("A", "G", "C", "T", "G"),
    other_allele = c("G", "T", "A", "C", "A"),
    chr = "1",
    pos = c(1000, 2000, 3000, 4000, 5000),
    n = 5000,
    phenotype = "test_out",
    stringsAsFactors = FALSE
  )
  list(exposure = exposure, outcome = outcome)
}

# Write a PLINK bfile for `exposure`'s SNPs to tempdir() with SNPs
# `ld_pair[1]` and `ld_pair[2]` in LD: the second is a copy of the first on
# a fraction `copy_prob` of haplotypes (copy_prob = 1 makes them identical,
# i.e. a singular LD matrix). Genotypes go through a text .ped/.map and
# `plink --make-bed`, so no .bed bytes are packed by hand. Skips if PLINK is
# unavailable. Returns the bfile prefix.
make_ld_panel <- function(
  exposure,
  ld_pair = c(1, 2),
  copy_prob = 0.85,
  n_individuals = 200,
  seed = 1
) {
  plink <- tryCatch(
    genetics.binaRies::get_plink_binary(),
    error = function(e) NULL
  )
  testthat::skip_if(
    is.null(plink) || !file.exists(plink),
    "PLINK binary not available"
  )

  set.seed(seed)
  m <- nrow(exposure)
  n_hap <- 2 * n_individuals
  # 1 = effect allele, 0 = other allele, per haplotype
  hap <- matrix(stats::rbinom(n_hap * m, 1, 0.4), nrow = n_hap, ncol = m)
  copy <- stats::rbinom(n_hap, 1, copy_prob) == 1
  hap[copy, ld_pair[2]] <- hap[copy, ld_pair[1]]

  allele <- function(j, h) {
    ifelse(
      h == 1,
      exposure$effect_allele.exposure[j],
      exposure$other_allele.exposure[j]
    )
  }
  geno <- vapply(
    seq_len(m),
    function(j) {
      a <- allele(j, hap[seq(1, n_hap, by = 2), j])
      b <- allele(j, hap[seq(2, n_hap, by = 2), j])
      paste(a, b)
    },
    character(n_individuals)
  )

  prefix <- file.path(
    tempdir(),
    sprintf("ld_panel_%s_%d", paste(ld_pair, collapse = "_"), seed)
  )
  ped <- data.frame(
    fid = "FAM",
    iid = seq_len(n_individuals),
    pat = 0,
    mat = 0,
    sex = 1,
    pheno = -9,
    geno,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    ped,
    paste0(prefix, ".ped"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  map <- data.frame(
    chr = exposure$chr.exposure,
    snp = exposure$SNP,
    cm = 0,
    pos = exposure$pos.exposure
  )
  utils::write.table(
    map,
    paste0(prefix, ".map"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  system2(
    plink,
    c("--file", prefix, "--make-bed", "--out", prefix),
    stdout = FALSE,
    stderr = FALSE
  )
  testthat::skip_if_not(
    file.exists(paste0(prefix, ".bed")),
    "PLINK failed to build the LD test panel"
  )
  prefix
}

# run_mr() with PLINK's console output captured and warnings collected, so
# tests can assert on them without nested expect_warning().
run_quiet <- function(exposure, outcome, ...) {
  warnings <- character()
  result <- NULL
  invisible(utils::capture.output(
    result <- withCallingHandlers(
      suppressMessages(run_mr(
        exposure = exposure,
        exposure_id = "test_exp",
        outcome = outcome,
        outcome_id = "test_out",
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
