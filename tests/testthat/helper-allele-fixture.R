# Shared fixtures for the allele orientation check (GitHub issue #18).
#
# An outcome file whose A1/A2 columns mean REF/ALT is read by format_gwas()
# with the wrong allele as the effect allele, while beta and eaf stay
# oriented to the true effect allele. `ok` is the correctly labelled
# outcome; `bug` is the same data with the effect/other allele columns
# swapped and beta/eaf untouched -- exactly what the mis-read file produces.

# Harmonise with TwoSampleMR's progress messages silenced, returning the
# analysed (filtered, deduplicated) frame. `harmonise_and_filter()` returns
# list(data =, raw =); tests that need the unfiltered frame call it directly.
hf <- function(...) suppressMessages(harmonise_and_filter(...)$data)

# 12 SNPs in TwoSampleMR outcome format: 10 non-palindromic, plus rs11 (C/G,
# eaf 0.05) and rs12 (A/T, eaf 0.20), palindromic with unambiguous
# frequencies. `noise_sd` scales a fixed +/- pattern added to the outcome EAF
# (0.02 = same population; 0.08 = cross-ancestry-like scatter) -- fixed rather
# than random so the fixture is reproducible without touching the RNG.
make_allele_fixture <- function(noise_sd = 0.02) {
  snps <- paste0("rs", 1:12)
  ea <- c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "C", "A")
  oa <- c("G", "T", "A", "C", "C", "A", "T", "G", "G", "T", "G", "T")
  eaf <- c(
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.65,
    0.75,
    0.80,
    0.90,
    0.05,
    0.20
  )

  exposure <- data.frame(
    SNP = snps,
    beta.exposure = seq(0.1, 1.2, by = 0.1),
    se.exposure = 0.05,
    effect_allele.exposure = ea,
    other_allele.exposure = oa,
    pval.exposure = 1e-8,
    eaf.exposure = eaf,
    exposure = "exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  noise <- noise_sd * rep(c(0.5, -1, 1.5, -0.3, 0.8, -1.2), 2)
  eaf_out <- pmin(pmax(eaf + noise, 0.01), 0.99)

  outcome_ok <- data.frame(
    SNP = snps,
    beta.outcome = seq(0.05, 0.60, by = 0.05),
    se.outcome = 0.02,
    effect_allele.outcome = ea,
    other_allele.outcome = oa,
    pval.outcome = 0.01,
    eaf.outcome = eaf_out,
    outcome = "out",
    id.outcome = "out1",
    mr_keep.outcome = TRUE,
    pval_origin.outcome = "reported",
    stringsAsFactors = FALSE
  )

  outcome_bug <- outcome_ok
  outcome_bug$effect_allele.outcome <- oa
  outcome_bug$other_allele.outcome <- ea

  list(exposure = exposure, ok = outcome_ok, bug = outcome_bug)
}

# 40 non-palindromic SNPs: TwoSampleMR-format exposure plus outcomes in the
# format_gwas() outcome schema (rsids, effect_allele, ...), for exercising
# check_allele_orientation_gwas() and run_mr() with only 3 instruments.
# `snps` overrides the rsIDs (the first three are the instruments).
make_allele_gwas_fixture <- function(snps = paste0("rs", 1:40)) {
  n <- length(snps)
  ea <- rep(c("A", "C", "G", "T"), length.out = n)
  oa <- rep(c("G", "T", "A", "C"), length.out = n)
  eaf <- seq(0.05, 0.95, length.out = n)

  exposure <- data.frame(
    SNP = snps,
    beta.exposure = rep(c(0.2, 0.3, 0.4, 0.5), length.out = n),
    se.exposure = 0.05,
    effect_allele.exposure = ea,
    other_allele.exposure = oa,
    pval.exposure = rep(c(1e-9, 1e-3), length.out = n),
    eaf.exposure = eaf,
    exposure = "exp",
    id.exposure = "exp1",
    mr_keep.exposure = TRUE,
    pval_origin.exposure = "reported",
    stringsAsFactors = FALSE
  )

  outcome_ok <- data.frame(
    rsids = snps,
    beta = rep(c(0.02, 0.03, 0.04, 0.05), length.out = n),
    se = 0.02,
    eaf = eaf,
    pval = 0.05,
    n = 10000,
    effect_allele = ea,
    other_allele = oa,
    phenotype = "out",
    stringsAsFactors = FALSE
  )

  outcome_bug <- outcome_ok
  outcome_bug$effect_allele <- oa
  outcome_bug$other_allele <- ea

  list(
    exposure = exposure,
    ok = outcome_ok,
    bug = outcome_bug,
    instruments = snps[1:3]
  )
}
