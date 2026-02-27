# Example: CD40 → Sjögren's Disease (MR + colocalization)
#
# Uses bundled test data and the minimal LD reference panel shipped with the
# package. All 50 SNPs overlap across exposure, outcome, and LD panel, so this
# runs without any external downloads.

library(mrpipeline)

# --- LD reference panel path ------------------------------------------------

bfile <- sub(
  "\\.bed$",
  "",
  system.file("extdata", "ld_ref.bed", package = "mrpipeline")
)

# --- MR analysis ------------------------------------------------------------

mr_res <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = bfile,
  methods = c("ivw", "egger", "weighted_median", "steiger")
)

print(mr_res)
summary(mr_res)

# --- Colocalization ---------------------------------------------------------

coloc_res <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20",
  gene_start = 44746911,
  gene_end = 44758502,
  exposure_n = 33655,
  outcome_n = 41420,
  bfile = bfile,
  methods = "abf"
)

print(coloc_res)
summary(coloc_res)

# Access ABF posterior probabilities directly
coloc_res$coloc_abf$summary
