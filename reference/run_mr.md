# Perform Mendelian randomisation analysis

Runs MR with automatic instrument selection (cis-MR, genome-wide, or
manual) and optional sensitivity analyses. Returns an `mr_result` S3
object.

## Usage

``` r
run_mr(
  exposure,
  exposure_id,
  outcome,
  outcome_id,
  instrument_region = NULL,
  window = 100000L,
  pval_thresh = 5e-08,
  rsq_thresh = 0.001,
  bfile = NULL,
  plink_bin = NULL,
  pop = "EUR",
  instruments = NULL,
  instruments_strict = FALSE,
  exclude_regions = NULL,
  methods = c("ivw", "egger", "weighted_median", "presso", "conmix", "steiger"),
  ld_correct = FALSE,
  exposure_n = NULL,
  presso_n_dist = 1000,
  verbose = TRUE
)
```

## Arguments

- exposure:

  Data frame of formatted exposure data (output of
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html)
  or `format_pqtl_*()` functions).

- exposure_id:

  Character. Identifier for the exposure (e.g. protein name).

- outcome:

  Data frame of outcome summary statistics with standardised columns:
  `rsids`, `chr`, `pos`, `beta`, `se`, `eaf`, `pval`, `n`,
  `effect_allele`, `other_allele`. Formatted internally via
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html).

- outcome_id:

  Character. Identifier for the outcome (e.g. disease name).

- instrument_region:

  List with elements `chromosome`, `start`, `end` defining the cis
  region. `NULL` for genome-wide or manual mode.

- window:

  Integer. Window (in bp) to extend either side of `instrument_region`.
  Default `100000L`.

- pval_thresh:

  Numeric. P-value threshold for instrument selection. Default `5e-8`.

- rsq_thresh:

  Numeric. R-squared clumping threshold. Default `0.001`.

- bfile:

  Character. Path to PLINK bfile prefix for local LD operations.
  Required when `ld_correct = TRUE`.

- plink_bin:

  Character. Path to PLINK binary. Auto-detected if `NULL`.

- pop:

  Character. Population for API-based LD clumping. Default `"EUR"`.

- instruments:

  Character vector of rsIDs for manual instrument mode, or `NULL`.

- instruments_strict:

  Logical. If `TRUE`, error when manual instruments are missing from
  exposure data. If `FALSE`, warn. Default `FALSE`.

- exclude_regions:

  Data frame with columns `chr`, `start`, `end` defining genomic regions
  to exclude instruments from, or `NULL`. For example, to exclude the
  MHC region: `data.frame(chr = "6", start = 26e6, end = 34e6)`.

- methods:

  Character vector of MR methods to run. Options: `"ivw"`, `"egger"`,
  `"weighted_median"`, `"presso"`, `"conmix"`, `"steiger"`.

- ld_correct:

  Logical. Use LD-corrected IVW/Egger via the `MendelianRandomization`
  package. Requires `bfile`. Default `FALSE`.

- exposure_n:

  Numeric. Exposure sample size. If `NULL`, inferred from
  `samplesize.exposure` column.

- presso_n_dist:

  Integer. Number of distributions for MR-PRESSO. Default `1000`.

- verbose:

  Logical. If `TRUE`, emit informational messages via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html).
  Warnings and errors are always emitted regardless. Default `TRUE`.

## Value

An `mr_result` object. Check `result$status` for `"success"` vs failure
reasons. The `$timing` field contains a named numeric vector of elapsed
seconds for each major step.

## Instrument selection modes

Exactly one of three modes is used, determined by the combination of
`instruments` and `instrument_region`:

- **Cis-MR** (`instrument_region` provided, `instruments = NULL`):
  filters `exposure` to the cis region defined by `instrument_region`
  +/- `window`, applies `pval_thresh`, then LD-clumps.

- **Genome-wide** (`instrument_region = NULL`, `instruments = NULL`):
  filters `exposure` by `pval_thresh` only, then LD-clumps.

- **Manual** (`instruments` provided): uses the supplied rsIDs directly.
  `instruments_strict` controls whether missing IDs are an error or
  warning.

## Method dispatch

Methods are dispatched based on the number of instruments after
clumping:

- 1 SNP: Wald ratio only

- 2 SNPs: IVW (+ ConMix/Steiger if requested);
  Egger/weighted_median/PRESSO skipped

- 3+ SNPs: all methods in `methods` are attempted

When `ld_correct = TRUE`, IVW and Egger use
[`MendelianRandomization::mr_ivw()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_ivw.html)
and
[`MendelianRandomization::mr_egger()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_egger.html)
with `correl = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Cis-MR using bundled CD40/Sjogren's data
bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")
result <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = bfile,
  methods = c("ivw", "egger", "weighted_median")
)
result
summary(result)
} # }
```
