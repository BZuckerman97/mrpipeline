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
  plink_threads = plink_option("threads"),
  plink_memory = plink_option("memory"),
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

  Character vector of MR methods to run. Named shortcuts: `"ivw"` (IVW
  random effects), `"ivw_fe"` (IVW fixed effects), `"egger"` (MR Egger),
  `"weighted_median"` (weighted median), `"presso"` (MR-PRESSO),
  `"conmix"` (ContMix), `"steiger"` (Steiger filtering), `"pleiotropy"`
  (Egger intercept test; result stored in `$pleiotropy`, not
  `$results`). You may also pass any method name from
  `TwoSampleMR::mr_method_list()$obj` directly (e.g.
  `"mr_simple_median"`, `"mr_raps"`). Note: `"ivw_fe"` does not support
  `ld_correct = TRUE`.

- ld_correct:

  Logical. Use LD-corrected IVW/Egger via the `MendelianRandomization`
  package. Requires `bfile`. Default `FALSE`.

- exposure_n:

  Numeric. Exposure sample size. If `NULL`, inferred from
  `samplesize.exposure` column.

- presso_n_dist:

  Integer. Number of distributions for MR-PRESSO. Default `1000`.

- plink_threads:

  Integer. Number of threads for PLINK. `NULL` (default) lets PLINK
  auto-detect. Read from `getOption("mrpipeline.plink_threads")` or the
  `MRPIPELINE_PLINK_THREADS` environment variable via
  [`plink_option()`](https://github.com/BZuckerman97/mrpipeline/reference/plink_option.md).

- plink_memory:

  Integer. Memory limit in MB for PLINK. `NULL` (default) lets PLINK
  auto-detect. Read from `getOption("mrpipeline.plink_memory")` or the
  `MRPIPELINE_PLINK_MEMORY` environment variable via
  [`plink_option()`](https://github.com/BZuckerman97/mrpipeline/reference/plink_option.md).

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

- 1 SNP: Wald ratio only; all other methods skipped

- 2+ SNPs: IVW, IVW-FE, ConMix, Steiger, and any raw TwoSampleMR methods
  are attempted; Egger, weighted median, and PRESSO require \>= 3 SNPs

- 3+ SNPs: all methods in `methods` are attempted

Generic TwoSampleMR methods (raw `mr_*` names) are dispatched via
[`TwoSampleMR::mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.html)
and errors are caught and reported as skipped.

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
