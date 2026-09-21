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
  methods = c("ivw_random", "egger", "weighted_median", "presso", "conmix", "steiger"),
  ld_correct = FALSE,
  exposure_n = NULL,
  presso_n_dist = 1000,
  plink_threads = plink_option("threads"),
  plink_memory = plink_option("memory"),
  allele_check = c("error", "warn", "none"),
  harmonise_action = 2,
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

  Character vector of methods to run. Named shortcuts (see
  [`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md)
  and the *Available methods* section): `"ivw_random"`, `"ivw_fixed"`,
  `"egger"`, `"weighted_median"`, `"presso"`, `"conmix"`, `"steiger"`,
  `"pleiotropy"`, `"heterogeneity"`, `"loo"`. You may also pass any raw
  method name from `TwoSampleMR::mr_method_list()$obj` that has no
  shortcut (e.g. `"mr_raps"`, `"mr_weighted_mode"`); raw names are never
  LD-corrected. The former shortcuts `"ivw"` and `"ivw_fe"` are accepted
  with a warning and mapped to `"ivw_random"` and `"ivw_fixed"`.

- ld_correct:

  Logical. Fit `ivw_random`, `ivw_fixed` and `egger` by GLS with the
  instruments' LD matrix (see the *LD correction* section). Requires
  `bfile`. Default `FALSE`.

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

- allele_check:

  Character. What to do when the allele orientation check finds that
  effect/other alleles look swapped between exposure and outcome – the
  signature of a GWAS file whose `A1`/`A2` mean REF/ALT, which silently
  inverts every beta (see
  [`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md),
  section *What does A1 mean?*). `"error"` (default) aborts, `"warn"`
  warns and continues, `"none"` runs the analysis regardless. The check
  harmonises the instruments together with up to 1000 further SNPs
  shared by the two datasets, so it works even for a cis-MR with a
  handful of instruments; it is skipped when fewer than 10 informative
  non-palindromic SNPs carry both allele frequencies. The full record is
  available afterwards from
  [`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)
  in every mode.

- harmonise_action:

  `1`, `2` (default) or `3`, passed to
  [`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html).
  `1` assumes every allele is on the forward strand; `2` infers the
  positive strand, resolving palindromic variants from their allele
  frequencies; `3` additionally drops every palindromic, ambiguous or
  incompatible SNP. Only palindrome handling differs – non-palindromic
  variants are aligned by allele letter at all three levels. Reach for
  `3` when the frequencies that level `2` relies on cannot be trusted: a
  failed allele orientation check makes every palindromic strand call in
  that pair unreliable, and a dataset without allele frequencies gives
  level `2` nothing to resolve them with.

- verbose:

  Logical. If `TRUE`, emit informational messages via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html).
  Warnings and errors are always emitted regardless. Default `TRUE`.

## Value

An `mr_result` object. Check `result$status` for `"success"` vs failure
reasons. Each `$results` row states the estimator that produced it:
`method` (the label), `model` (`"random"`, `"fixed"`, or `NA` where the
distinction does not apply) and `ld_corrected` (whether the LD matrix
was used for that fit), alongside `nsnp`, `b`, `se`, `pval` and the
`or`, `or_lci95`, `or_uci95` columns from
[`TwoSampleMR::generate_odds_ratios()`](https://mrcieu.github.io/TwoSampleMR/reference/generate_odds_ratios.html).
The `$timing` field contains a named numeric vector of elapsed seconds
for each major step.

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

Every method is looked up in the package's method registry – see
[`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md)
and the *Available methods* section below for what each shortcut runs,
whether it can be LD-corrected, and how many instruments it needs.
Dispatch depends on the number of instruments after clumping:

- 1 SNP: Wald ratio only; every multi-SNP method is skipped

- 2+ SNPs: `ivw_random`, `ivw_fixed`, `conmix`, `heterogeneity` and any
  raw TwoSampleMR method are attempted; `egger`, `weighted_median`,
  `presso`, `pleiotropy` and `loo` require \>= 3 SNPs

Raw TwoSampleMR methods (`mr_*` names from
`TwoSampleMR::mr_method_list()$obj` that have no shortcut, e.g.
`"mr_raps"`) are dispatched via
[`TwoSampleMR::mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.html)
under TwoSampleMR's own label; errors are caught and reported in
`$methods_skipped`. Names that are the engine behind a shortcut
(`mr_ivw`, `mr_ivw_fe`, `mr_egger_regression`, `mr_weighted_median`) are
refused with a pointer to the shortcut, so the same estimator cannot
enter unlabelled and uncorrected.

ConMix reports `se = NA`:
[`MendelianRandomization::mr_conmix()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_conmix.html)
returns a confidence interval that may be asymmetric or multi-modal
rather than a standard error, so the `or_lci95`/`or_uci95` columns are
`NA` for it too.

When `"egger"` is in `methods` and there are \>= 3 instruments,
[`TwoSampleMR::mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.html)
(the Egger intercept test) is always run automatically and its result
stored in `$pleiotropy`. You do not need to add `"pleiotropy"` to
`methods` separately. The `"pleiotropy"` shortcut remains available for
running the intercept test without Egger.

## LD correction

`ld_correct = TRUE` computes a signed LD matrix for the instruments from
`bfile` and re-orients it to the exposure's effect alleles. Instruments
absent from the reference panel, or with ambiguous palindromic alleles,
are dropped at that step, so an LD-corrected run can have fewer
instruments than the same call uncorrected.

Methods with a correlated form – `ivw_random`, `ivw_fixed` and `egger` –
are then fitted by generalised least squares through
`MendelianRandomization` with `correl = TRUE`: the weight matrix is
`diag(se_y) %*% R %*% diag(se_y)` in place of `diag(se_y^2)`, so two
instruments in LD are no longer counted as two independent looks at the
causal effect. The IVW estimator is pinned explicitly
(`model = "random"` for `ivw_random`, `"fixed"` for `ivw_fixed`) so that
each shortcut means the same thing at every instrument count;
`MendelianRandomization`'s own default would switch to fixed effects
below 4 instruments.

Every other `$results` method has no correlated form and runs on the
uncorrected data: a warning names each such method, and its row carries
`ld_corrected = FALSE`. The diagnostics come from the correlated fits
too: `$heterogeneity` holds the generalised Cochran Q for correlated
instruments (`Q = r' O^-1 r` on the GLS residuals, from
`mr_ivw()@Heter.Stat` and `mr_egger()@Heter.Stat`), `$pleiotropy` the
correlated Egger intercept (`mr_egger()@Intercept`), and `$loo` a
per-SNP correlated random-effects refit (a block-inverse update, so it
stays O(n^3)). Each of those frames carries an `ld_corrected` column on
both arms. `steiger` is the one thing left as-is: the Steiger direction
test compares per-SNP r^2 values and involves no weight matrix.
`ld_correct` is never silently ignored – it is either applied, or
visibly not applied. To compare corrected and uncorrected estimates,
call `run_mr()` twice and pass both results to
[`forest_plot()`](https://github.com/BZuckerman97/mrpipeline/reference/forest_plot.md)
as a named list: one `mr_result` is always one instrument set under one
weight matrix.

If the GLS weight matrix is near-singular (reciprocal condition number
below `1e-10`) `run_mr()` warns that every LD-corrected estimate is
unstable. The usual causes are identical or near-identical instruments
(r^2 ~ 1, from absent clumping or a manual set with a duplicated
variant) and more instruments than reference-panel individuals, which
makes the sample correlation matrix singular; clump more stringently or
drop the duplicate.

Random effects are multiplicative: the standard error is inflated by
`max(RSE, 1)`, never deflated, so when the instruments are
under-dispersed the random- and fixed-effect results coincide exactly.

## Available methods

Rendered from
[`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md).
`label` is the `method` value written to `$results`; `model` the effects
model; `LD` whether `ld_correct = TRUE` applies; `min n` the minimum
number of instruments.

|  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| shortcut | description | label | output | model | LD | min n |
| `ivw_random` | IVW, multiplicative random effects | IVW (random effects) | `$results` | random | yes | 2 |
| `ivw_fixed` | IVW, fixed effects | IVW (fixed effects) | `$results` | fixed | yes | 2 |
| `egger` | MR Egger regression | MR Egger | `$results` | random | yes | 3 |
| `weighted_median` | Weighted median | Weighted median | `$results` | – | no | 3 |
| `presso` | MR-PRESSO outlier test | MR-PRESSO | `$results` | – | no | 3 |
| `conmix` | Contamination mixture | ConMix | `$results` | – | no | 2 |
| `steiger` | Steiger directionality test | – | `$steiger` | – | no | 1 |
| `pleiotropy` | Egger intercept (pleiotropy) test | – | `$pleiotropy` | – | yes | 3 |
| `heterogeneity` | Cochran's Q heterogeneity test | – | `$heterogeneity` | – | yes | 2 |
| `loo` | Leave-one-out IVW | – | `$loo` | – | yes | 3 |

## See also

[`mr_methods()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_methods.md)
for the table of methods and what each supports.

## Examples

``` r
if (FALSE) { # \dontrun{
# Cis-MR using bundled CD40/Sjogren's data
bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
result <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  rsq_thresh = 0.3,
  bfile = bfile,
  methods = c("ivw_random", "egger", "weighted_median")
)
result
summary(result)

# LD-corrected: both IVW estimators by GLS; weighted median warns and
# runs uncorrected. Compare arms by passing both results to forest_plot().
corrected <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  rsq_thresh = 0.3,
  bfile = bfile,
  ld_correct = TRUE,
  methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
)
corrected$results[, c("method", "model", "ld_corrected", "b", "se")]
forest_plot(list("Uncorrected" = result, "LD-corrected" = corrected))
} # }
```
