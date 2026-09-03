# Perform colocalization analysis

Runs colocalization between a pre-formatted exposure and a raw outcome
dataset within a specified gene region. Supports coloc.abf, SuSiE-based
coloc (coloc.susie), coloc.signals, and colocPropTest. Returns a
`coloc_result` S3 object.

## Usage

``` r
run_coloc(
  exposure,
  exposure_id,
  outcome,
  outcome_id,
  gene_chr,
  gene_start,
  gene_end,
  coloc_window = 10000L,
  exposure_n = NULL,
  outcome_n = NULL,
  exposure_type = c("quant", "cc"),
  outcome_type = c("quant", "cc"),
  exposure_s = NULL,
  outcome_s = NULL,
  exposure_sdY = 1,
  outcome_sdY = 1,
  bfile,
  plink_bin = NULL,
  methods = c("abf", "susie", "signals"),
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 1e-05,
  plink_threads = plink_option("threads"),
  plink_memory = plink_option("memory"),
  susie_maxit = 10000L,
  susie_repeat_until_convergence = FALSE,
  exclude_regions = NULL,
  ref_frq = NULL,
  verbose = TRUE
)
```

## Arguments

- exposure:

  Data frame of formatted exposure data (output of
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html)
  or `format_pqtl_*()` functions).

- exposure_id:

  Character. Label for the exposure (e.g. protein name or trait). Used
  in printed output and stored in the returned object.

- outcome:

  Data frame of outcome summary statistics with standardised columns:
  `rsids`, `chr`, `pos`, `beta`, `se`, `eaf`, `pval`, `n`,
  `effect_allele`, `other_allele`.

- outcome_id:

  Character. Label for the outcome (e.g. disease name or trait). Used in
  printed output and stored in the returned object.

- gene_chr:

  Character or integer. Chromosome of the gene/region.

- gene_start:

  Integer. Start position (bp) of the gene/region.

- gene_end:

  Integer. End position (bp) of the gene/region.

- coloc_window:

  Integer. Window in base pairs to extend around the gene region.
  Default `10000L` (10 kb).

- exposure_n:

  Integer. Exposure sample size. If `NULL`, inferred from
  `samplesize.exposure` column.

- outcome_n:

  Integer. Outcome sample size. If `NULL`, inferred from `n` column of
  `outcome`.

- exposure_type:

  Character. Type of exposure trait: `"quant"` or `"cc"`. Default
  `"quant"`.

- outcome_type:

  Character. Type of outcome trait: `"quant"` or `"cc"`. Default
  `"quant"`.

- exposure_s:

  Numeric. Proportion of cases for case-control exposure. Required when
  `exposure_type = "cc"`.

- outcome_s:

  Numeric. Proportion of cases for case-control outcome. Required when
  `outcome_type = "cc"`.

- exposure_sdY:

  Numeric. Standard deviation of the exposure trait (for quantitative
  traits). Default `1`.

- outcome_sdY:

  Numeric. Standard deviation of the outcome trait (for quantitative
  traits). Default `1`.

- bfile:

  Character. Path to PLINK bfile prefix (without .bed/.bim/.fam) for LD
  reference. Required.

- plink_bin:

  Character. Path to PLINK binary. Auto-detected if `NULL`.

- methods:

  Character vector of colocalization methods to run. Options: `"abf"`,
  `"susie"`, `"signals"`, `"prop_test"`. Default
  `c("abf", "susie", "signals")`.

- p1:

  Numeric. Prior probability a SNP is associated with trait 1. Default
  `1e-4`.

- p2:

  Numeric. Prior probability a SNP is associated with trait 2. Default
  `1e-4`.

- p12:

  Numeric. Prior probability a SNP is associated with both traits.
  Default `1e-5`.

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

- susie_maxit:

  Integer. Maximum iterations for
  [`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html).
  Default `10000L`. Increase if SuSiE warns about non-convergence;
  decrease for faster exploratory runs.

- susie_repeat_until_convergence:

  Logical. Passed to
  [`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html).
  Default `FALSE` – prevents infinite loops when SuSiE has not converged
  within `susie_maxit` iterations.

- exclude_regions:

  Data frame with columns `chr`, `start`, `end` defining genomic regions
  to exclude SNPs from before colocalization, or `NULL`. SNPs in both
  the exposure and outcome that fall within any listed region are
  dropped before harmonisation. For example, to exclude the MHC:
  `data.frame(chr = "6", start = 26e6, end = 34e6)`.

- ref_frq:

  Character. Path to a PLINK `.frq` file (produced by `plink --freq`)
  used to patch EAF for SNPs missing it in both the exposure and
  outcome. Columns expected: `SNP`, `A1`, `A2`, `MAF`. The A1 allele is
  matched against `effect_allele.outcome` (post-harmonisation) to orient
  the frequency correctly. `NULL` (default) disables the lookup.

- verbose:

  Logical. If `TRUE`, emit informational messages via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html).
  Warnings and errors are always emitted regardless. Default `TRUE`.

## Value

A `coloc_result` object. Check `result$status` for `"success"` vs
failure reasons. The `$timing` field contains a named numeric vector of
elapsed seconds for each major step.

## Methods

- `"abf"` – Approximate Bayes Factor colocalization via
  [`coloc::coloc.abf()`](https://rdrr.io/pkg/coloc/man/coloc.abf.html)

- `"susie"` – SuSiE fine-mapping + colocalization via
  [`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html) and
  [`coloc::coloc.susie()`](https://rdrr.io/pkg/coloc/man/coloc.susie.html)

- `"signals"` – Multi-signal colocalization via
  [`coloc::coloc.signals()`](https://rdrr.io/pkg/coloc/man/coloc.signals.html),
  called with `method = "cond"` (LD-based conditioning) on the same
  `dataset_exp`/`dataset_out` objects used by `"abf"`.
  [`coloc::coloc.signals()`](https://rdrr.io/pkg/coloc/man/coloc.signals.html)
  performs its own internal signal detection (`finemap.signals()`) and
  does not consume `"susie"`'s `runsusie()` output – despite the name,
  it's an independent multi-signal method, not a wrapper around SuSiE's
  credible sets. It is still gated on `"susie"` having been requested
  and having run successfully (skipped otherwise), for consistency with
  this function's `methods` design rather than because it needs SuSiE's
  output.

- `"prop_test"` – Proportionality test via
  `colocPropTest::coloc.prop.test()`. Requires `"signals"` to have run
  successfully and the `colocPropTest` package to be installed.

## MAF / EAF handling

[`coloc::coloc.abf()`](https://rdrr.io/pkg/coloc/man/coloc.abf.html)
takes a single MAF per SNP, shared between the exposure and outcome
datasets, and errors out for the *entire* dataset if even one SNP has a
missing MAF. Some GWASes (e.g. the CRP exposure used in this pipeline)
report no EAF column at all. To avoid losing every SNP over a single
dataset's missing EAF, `run_coloc()` builds the MAF from `eaf.exposure`,
falling back to `eaf.outcome` wherever the exposure value is missing
(and vice versa). A
[`cli::cli_warn()`](https://cli.r-lib.org/reference/cli_abort.html) is
emitted whenever either side has missing EAF, so the fallback is visible
in the warnings captured by `targets` metadata
(`targets::tar_meta(fields = warnings)`) even though it doesn't print
live during a `tar_make()` run. Any SNP missing EAF on *both* sides is
dropped (with a warning) before colocalization, since there is no value
left to fall back to.

## LD matrix orientation

[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md)
(via
[`ieugwasr::ld_matrix()`](https://mrcieu.github.io/ieugwasr/reference/ld_matrix.html),
`--r square --keep-allele-order`) returns a **signed** correlation
matrix whose sign is anchored to the reference panel's own, arbitrary A1
allele for each SNP – not to
`effect_allele.exposure`/`effect_allele.outcome`. Roughly half of SNPs,
at random, will have the panel's A1 differ from the harmonised effect
allele.
[`coloc::coloc.abf()`](https://rdrr.io/pkg/coloc/man/coloc.abf.html) is
unaffected by this: it never uses `LD`, and each SNP's Bayes factor
depends only on `beta^2`, so a per-SNP sign flip changes nothing.
[`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html) (and
the
[`susieR::susie_rss()`](https://rdrr.io/pkg/susieR/man/susie_rss.html)
it wraps), however, fits a joint model across all SNPs from `LD` and
`beta`/`varbeta` together – if even a subset of SNPs have a sign
inconsistent with `LD`'s convention, the fit becomes internally
incoherent, which is exactly what `susie_rss()`'s `check_prior` safety
check (comparing the fitted prior variance against
`100 * max(abs(z))^2`) is designed to catch, surfacing as the error
`"the estimated prior variance is unreasonably large"`. If you ever see
that error and `"abf"` succeeded on the same data, an LD-orientation bug
– not LD-reference-panel quality – is the first thing to suspect.
[`align_to_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/align_to_ld_matrix.md)
corrects for this automatically: it classifies each SNP as matching,
needing a sign flip, or unresolvable against the reference panel (also
dropping palindromic SNPs with an EAF in the ambiguous 0.42-0.58 zone,
which can't be resolved by allele identity alone), and returns a
`ld_sign` vector multiplied into `dataset_exp$beta`/`dataset_out$beta`
immediately before they're passed to
[`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html) – see
that function's documentation for the full algorithm.

## See also

[`new_coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/new_coloc_result.md)
for the S3 class structure,
[`print.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/print.coloc_result.md)
and
[`summary.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.coloc_result.md)
for display methods.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_coloc(
  exposure = formatted_exposure,
  outcome = outcome_data,
  gene_chr = 20, gene_start = 44746911, gene_end = 44758502,
  bfile = "/path/to/ld_reference",
  methods = "abf"
)
print(result)
summary(result)
} # }
```
