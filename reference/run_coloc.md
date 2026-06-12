# Perform colocalization analysis

Runs colocalization between a pre-formatted exposure and a raw outcome
dataset within a specified gene region. Supports coloc.abf, SuSiE-based
coloc (coloc.susie), coloc.signals, and colocPropTest. Returns a
`coloc_result` S3 object.

## Usage

``` r
run_coloc(
  exposure,
  outcome,
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
  verbose = TRUE
)
```

## Arguments

- exposure:

  Data frame of formatted exposure data (output of
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html)
  or `format_pqtl_*()` functions).

- outcome:

  Data frame of outcome summary statistics with standardised columns:
  `rsids`, `chr`, `pos`, `beta`, `se`, `eaf`, `pval`, `n`,
  `effect_allele`, `other_allele`.

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
  [`coloc::coloc.signals()`](https://rdrr.io/pkg/coloc/man/coloc.signals.html).
  Requires `"susie"` to have run successfully.

- `"prop_test"` – Proportionality test via
  `colocPropTest::coloc.prop.test()`. Requires `"signals"` to have run
  successfully and the `colocPropTest` package to be installed.

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
