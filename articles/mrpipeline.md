# Getting Started with mrpipeline

mrpipeline provides a streamlined interface for Mendelian randomisation
(MR) and colocalization analysis, with a focus on proteomic GWAS data
from deCODE and UKB-PPP. It wraps TwoSampleMR, coloc, and
MendelianRandomization into a consistent workflow with S3 result objects
and built-in sensitivity analyses.

## Installation

Install from GitHub:

``` r

# install.packages("pak")
pak::pak("BZuckerman97/mrpipeline")
```

## Quick start

``` r

library(mrpipeline)
```

### Mendelian randomisation

mrpipeline ships with bundled test datasets for CD40 protein and
Sjogren’s disease. Use these to explore the package without any external
data.

``` r

# Bundled datasets: cd40_exposure (formatted exposure), sjogren_outcome (outcome)
bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")

# Run cis-MR
mr_res <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,

  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = bfile,
  methods = c("ivw", "egger", "weighted_median")
)

# Inspect results
mr_res
summary(mr_res)

# Plot (requires ggplot2)
plot(mr_res, type = "scatter")
plot(mr_res, type = "forest")
```

The `mr_result` object stores the full results table, harmonised
instruments, F-statistics, Steiger filtering output, and any skipped
methods — accessible via `mr_res$results`, `mr_res$instruments`, etc.

### Colocalization

``` r

coloc_res <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20",
  gene_start = 44746911,
  gene_end = 44758502,
  bfile = bfile,
  methods = "abf"
)

coloc_res
summary(coloc_res)
plot(coloc_res, type = "pp_bar")
```

### Gene coordinate lookup

Look up genomic coordinates for HGNC gene symbols via Ensembl (requires
the `biomaRt` package):

``` r

coords <- get_gene_coords(c("CD40", "APOE"), build = "grch38")
coords
```

These coordinates can be passed directly to
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
via the `instrument_region` and `gene_*` arguments.

## Formatting exposure data

mrpipeline includes formatters for common proteomic GWAS sources:

- [`format_pqtl_decode()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_decode.md)
  — deCODE genetics
- [`format_pqtl_ukbppp()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_ukbppp.md)
  — UKB-PPP (Olink)
- [`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)
  — OneK1K single-cell eQTL

Each returns data formatted for TwoSampleMR, ready to pass to
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md).

## Further reading

- [`vignette("mrpipeline-user-guide")`](https://github.com/BZuckerman97/mrpipeline/articles/mrpipeline-user-guide.md)
  — detailed usage examples
- [`vignette("mrpipeline-developer-guide")`](https://github.com/BZuckerman97/mrpipeline/articles/mrpipeline-developer-guide.md)
  — architecture and internals
