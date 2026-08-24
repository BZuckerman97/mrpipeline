
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mrpipeline

<!-- badges: start -->

[![R-CMD-check](https://github.com/BZuckerman97/mrpipeline/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BZuckerman97/mrpipeline/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/BZuckerman97/mrpipeline/graph/badge.svg)](https://app.codecov.io/gh/BZuckerman97/mrpipeline)
<!-- badges: end -->

mrpipeline provides a streamlined interface for Mendelian randomisation
(MR) and colocalization analysis, with a focus on proteomic GWAS data
(deCODE, UKB-PPP). It wraps TwoSampleMR, coloc, and
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

# mrpipeline ships with bundled test datasets for CD40 protein and Sjogren's
# disease -- cd40_exposure (formatted exposure) and sjogren_outcome (outcome)
# -- plus a minimal LD reference panel, so this runs without any external
# data.
bfile <- sub(
  "\\.bed$", "",
  system.file("extdata", "ld_ref.bed", package = "mrpipeline")
)

# Run cis-MR
mr_res <- run_mr(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = bfile
)
mr_res
summary(mr_res)

# Run colocalization
coloc_res <- run_coloc(
  exposure = cd40_exposure,
  exposure_id = "CD40",
  outcome = sjogren_outcome,
  outcome_id = "SjD",
  gene_chr = "20",
  gene_start = 44746911,
  gene_end = 44758502,
  bfile = bfile
)
coloc_res
summary(coloc_res)
```

## Vignettes

- `vignette("mrpipeline-user-guide")` — end-to-end usage examples
- `vignette("mrpipeline-developer-guide")` — architecture and internals
