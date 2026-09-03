# Sjogren's disease outcome data (standardised format)

50 SNPs in the CD40 cis region on chr20, in standardised outcome format,
matching
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)'s
`type = "outcome"` schema. One real SNP (`rs9074`) from
`sjogren_sumstats`; 49 synthetic SNPs with realistic effect sizes to
match the bundled LD reference panel (`inst/extdata/ld_ref.*`) and
`cd40_exposure`.

## Usage

``` r
sjogren_outcome
```

## Format

A data frame with 50 rows and columns `rsids`, `chr`, `pos`, `beta`,
`se`, `eaf`, `pval`, `n`, `effect_allele`, `other_allele`, `phenotype`.
