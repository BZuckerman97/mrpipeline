# CD40 exposure data (TwoSampleMR format)

50 SNPs from the CD40 cis region on chr20 (44646911-44858502),
pre-formatted as a TwoSampleMR exposure data frame. Derived from
`cd40_sumstats`. All SNPs overlap with the bundled LD reference panel
(`inst/extdata/ld_ref.*`) and with `sjogren_outcome`, so the two can be
used together in
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)/[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
examples without any external data.

## Usage

``` r
cd40_exposure
```

## Format

A data frame with 50 rows and the standard TwoSampleMR exposure columns
(`SNP`, `chr.exposure`, `pos.exposure`, `effect_allele.exposure`,
`other_allele.exposure`, `beta.exposure`, `se.exposure`,
`pval.exposure`, `eaf.exposure`, `samplesize.exposure`, `exposure`,
`id.exposure`, `mr_keep.exposure`, `pval_origin.exposure`).
