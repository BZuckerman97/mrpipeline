# Create a coloc_result object

Create a coloc_result object

## Usage

``` r
new_coloc_result(
  coloc_abf = NULL,
  coloc_susie = NULL,
  coloc_signals = NULL,
  coloc_prop_test = NULL,
  n_snps = 0L,
  harmonised_data = data.frame(),
  methods_skipped = character(),
  params = list(),
  status = "success",
  status_reason = NULL,
  timing = numeric(0)
)
```

## Arguments

- coloc_abf:

  Output of
  [`coloc::coloc.abf()`](https://rdrr.io/pkg/coloc/man/coloc.abf.html),
  or `NULL`.

- coloc_susie:

  Output of
  [`coloc::coloc.susie()`](https://rdrr.io/pkg/coloc/man/coloc.susie.html),
  or `NULL`.

- coloc_signals:

  Output of
  [`coloc::coloc.signals()`](https://rdrr.io/pkg/coloc/man/coloc.signals.html),
  or `NULL`.

- coloc_prop_test:

  Output of `colocPropTest::coloc.prop.test()`, or `NULL`.

- n_snps:

  Integer. Number of SNPs used in the analysis.

- harmonised_data:

  Data frame of harmonised data used for analysis.

- methods_skipped:

  Named character vector: names are method names, values are reasons for
  skipping.

- params:

  List of all input parameters to
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

- status:

  Character. One of `"success"`, `"no_snps_in_region"`,
  `"too_few_snps"`, `"no_harmonised_variants"`. Default `"success"`.

- status_reason:

  Character or `NULL`. Human-readable explanation when
  `status != "success"`.

- timing:

  Named numeric vector of elapsed times (seconds) for each major step
  inside
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).
  Empty by default.

## Value

An object of class `coloc_result`.
