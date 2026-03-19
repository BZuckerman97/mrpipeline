# Create an mr_result object

Create an mr_result object

## Usage

``` r
new_mr_result(
  results = data.frame(),
  instruments = data.frame(),
  f_stats = list(per_snp = numeric(), mean = NA_real_, min = NA_real_),
  steiger = NULL,
  methods_skipped = character(),
  ld_matrix = NULL,
  params = list(),
  status = "success",
  status_reason = NULL,
  timing = numeric(0)
)
```

## Arguments

- results:

  Data frame with columns: exposure, outcome, method, nsnp, b, se, pval.

- instruments:

  Data frame of harmonised (and clumped) instrument data.

- f_stats:

  List with elements `per_snp` (numeric vector), `mean` (numeric
  scalar), `min` (numeric scalar).

- steiger:

  Output of
  [`TwoSampleMR::steiger_filtering()`](https://mrcieu.github.io/TwoSampleMR/reference/steiger_filtering.html),
  or `NULL`.

- methods_skipped:

  Named character vector: names are method names, values are reasons for
  skipping.

- ld_matrix:

  LD correlation matrix if `ld_correct = TRUE`, or `NULL`.

- params:

  List of all input parameters to
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md).

- status:

  Character. One of `"success"`, `"no_instruments"`,
  `"no_harmonised_variants"`. Default `"success"`.

- status_reason:

  Character or `NULL`. Human-readable explanation when
  `status != "success"`.

- timing:

  Named numeric vector of elapsed times (seconds) for each major step
  inside
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md).
  Empty by default.

## Value

An object of class `mr_result`.
