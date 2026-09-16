# Create an mr_result object

Create an mr_result object

## Usage

``` r
new_mr_result(
  results = empty_mr_results(),
  instruments = data.frame(),
  harmonisation = data.frame(),
  f_stats = list(per_snp = numeric(), mean = NA_real_, min = NA_real_),
  steiger = NULL,
  pleiotropy = NULL,
  heterogeneity = NULL,
  loo = NULL,
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

  Data frame with columns: exposure, outcome, method, nsnp, b, se, pval,
  ld_corrected (logical: whether the LD matrix was used for that fit),
  model (`"random"`, `"fixed"` or `NA`), plus or, or_lci95, or_uci95
  (and lo_ci, up_ci) from
  [`TwoSampleMR::generate_odds_ratios()`](https://mrcieu.github.io/TwoSampleMR/reference/generate_odds_ratios.html)
  on a successful run. Defaults to
  [`empty_mr_results()`](https://github.com/BZuckerman97/mrpipeline/reference/empty_mr_results.md).

- instruments:

  Data frame of harmonised (and clumped) instrument data: the kept
  variants the MR estimates are computed from.

- harmonisation:

  Data frame. The complete, unfiltered
  [`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
  output – every candidate variant, with the `mr_keep`, `palindromic`,
  `ambiguous` and `remove` flags that explain why each one was or was
  not carried forward. `instruments` is the subset of this that
  survived; see
  [`harmonisation_summary()`](https://github.com/BZuckerman97/mrpipeline/reference/harmonisation_summary.md).

- f_stats:

  List with elements `per_snp` (numeric vector), `mean` (numeric
  scalar), `min` (numeric scalar).

- steiger:

  Output of
  [`TwoSampleMR::steiger_filtering()`](https://mrcieu.github.io/TwoSampleMR/reference/steiger_filtering.html),
  or `NULL`.

- pleiotropy:

  Output of
  [`TwoSampleMR::mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.html),
  or `NULL`.

- heterogeneity:

  Output of
  [`TwoSampleMR::mr_heterogeneity()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_heterogeneity.html)
  (Cochran's Q per method), or `NULL`.

- loo:

  Output of
  [`TwoSampleMR::mr_leaveoneout()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout.html),
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
