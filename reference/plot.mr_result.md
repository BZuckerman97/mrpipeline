# Plot an mr_result object

Creates diagnostic plots for Mendelian randomisation results using
TwoSampleMR plotting functions. Requires the `ggplot2` package.

## Usage

``` r
# S3 method for class 'mr_result'
plot(x, type = c("scatter", "forest", "funnel", "loo"), ...)
```

## Arguments

- x:

  An `mr_result` object.

- type:

  Character. Plot type: `"scatter"` (default), `"forest"`, `"funnel"`,
  or `"loo"` (leave-one-out, via
  [`TwoSampleMR::mr_leaveoneout_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout_plot.html)
  on `x$loo`; requires `"loo"` to have been in `methods` when
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  was called).

- ...:

  Ignored.

## Value

A ggplot object (or list of ggplot objects for `"scatter"` and `"loo"`).
Returns `NULL` invisibly if the result status is not `"success"`, if no
results are available, or (for `type = "loo"`) if `x$loo` was not
computed.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline"),
  methods = c("ivw", "egger", "weighted_median", "loo")
)
plot(result, type = "scatter")
plot(result, type = "forest")
plot(result, type = "funnel")
plot(result, type = "loo")
} # }
```
