# Plot an mr_result object

Creates diagnostic plots for Mendelian randomisation results using
TwoSampleMR plotting functions. Requires the `ggplot2` package.

## Usage

``` r
# S3 method for class 'mr_result'
plot(x, type = c("scatter", "forest", "funnel"), ...)
```

## Arguments

- x:

  An `mr_result` object.

- type:

  Character. Plot type: `"scatter"` (default), `"forest"`, or
  `"funnel"`.

- ...:

  Ignored.

## Value

A ggplot object (or list of ggplot objects for `"scatter"`). Returns
`NULL` invisibly if the result status is not `"success"` or if no
results are available.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
plot(result, type = "scatter")
plot(result, type = "forest")
plot(result, type = "funnel")
} # }
```
