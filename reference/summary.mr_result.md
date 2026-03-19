# Summarise an mr_result object

Displays the full results table, F-statistics, Steiger summary, Egger
intercept (if available), and skipped methods.

## Usage

``` r
# S3 method for class 'mr_result'
summary(object, ...)
```

## Arguments

- object:

  An `mr_result` object.

- ...:

  Ignored.

## Value

Invisibly returns `object`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
summary(result)
} # }
```
