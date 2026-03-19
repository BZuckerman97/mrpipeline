# Print an mr_result object

Displays a one-line summary: exposure -\> outcome, primary estimate,
number of SNPs, and mean F-statistic.

## Usage

``` r
# S3 method for class 'mr_result'
print(x, ...)
```

## Arguments

- x:

  An `mr_result` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(
  exposure = cd40_exposure, exposure_id = "CD40",
  outcome = sjogren_outcome, outcome_id = "SjD",
  instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
print(result)
} # }
```
