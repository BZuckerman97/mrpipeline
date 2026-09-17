# Identifier columns every diagnostics frame starts with

The `id.exposure`, `id.outcome`, `outcome`, `exposure` values
TwoSampleMR puts at the front of `mr_heterogeneity()` /
`mr_pleiotropy_test()` output, so the LD-corrected frames built here
have the same shape as the uncorrected ones.

## Usage

``` r
diag_ids(harmonised)
```

## Arguments

- harmonised:

  Harmonised instrument data frame.

## Value

A one-row data frame.
