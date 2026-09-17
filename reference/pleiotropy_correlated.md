# Egger intercept from the correlated fit

`mr_egger(correl = TRUE)@Intercept` with its standard error and p-value,
in the columns of
[`TwoSampleMR::mr_pleiotropy_test()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_pleiotropy_test.html)
plus `ld_corrected`.

## Usage

``` r
pleiotropy_correlated(ld_input, harmonised)
```

## Arguments

- ld_input:

  An `MRInput` carrying the correlation matrix.

- harmonised:

  Harmonised instrument data frame (for the id columns).

## Value

A one-row data frame with columns `id.exposure`, `id.outcome`,
`outcome`, `exposure`, `egger_intercept`, `se`, `pval`, `ld_corrected`.
