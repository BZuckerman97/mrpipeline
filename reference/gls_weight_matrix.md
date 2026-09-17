# GLS weight matrix for correlated instruments

`O = diag(se_y) R diag(se_y)`: the matrix
[`MendelianRandomization::mr_ivw()`](https://rdrr.io/pkg/MendelianRandomization/man/mr_ivw.html)
and `mr_egger()` solve when `correl = TRUE`, and the one whose inverse
[`loo_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/loo_correlated.md)
updates.

## Usage

``` r
gls_weight_matrix(ld_matrix, se_outcome)
```

## Arguments

- ld_matrix:

  Signed, aligned LD correlation matrix.

- se_outcome:

  Outcome standard errors, in the matrix's row order.

## Value

A numeric matrix.
