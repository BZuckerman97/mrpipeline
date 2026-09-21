# Instrument pairs that are (near-)perfectly correlated

The usual reason a GLS weight matrix is singular: two instruments in
complete LD carry the same column, so `R` has a zero eigenvalue and
[`solve()`](https://rdrr.io/r/base/solve.html) fails. Used by
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
to name the offending pairs rather than report the condition number
alone (issue \#36).

## Usage

``` r
collinear_pairs(ld_matrix, threshold = 0.999)
```

## Arguments

- ld_matrix:

  Signed, aligned LD correlation matrix.

- threshold:

  Absolute correlation at or above which a pair is reported. Default
  `0.999`.

## Value

A data frame with columns `snp_a`, `snp_b` and `r`, ordered by
decreasing `abs(r)`; zero rows when no pair reaches `threshold`.
