# Warn that a requested method has no LD-corrected form

Called by
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
for every `$results`-producing method that runs while
`ld_correct = TRUE` but has no correlated implementation (registry
`ld_correctable = FALSE`), so `ld_correct` is never silently ignored: it
is either applied, or visibly not applied. The diagnostics `pleiotropy`,
`heterogeneity` and `loo` never reach here: they are computed from the
correlated fits (see
[`heterogeneity_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/heterogeneity_correlated.md),
[`pleiotropy_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/pleiotropy_correlated.md),
[`loo_correlated()`](https://github.com/BZuckerman97/mrpipeline/reference/loo_correlated.md)).
`steiger` is the one method exempt on purpose – the Steiger direction
test compares per-SNP r^2 values and involves no weight matrix, so LD
correction is not a concept for it.

## Usage

``` r
warn_no_ld_correction(method)
```

## Arguments

- method:

  The shortcut or raw method name, as the user passed it.

## Value

`NULL`, invisibly.
