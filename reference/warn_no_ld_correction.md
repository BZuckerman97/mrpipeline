# Warn that a requested method has no LD-corrected form

Called by
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
for every `$results`-producing method that runs while
`ld_correct = TRUE` but has no correlated implementation (registry
`ld_correctable = FALSE`), so `ld_correct` is never silently ignored: it
is either applied, or visibly not applied. Diagnostics (`steiger`,
`pleiotropy`, `heterogeneity`, `loo`) do not warn – they produce no
estimate row.

## Usage

``` r
warn_no_ld_correction(method)
```

## Arguments

- method:

  The shortcut or raw method name, as the user passed it.

## Value

`NULL`, invisibly.
