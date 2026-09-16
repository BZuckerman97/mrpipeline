# Empty `$results` frame carrying the full column schema

Used as the `results` default of
[`new_mr_result()`](https://github.com/BZuckerman97/mrpipeline/reference/new_mr_result.md)
(so every early return has the same shape as a successful run) and by
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
when no method produced an estimate. The OR columns are appended by
[`TwoSampleMR::generate_odds_ratios()`](https://mrcieu.github.io/TwoSampleMR/reference/generate_odds_ratios.html)
in
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
not here.

## Usage

``` r
empty_mr_results()
```

## Value

A zero-row data frame with columns `exposure`, `outcome`, `method`,
`nsnp`, `b`, `se`, `pval`, `ld_corrected`, `model`.
