# Methods available to run_mr()

Returns the table of methods
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
understands: what each shortcut runs, whether it is a fixed- or
random-effects estimator, whether `ld_correct = TRUE` applies to it, how
many instruments it needs, and where its result lands on the returned
`mr_result`. The table is the same object
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
dispatches from (see
[`mr_method_registry()`](https://github.com/BZuckerman97/mrpipeline/reference/mr_method_registry.md)),
so it is always current.

## Usage

``` r
mr_methods(detail = c("concise", "full"))
```

## Arguments

- detail:

  `"concise"` (default) lists the ten shortcuts that can be passed to
  `methods=`. `"full"` adds the two paths that have no single name – the
  automatic Wald ratio used at exactly one instrument, and the raw
  TwoSampleMR passthrough – plus the `engine` and `engine_ld` columns
  naming the function that runs on each LD path.

## Value

A data frame with columns:

- `shortcut`:

  Name accepted by `run_mr(methods = )`. `NA` for the two rows that
  cannot be requested by a single name.

- `description`:

  What the method does.

- `label`:

  The `method` value written to `$results`; `NA` for diagnostics, which
  produce no results row.

- `output`:

  Where the result lands: `$results`, `$steiger`, `$pleiotropy`,
  `$heterogeneity` or `$loo`.

- `model`:

  `"random"`, `"fixed"`, or `NA` where the distinction does not apply.
  Random effects are multiplicative: the standard error is inflated by
  `max(RSE, 1)` and never deflated.

- `ld_correctable`:

  Whether `ld_correct = TRUE` changes this method. Methods with `FALSE`
  run on uncorrected data and warn.

- `min_instruments`:

  Minimum number of instruments; below it the method is skipped and
  `$methods_skipped` records why.

- `engine`, `engine_ld`:

  (`detail = "full"` only) The function mrpipeline calls when
  `ld_correct` is `FALSE` / `TRUE`. Documentary – users never pass an
  engine;
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  picks it from the shortcut and `ld_correct`.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)

## Examples

``` r
mr_methods()
#>           shortcut                        description                label
#> 1       ivw_random IVW, multiplicative random effects IVW (random effects)
#> 2        ivw_fixed                 IVW, fixed effects  IVW (fixed effects)
#> 3            egger                MR Egger regression             MR Egger
#> 4  weighted_median                    Weighted median      Weighted median
#> 5           presso             MR-PRESSO outlier test            MR-PRESSO
#> 6           conmix              Contamination mixture               ConMix
#> 7          steiger        Steiger directionality test                 <NA>
#> 8       pleiotropy  Egger intercept (pleiotropy) test                 <NA>
#> 9    heterogeneity     Cochran's Q heterogeneity test                 <NA>
#> 10             loo                  Leave-one-out IVW                 <NA>
#>            output  model ld_correctable min_instruments
#> 1        $results random           TRUE               2
#> 2        $results  fixed           TRUE               2
#> 3        $results random           TRUE               3
#> 4        $results   <NA>          FALSE               3
#> 5        $results   <NA>          FALSE               3
#> 6        $results   <NA>          FALSE               2
#> 7        $steiger   <NA>          FALSE               1
#> 8     $pleiotropy   <NA>          FALSE               3
#> 9  $heterogeneity   <NA>          FALSE               2
#> 10           $loo   <NA>          FALSE               3
mr_methods(detail = "full")
#>           shortcut                        description                label
#> 1       ivw_random IVW, multiplicative random effects IVW (random effects)
#> 2        ivw_fixed                 IVW, fixed effects  IVW (fixed effects)
#> 3            egger                MR Egger regression             MR Egger
#> 4  weighted_median                    Weighted median      Weighted median
#> 5           presso             MR-PRESSO outlier test            MR-PRESSO
#> 6           conmix              Contamination mixture               ConMix
#> 7             <NA>      Wald ratio, single instrument           Wald ratio
#> 8             <NA>       Any other TwoSampleMR method                 <NA>
#> 9          steiger        Steiger directionality test                 <NA>
#> 10      pleiotropy  Egger intercept (pleiotropy) test                 <NA>
#> 11   heterogeneity     Cochran's Q heterogeneity test                 <NA>
#> 12             loo                  Leave-one-out IVW                 <NA>
#>            output  model ld_correctable min_instruments
#> 1        $results random           TRUE               2
#> 2        $results  fixed           TRUE               2
#> 3        $results random           TRUE               3
#> 4        $results   <NA>          FALSE               3
#> 5        $results   <NA>          FALSE               3
#> 6        $results   <NA>          FALSE               2
#> 7        $results   <NA>          FALSE               1
#> 8        $results   <NA>          FALSE               2
#> 9        $steiger   <NA>          FALSE               1
#> 10    $pleiotropy   <NA>          FALSE               3
#> 11 $heterogeneity   <NA>          FALSE               2
#> 12           $loo   <NA>          FALSE               3
#>                               engine
#> 1                TwoSampleMR::mr_ivw
#> 2             TwoSampleMR::mr_ivw_fe
#> 3   TwoSampleMR::mr_egger_regression
#> 4    TwoSampleMR::mr_weighted_median
#> 5         TwoSampleMR::run_mr_presso
#> 6  MendelianRandomization::mr_conmix
#> 7         TwoSampleMR::mr_wald_ratio
#> 8                    TwoSampleMR::mr
#> 9     TwoSampleMR::steiger_filtering
#> 10   TwoSampleMR::mr_pleiotropy_test
#> 11     TwoSampleMR::mr_heterogeneity
#> 12       TwoSampleMR::mr_leaveoneout
#>                                           engine_ld
#> 1  MendelianRandomization::mr_ivw(model = "random")
#> 2   MendelianRandomization::mr_ivw(model = "fixed")
#> 3                  MendelianRandomization::mr_egger
#> 4                                              <NA>
#> 5                                              <NA>
#> 6                                              <NA>
#> 7                                              <NA>
#> 8                                              <NA>
#> 9                                              <NA>
#> 10                                             <NA>
#> 11                                             <NA>
#> 12                                             <NA>
```
