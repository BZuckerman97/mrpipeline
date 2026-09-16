# The MR method registry

The single source of truth for every method
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
can run: its shortcut name, the label written to `$results$method`,
whether it is a fixed- or random-effects estimator, whether
`ld_correct = TRUE` applies to it, the minimum number of instruments it
needs, where its output lands on the `mr_result`, and the function that
actually runs on each LD path.
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
reads validation, dispatch, skip reasons and the LD warnings from this
table, and the method tables in
[`?run_mr`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and the vignettes are rendered from it, so the documentation cannot
drift from the code. Adding a method means adding a row here.

## Usage

``` r
mr_method_registry()
```

## Value

A data frame with one row per method path and columns `shortcut`,
`description`, `label`, `output`, `model`, `ld_correctable`,
`min_instruments`, `engine`, `engine_ld`.

## Details

Two rows have no `shortcut`: the Wald ratio, which
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
uses automatically when exactly one instrument survives, and the raw
TwoSampleMR passthrough, one row standing for every
`TwoSampleMR::mr_method_list()$obj` name that has no shortcut.

`ld_correctable` is declared, not inferred: nothing inspects upstream
function signatures at run time. The row is kept honest by `engine_ld` –
a method claims LD support only by naming the function that provides it,
and the tests assert `ld_correctable == !is.na(engine_ld)`. There is
deliberately no free-text notes column: notes about upstream internals
go stale silently when a dependency changes, whereas
`engine`/`engine_ld` name what mrpipeline calls and so stay true as long
as the dispatch does.
