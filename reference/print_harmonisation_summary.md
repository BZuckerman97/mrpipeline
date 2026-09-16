# Print a harmonisation breakdown

Shared by
[`summary.mr_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.mr_result.md)
and
[`summary.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.coloc_result.md).
Prints nothing when there is no harmonisation record to describe (older
result objects, or a run that never reached harmonisation).

## Usage

``` r
print_harmonisation_summary(harmonisation)
```

## Arguments

- harmonisation:

  The `harmonisation` field of an `mr_result` or `coloc_result`.

## Value

`invisible(NULL)`, called for its output.
