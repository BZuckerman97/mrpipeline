# Check that compressed files can be read

Check that compressed files can be read

## Usage

``` r
check_gz_support(path, has_rutils = NULL)
```

## Arguments

- path:

  Character scalar file path. Non-character input is ignored.

- has_rutils:

  Logical, or `NULL` (default) to look `R.utils` up with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html). Only
  consulted when `path` is compressed. Exposed so tests can exercise the
  failure branch without uninstalling the package.

## Value

`invisible(TRUE)` if `path` is readable, otherwise aborts.
