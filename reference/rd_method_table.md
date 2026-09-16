# Render the concise method table as roxygen lines

Used through `@eval` in
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)'s
documentation, so
[`?run_mr`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
is regenerated from the registry on every `devtools::document()` and
cannot drift from what
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
accepts.

## Usage

``` r
rd_method_table()
```

## Value

Character vector of roxygen lines (a markdown table inside a section).
