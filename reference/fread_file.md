# Read a delimited file, checking compression support first

Thin wrapper around
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
that fails early, with an actionable message, when `path` is compressed
but the `R.utils` package (which `fread()` needs to decompress
`.gz`/`.bz2` files) is unavailable. `R.utils` is only a `Suggests` of
`data.table`, so declaring it in `mrpipeline`'s `Imports` is what
actually guarantees the gzipped inputs that GWAS summary statistics are
routinely distributed as (issue \#20). This guard is a backstop for a
broken library rather than an expected failure mode, and keeps the error
at the top of the call stack instead of deep inside `fread()`.

## Usage

``` r
fread_file(path, ...)
```

## Arguments

- path:

  Character scalar file path.

- ...:

  Passed to
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).

## Value

The value of
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).
