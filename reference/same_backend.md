# Are two lazy tables on the same database connection?

`TRUE` only for two dbplyr tables sharing one connection, the case in
which
[`shared_rsids()`](https://github.com/BZuckerman97/mrpipeline/reference/shared_rsids.md)
can leave the intersection to the database.

## Usage

``` r
same_backend(x, y)
```

## Arguments

- x, y:

  Data frames or lazy tables.

## Value

Logical scalar.
