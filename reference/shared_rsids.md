# rsIDs shared by an exposure and an outcome, minus the instruments

The shared-SNP step of
[`check_allele_orientation_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation_gwas.md).
Each input is cut to its one rsID column and the two are combined with
[`dplyr::intersect()`](https://generics.r-lib.org/reference/setops.html),
which on two tables from the same database connection compiles to a SQL
`INTERSECT`, so only the shared rsIDs leave the database (issue \#42).
Any other pairing – data frames, a data frame with a lazy table, tables
on different connections, Arrow queries – has its rsID columns collected
first and intersected in R. The set operations run on one-column tibbles
rather than vectors: on tibbles dplyr uses vctrs, about 6x faster than
base
[`intersect()`](https://rdrr.io/r/base/sets.html)/[`setdiff()`](https://rdrr.io/r/base/sets.html)
on 10M rsIDs; on bare character vectors dplyr's generics fall through to
base.

## Usage

``` r
shared_rsids(exposure, outcome, instrument_snps, verbose = FALSE)
```

## Arguments

- exposure:

  Data frame or lazy table with an `SNP` column.

- outcome:

  Data frame or lazy table with an `rsids` column.

- instrument_snps:

  Character vector of rsIDs to leave out.

- verbose:

  Logical. If `TRUE`, say when a lazy input's rsID column has to be
  collected because the two tables are not on one connection.

## Value

Character vector of unique, non-missing shared rsIDs, in no particular
order.
