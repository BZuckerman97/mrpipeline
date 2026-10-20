# Changelog

## mrpipeline 0.8.0

### New features

- [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  and
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
  accept a lazy dplyr table for `exposure` and `outcome` as well as a
  data frame: a
  [`dplyr::tbl()`](https://dplyr.tidyverse.org/reference/tbl.html) on a
  DuckDB connection (via dbplyr), or any backend that supports
  [`filter()`](https://rdrr.io/pkg/ensembldb/man/global-filters.html)
  and `collect()`, such as an Arrow dataset query. Each function filters
  in the backend and collects only the rows it uses: the instruments,
  the outcome at those rsIDs, about 1000 SNPs for the allele orientation
  check, or the colocalization window. A genome-wide GWAS can therefore
  stay on disk instead of being loaded into memory on every call.
  Results are identical to the data-frame path, and a test pins this on
  the bundled data (#42).

- For two lazy tables on the same database connection, the allele
  orientation check intersects the exposure and outcome rsIDs inside the
  database. Tables on different connections, or a data frame paired with
  a lazy table, still work, but each rsID column is collected into R for
  the intersect. `verbose = TRUE` reports when that happens.

- DBI, dbplyr and duckdb are suggested packages only. mrpipeline
  installs and runs on data frames without them. Opening the database is
  up to the caller; see *Large GWAS: lazy tables* in the user guide.

### Performance

- The allele orientation check intersects the shared rsIDs with
  [`dplyr::intersect()`](https://generics.r-lib.org/reference/setops.html)/[`dplyr::setdiff()`](https://generics.r-lib.org/reference/setops.html)
  on one-column tibbles rather than base
  [`intersect()`](https://rdrr.io/r/base/sets.html)/[`setdiff()`](https://rdrr.io/r/base/sets.html)
  on vectors. That is about 6x faster for data-frame inputs too (vctrs
  set operations; benchmarked on 10M rsIDs).

### Bug fixes

- Seeded analyses no longer depend on whether a call is the first in the
  session. Loading ieugwasr draws from the random-number generator, and
  ieugwasr used to load inside the first
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  call, after the caller’s
  [`set.seed()`](https://rdrr.io/r/base/Random.html). That call’s
  weighted-median bootstrap SE and TwoSampleMR’s `id.outcome` therefore
  differed from an identical later call. mrpipeline now loads ieugwasr
  when it loads, and puts the RNG state back as it was. withr moves from
  an indirect dependency (via dplyr’s tidyselect) to Imports (#43).
