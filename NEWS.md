# mrpipeline 0.10.0

## New features

* `plot(type = "regional")` for `coloc_result` objects colours both panels by
  LD (r2) with an index SNP and marks that SNP, with SNPs that have no r2 in
  grey. `run_coloc()` now stores the r2 when given a `bfile`, so plotting needs
  no PLINK or LD panel. It uses the exposure lead, with ties in `pval.exposure`
  going to the largest |beta/se|, or the SNP given by the new `index_snp`
  argument to `run_coloc()`, which affects only plotting. `index_snp` in
  `plot()` overrides it for one plot (r2 is then computed, so `bfile` is
  needed). Results with no stored r2, from
  `run_coloc()` with no `bfile` or saved by an earlier version, still plot:
  they are uncoloured unless `bfile` is passed to `plot()`.
* `plot(type = "regional", ens_db = )` adds a gene track under the panels and
  plots position in Mb. Without `ens_db` the plot needs no `locuszoomr` or
  `ensembldb`, and stays a plain ggplot. `cowplot` is added to Suggests for the
  gene track.
* `coloc_result` gains `ld_index_snp`, and on a successful run with a `bfile`,
  `harmonised_data` gains a column `ld_r2_index`. Coloc results are unchanged.

## Deprecated

* `plot(type = "locuszoom")` for `coloc_result` objects is deprecated in favour
  of `plot(type = "regional", ens_db = )`, which returns a ggplot. It still
  works, and now warns.

# mrpipeline 0.9.0

## New features

* `run_coloc()` no longer requires `bfile` for `methods = "abf"`, which uses
  no LD. With no `bfile` it runs on the full harmonised SNP set and PLINK is
  not needed; `"susie"`, `"signals"` and `"prop_test"` still need `bfile`, and
  the error names them. Supplying a `bfile` is unchanged: every method,
  `"abf"` included, runs on the SNPs aligned to the panel. Without `bfile`,
  `plot(type = "locuszoom")` now warns and draws the plots with no LD
  colouring, and SNPs absent from the panel are kept rather than dropped
  (#48).

# mrpipeline 0.8.0

## New features

* `run_mr()` and `run_coloc()` accept a lazy dplyr table for `exposure` and
  `outcome` as well as a data frame: a `dplyr::tbl()` on a DuckDB
  connection (via dbplyr), or any backend that supports `filter()` and
  `collect()`, such as an Arrow dataset query. Each function filters in the
  backend and collects only the rows it uses: the instruments, the outcome
  at those rsIDs, about 1000 SNPs for the allele orientation check, or the
  colocalization window. A genome-wide GWAS can therefore stay on disk
  instead of being loaded into memory on every call. Results are identical
  to the data-frame path, and a test pins this on the bundled data (#42).

* For two lazy tables on the same database connection, the allele
  orientation check intersects the exposure and outcome rsIDs inside the
  database. Tables on different connections, or a data frame paired with a
  lazy table, still work, but each rsID column is collected into R for the
  intersect. `verbose = TRUE` reports when that happens.

* DBI, dbplyr and duckdb are suggested packages only. mrpipeline installs
  and runs on data frames without them. Opening the database is up to the
  caller; see *Large GWAS: lazy tables* in the user guide.

## Performance

* The allele orientation check intersects the shared rsIDs with
  `dplyr::intersect()`/`dplyr::setdiff()` on one-column tibbles rather than
  base `intersect()`/`setdiff()` on vectors. That is about 6x faster for
  data-frame inputs too (vctrs set operations; benchmarked on 10M rsIDs).

## Bug fixes

* Seeded analyses no longer depend on whether a call is the first in the
  session. Loading ieugwasr draws from the random-number generator, and
  ieugwasr used to load inside the first `run_mr()` call, after the
  caller's `set.seed()`. That call's weighted-median bootstrap SE and
  TwoSampleMR's `id.outcome` therefore differed from an identical later
  call. mrpipeline now loads ieugwasr when it loads, and puts the RNG state
  back as it was. withr moves from an indirect dependency (via dplyr's
  tidyselect) to Imports (#43).
