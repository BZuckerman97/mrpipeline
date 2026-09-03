# LD-coloured regional plots for coloc results, via locuszoomr

Builds exposure and outcome regional plots using
[`locuszoomr::locus()`](https://rdrr.io/pkg/locuszoomr/man/locus.html),
coloured by LD (r2) to a chosen index SNP, and stacks them with
[`locuszoomr::multi_layout()`](https://rdrr.io/pkg/locuszoomr/man/multi_layout.html).
LD is computed from a local reference panel via
[`compute_ld_to_index()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_to_index.md)
(reusing
[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md)
– no network lookups).

## Usage

``` r
plot_coloc_locuszoom(x, ens_db, bfile, plink_bin = NULL, index_snp = NULL, ...)
```

## Arguments

- x:

  A `coloc_result` object.

- ens_db:

  Either an
  [`ensembldb::EnsDb`](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)
  object, or the name of an Ensembl annotation package as a string (e.g.
  `"EnsDb.Hsapiens.v75"`), matching the genome build of `x`'s data. If
  given as a string, the package must be **attached** first
  ([`library(EnsDb.Hsapiens.v75)`](https://rdrr.io/r/base/library.html)),
  not merely installed –
  [`locuszoomr::locus()`](https://rdrr.io/pkg/locuszoomr/man/locus.html)
  resolves a character `ens_db` via `get(ens_db)` on the search path,
  and errors ("Ensembl database not loaded") if it isn't attached.

- bfile:

  Path to PLINK bfile prefix (without .bed/.bim/.fam), used as the local
  LD reference panel.

- plink_bin:

  Path to PLINK binary. If `NULL`, auto-detected via
  [`genetics.binaRies::get_plink_binary()`](https://rdrr.io/pkg/genetics.binaRies/man/get_plink_binary.html).

- index_snp:

  rsID to colour LD against and centre the plotted region on. Defaults
  to the SNP with the lowest `pval.exposure`, matching
  [`locuszoomr::locus()`](https://rdrr.io/pkg/locuszoomr/man/locus.html)'s
  own default when no index SNP is given.

- ...:

  Forwarded to
  [`compute_ld_to_index()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_to_index.md)
  (`plink_threads`, `plink_memory`).

## Value

`NULL`, invisibly. Like base/grid plotting functions (and unlike this
package's other `plot_coloc_*()` functions, which return a ggplot
object), this draws directly to the current graphics device as a side
effect – open a device ([`pdf()`](https://rdrr.io/r/grDevices/pdf.html),
[`png()`](https://rdrr.io/r/grDevices/png.html), ...) before calling,
and [`dev.off()`](https://rdrr.io/r/grDevices/dev.html) after, to save
the result.

## Details

There is no genome-build field stored on a `coloc_result` – `ens_db`
must match whatever build the exposure/outcome data (and `bfile`)
actually use (e.g. `EnsDb.Hsapiens.v75` for GRCh37, `EnsDb.Hsapiens.v86`
for GRCh38).
