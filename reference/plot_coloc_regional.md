# Regional association plot for coloc results

Plots -log10(p) against position, with exposure and outcome panels.
Points are coloured by LD (r2) with an index SNP, using the r2 stored by
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md),
so no PLINK or panel is needed at plot time. With no stored r2 (a
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
call with no `bfile`, or an object saved before mrpipeline 0.10.0),
passing `bfile` computes it for this plot, and otherwise the plot is
uncoloured, with any recorded index SNP still marked.

## Usage

``` r
plot_coloc_regional(
  x,
  ens_db = NULL,
  bfile = NULL,
  plink_bin = NULL,
  index_snp = NULL,
  ...
)
```

## Arguments

- x:

  A `coloc_result` object.

- ens_db:

  `NULL` (the default) for no gene track, or an
  [`ensembldb::EnsDb`](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)
  object, or the name of an Ensembl annotation package as a string (e.g.
  `"EnsDb.Hsapiens.v75"`), matching the genome build of `x`'s data. If
  given as a string, the package must be **attached** first
  ([`library(EnsDb.Hsapiens.v75)`](https://rdrr.io/r/base/library.html)),
  not merely installed.

- bfile:

  Path to PLINK bfile prefix (without .bed/.bim/.fam), the local LD
  reference panel. Only used when r2 has to be computed: `x` has no
  stored r2, or `index_snp` differs from the stored index SNP. SNPs
  absent from the panel stay on the plot in grey.

- plink_bin:

  Path to PLINK binary. If `NULL`, auto-detected via
  [`genetics.binaRies::get_plink_binary()`](https://rdrr.io/pkg/genetics.binaRies/man/get_plink_binary.html).

- index_snp:

  rsID to colour LD against, overriding any index SNP stored by
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
  for this plot. Defaults to the stored index SNP, or the exposure lead
  when there is none: the lowest `pval.exposure`, with ties going to the
  largest \|beta/se\|. Must be one of the SNPs in `x`.

- ...:

  Forwarded to
  [`compute_ld_to_index()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_to_index.md)
  (`plink_threads`, `plink_memory`).

## Value

A ggplot object, or `NULL` invisibly if there is no harmonised data.

## Details

The same index SNP and r2 are used in both panels, so the outcome panel
shows LD with the exposure lead by default, not with the outcome lead.

With `ens_db`, a gene track from `locuszoomr` is added under the panels,
the x axis is in Mb rather than bp, and the plots are composed with
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).
Adding `ggplot2::labs(title = ...)` still sets the title, but other
`ggplot2` additions, such as `theme()` changes to the panels, do not
reach a composed plot.
