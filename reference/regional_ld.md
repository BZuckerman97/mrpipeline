# r2 and index SNP for a regional plot

The index SNP is `index_snp`, or else the one stored by
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md),
or else the exposure lead. Stored r2 is used when it is relative to that
SNP, and otherwise r2 is computed from `bfile`.

## Usage

``` r
regional_ld(x, index_snp = NULL, bfile = NULL, plink_bin = NULL, ...)
```

## Arguments

- x:

  A `coloc_result` object.

- index_snp:

  rsID to colour LD against, overriding any index SNP stored by
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
  for this plot. Defaults to the stored index SNP, or the exposure lead
  when there is none: the lowest `pval.exposure`, with ties going to the
  largest \|beta/se\|. Must be one of the SNPs in `x`.

- bfile:

  Path to PLINK bfile prefix (without .bed/.bim/.fam), the local LD
  reference panel. Only used when r2 has to be computed: `x` has no
  stored r2, or `index_snp` differs from the stored index SNP. SNPs
  absent from the panel stay on the plot in grey.

- plink_bin:

  Path to PLINK binary. If `NULL`, auto-detected via
  [`genetics.binaRies::get_plink_binary()`](https://rdrr.io/pkg/genetics.binaRies/man/get_plink_binary.html).

- ...:

  Forwarded to
  [`compute_ld_to_index()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_to_index.md)
  (`plink_threads`, `plink_memory`).

## Value

`NULL` for no LD colouring, or a list with `index_snp` and `r2`, aligned
to the rows of `harmonised_data`. `r2` is `NULL` when an index SNP is
recorded but there is no r2 for it and no `bfile`.
