# Compute LD correlation matrix from a local reference panel

Calls
[`ieugwasr::ld_matrix()`](https://mrcieu.github.io/ieugwasr/reference/ld_matrix.html),
which (for a local `bfile`) runs `plink --r square --keep-allele-order`
and names rows/columns `rsid_A1_A2`, where A1/A2 are the reference
panel's own `.bim` allele coding. This is a **signed** correlation
matrix (`--r`, not `--r2`) and the sign of every entry is anchored to
that specific, arbitrary A1 allele – which need not match the effect
allele used anywhere else in the pipeline. That allele identity is only
recoverable from the `rsid_A1_A2` name suffix, so it is parsed out here
and returned alongside the matrix (rather than discarded) so
[`align_to_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/align_to_ld_matrix.md)
can re-orient betas/z-scores to this matrix's sign convention before
they are used together in
[`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html)/[`susieR::susie_rss()`](https://rdrr.io/pkg/susieR/man/susie_rss.html).
See
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)'s
"LD matrix orientation" section for the full rationale.

## Usage

``` r
compute_ld_matrix(
  snps,
  bfile,
  plink_bin = NULL,
  plink_threads = NULL,
  plink_memory = NULL
)
```

## Arguments

- snps:

  Character vector of rsIDs.

- bfile:

  Path to PLINK bfile prefix (without .bed/.bim/.fam).

- plink_bin:

  Path to PLINK binary. If `NULL`, auto-detected via
  [`genetics.binaRies::get_plink_binary()`](https://rdrr.io/pkg/genetics.binaRies/man/get_plink_binary.html).

- plink_threads:

  Number of threads for PLINK. `NULL` = auto-detect.

- plink_memory:

  Memory (MB) for PLINK. `NULL` = auto-detect.

## Value

A named list:

- `ld`: square signed correlation matrix with rsID-only row/column names

- `alleles`: data frame with columns `SNP`, `ld_a1`, `ld_a2` giving the
  reference panel's allele coding for each SNP in `ld` – the allele that
  `ld`'s sign is anchored to (`ld_a1`) and its counterpart
