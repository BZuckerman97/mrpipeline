# LD (r2) of every SNP to a specified index SNP, via local reference panel

Reuses
[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md)
– no separate LD mechanism, no network calls, and consistent sign/allele
handling with the rest of the package.
[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md)'s
matrix is signed (`--r`, not `--r2`) but squaring removes the sign,
which is all a colocalization regional plot's LD-colour scale needs.

## Usage

``` r
compute_ld_to_index(
  snps,
  index_snp,
  bfile,
  plink_bin = NULL,
  plink_threads = NULL,
  plink_memory = NULL
)
```

## Arguments

- snps:

  Character vector of rsIDs.

- index_snp:

  rsID of the SNP to compute LD against. Must be one of `snps`.

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

A data frame with columns `SNP` and `r2`.
