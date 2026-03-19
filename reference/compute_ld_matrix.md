# Compute LD correlation matrix from a local reference panel

Calls
[`ieugwasr::ld_matrix()`](https://mrcieu.github.io/ieugwasr/reference/ld_matrix.html)
and strips allele suffixes (e.g. `rs123_A_G` -\> `rs123`) from row and
column names.

## Usage

``` r
compute_ld_matrix(snps, bfile, plink_bin = NULL)
```

## Arguments

- snps:

  Character vector of rsIDs.

- bfile:

  Path to PLINK bfile prefix (without .bed/.bim/.fam).

- plink_bin:

  Path to PLINK binary. If `NULL`, auto-detected via
  [`genetics.binaRies::get_plink_binary()`](https://rdrr.io/pkg/genetics.binaRies/man/get_plink_binary.html).

## Value

A square correlation matrix with rsID-only row/column names.
