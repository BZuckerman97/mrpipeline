# Align harmonised data to an LD matrix

Subsets both the harmonised data frame and LD matrix to their shared
SNPs, and reorders both to match.

## Usage

``` r
align_to_ld_matrix(harmonised_data, ld_matrix)
```

## Arguments

- harmonised_data:

  Data frame with a `SNP` column.

- ld_matrix:

  Square matrix with rsID row/column names.

## Value

A named list with elements:

- `data`: subset and reordered data frame

- `ld_matrix`: subset and reordered matrix
