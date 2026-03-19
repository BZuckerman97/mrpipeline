# Formats pQTL data from UKB-PPP for MR analysis.

This function processes pQTL summary statistics from UKB-PPP, matches
them with rsID information, handles non-Mendelian chromosomes,
standardizes column names, and formats the data for use with the
TwoSampleMR package.

## Usage

``` r
format_pqtl_ukbppp(ukbppp, ukbppp_rsid, pqtl_assay, x_y_chr_file = NULL)
```

## Arguments

- ukbppp:

  Dataframe, file path to the file containing the ukbppp GWAS data

- ukbppp_rsid:

  Dataframe, file path to the file containing the ukbppp GWAS data rsids

- pqtl_assay:

  String, of the ukbppp protein assay

- x_y_chr_file:

  Data frame or string file path containing rsIDs for X and Y
  chromosomes, or NULL. When a data frame is supplied, it is used
  directly (skipping `fread()`). When a string path is supplied, the
  file is read via
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).

## Value

A list with two elements:

- `exposure`: Formatted exposure data frame (output of
  TwoSampleMR::format_data).

## Examples

``` r
# See the test script for example usage.
```
