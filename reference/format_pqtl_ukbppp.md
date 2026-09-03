# Formats pQTL data from UKB-PPP for MR analysis.

This function processes pQTL summary statistics from UKB-PPP, matches
them with rsID information, handles non-Mendelian chromosomes,
standardizes column names, and formats the data for use with the
TwoSampleMR package.

## Usage

``` r
format_pqtl_ukbppp(
  ukbppp,
  ukbppp_rsid,
  pqtl_assay,
  x_y_chr_file = NULL,
  pos_build = c("b38", "b37"),
  type = c("exposure", "outcome")
)
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

- pos_build:

  Character, genome build to use for position matching against the rsid
  file. Either `"b38"` (GRCh38, uses `POS38` column; default) or `"b37"`
  (GRCh37, uses `POS19` column). Must match the build used in the
  UKB-PPP GWAS `GENPOS` column.

- type:

  Character, either `"exposure"` (default) or `"outcome"`. `"exposure"`
  runs the data through
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html).
  `"outcome"` skips `format_data()` and returns a plain data frame
  normalised to the same schema as `format_gwas(type = "outcome")`, so
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  can pre-filter by rsids before calling `format_data()` internally.

## Value

A plain data frame (not a list) in both cases: for `type = "exposure"`,
the output of TwoSampleMR::format_data(); for `type = "outcome"`, data
normalised to `format_gwas(type = "outcome")`'s schema.

## Examples

``` r
# See the test script for example usage.
```
