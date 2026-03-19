# Formats deCODE genetics proteomic data for Mendelian Randomization analysis

This function filters deCODE genetics GWAS by the deCODE included
variants file, processes the data to clean the column headings and uses
TwoSampleMR::format_data() to prepare the exposure dataset before
run_mr()

## Usage

``` r
format_pqtl_decode(
  decode_proteomic_gwas_file_path,
  decode_included_variants_file_path,
  pqtl_assay,
  x_y_chr_file = NULL
)
```

## Arguments

- decode_proteomic_gwas_file_path:

  Character vector of file path(s) to deCODE GWAS data, or a single
  pre-loaded data.frame. If multiple paths are provided, data will be
  combined.

- decode_included_variants_file_path:

  Character vector of file path(s) to the deCODE included variants data,
  or a single pre-loaded data.frame. If multiple paths are provided,
  data will be combined.

- pqtl_assay:

  String, name of the deCODE genetics protein assayed

- x_y_chr_file:

  String, optional file path to a file containing rsids for X and Y
  chromosomes. This file should be tab-separated with columns:
  Chromosome, Position, RSID.

## Value

A list with two elements:

- `exposure`: Formatted exposure data frame (output of
  TwoSampleMR::format_data).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- format_pqtl_decode(
  decode_proteomic_gwas_file_path = "path/to/decode_gwas.txt.gz",
  decode_included_variants_file_path = "path/to/included_variants.txt.gz",
  pqtl_assay = "CD40"
)
head(result$exposure)
} # }
```
