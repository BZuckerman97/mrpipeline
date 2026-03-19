# Package index

## MR Analysis

Run Mendelian randomisation and inspect results.

- [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  : Perform Mendelian randomisation analysis
- [`print(`*`<mr_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/print.mr_result.md)
  : Print an mr_result object
- [`summary(`*`<mr_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/summary.mr_result.md)
  : Summarise an mr_result object
- [`plot(`*`<mr_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/plot.mr_result.md)
  : Plot an mr_result object

## Colocalization

Run colocalization analysis and inspect results.

- [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
  : Perform colocalization analysis
- [`print(`*`<coloc_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/print.coloc_result.md)
  : Print a coloc_result object
- [`summary(`*`<coloc_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/summary.coloc_result.md)
  : Summarise a coloc_result object
- [`plot(`*`<coloc_result>`*`)`](https://github.com/BZuckerman97/mrpipeline/reference/plot.coloc_result.md)
  : Plot a coloc_result object

## Data Formatting

Format proteomics GWAS data for use with mrpipeline.

- [`format_pqtl_decode()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_decode.md)
  : Formats deCODE genetics proteomic data for Mendelian Randomization
  analysis
- [`format_pqtl_ukbppp()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_ukbppp.md)
  : Formats pQTL data from UKB-PPP for MR analysis.
- [`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)
  : Format single-cell RNA eQTL data for analysis
- [`decode_pqtl_file_name()`](https://github.com/BZuckerman97/mrpipeline/reference/decode_pqtl_file_name.md)
  : DECODE_PQTL_FILE_NAME
- [`ukbppp_pqtl_file_name()`](https://github.com/BZuckerman97/mrpipeline/reference/ukbppp_pqtl_file_name.md)
  : UKBPPP_PQTL_FILE_NAME

## Utilities

Helper functions for analysis workflows.

- [`get_gene_coords()`](https://github.com/BZuckerman97/mrpipeline/reference/get_gene_coords.md)
  : Look up gene coordinates from Ensembl via biomaRt
