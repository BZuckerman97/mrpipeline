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

## Plotting

Plot MR results across multiple exposures/outcomes.

- [`forest_plot()`](https://github.com/BZuckerman97/mrpipeline/reference/forest_plot.md)
  : Forest plot of MR methods for one or more results
- [`outcome_forest_plot()`](https://github.com/BZuckerman97/mrpipeline/reference/outcome_forest_plot.md)
  : Forest plot of MR outcomes, optionally coloured/shaped by a grouping
  column
- [`table_forest_plot()`](https://github.com/BZuckerman97/mrpipeline/reference/table_forest_plot.md)
  : Table-style forest plot (forestplot package backend)

## Data Formatting

Format GWAS summary statistics for use with mrpipeline.

- [`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
  : Format GWAS summary statistics to a standard schema
- [`format_pqtl_decode()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_decode.md)
  : Formats deCODE genetics proteomic data for Mendelian Randomization
  analysis
- [`format_pqtl_ukbppp()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_ukbppp.md)
  : Formats pQTL data from UKB-PPP for MR analysis.
- [`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)
  : Format single-cell RNA eQTL data from OneK1K for MR analysis
- [`format_sceqtl_1m_scbloodnl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_1m_scbloodnl.md)
  : Format single-cell RNA eQTL data from 1M-scBloodNL for MR analysis
- [`format_sceqtl_dice()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dice.md)
  : Format single-cell RNA eQTL data from DICE for MR analysis
- [`format_sceqtl_dynamic_cseqtl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dynamic_cseqtl.md)
  : Format single-cell RNA eQTL data from dynamic_cseqtl for MR analysis
- [`decode_pqtl_file_name()`](https://github.com/BZuckerman97/mrpipeline/reference/decode_pqtl_file_name.md)
  : DECODE_PQTL_FILE_NAME
- [`ukbppp_pqtl_file_name()`](https://github.com/BZuckerman97/mrpipeline/reference/ukbppp_pqtl_file_name.md)
  : UKBPPP_PQTL_FILE_NAME

## Utilities

Helper functions for analysis workflows.

- [`get_gene_coords()`](https://github.com/BZuckerman97/mrpipeline/reference/get_gene_coords.md)
  : Look up gene coordinates from Ensembl via biomaRt

## Example Data

Bundled datasets used in examples and vignettes.

- [`cd40_exposure`](https://github.com/BZuckerman97/mrpipeline/reference/cd40_exposure.md)
  : CD40 exposure data (TwoSampleMR format)
- [`sjogren_outcome`](https://github.com/BZuckerman97/mrpipeline/reference/sjogren_outcome.md)
  : Sjogren's disease outcome data (standardised format)
