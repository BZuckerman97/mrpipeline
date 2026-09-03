# Format single-cell RNA eQTL data from OneK1K for MR analysis

Reads and formats eQTL summary statistics from the OneK1K project for a
specific cell type, producing a data frame ready for
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
or
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Usage

``` r
format_single_cell_onek1k(
  onek1k_mapping,
  onek1k_cell_type,
  sample_size = NULL,
  pval_thresh = 1e-05,
  mhc_region = NULL
)
```

## Arguments

- onek1k_mapping:

  A data frame with columns `cell_type` and `path_to_eqtl_file`, mapping
  each cell type to its eQTL summary statistics file.

- onek1k_cell_type:

  Character. The cell type to format (e.g., `"cd4nc"`). Must match an
  entry in the `cell_type` column of `onek1k_mapping`. When the file
  contains a `CELL_ID` column, rows are filtered to this value.

- sample_size:

  Integer or `NULL`. Sample size for this cell type, stored as
  `samplesize.exposure`. Required by
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md);
  if `NULL` (default), no sample size column is added.

- pval_thresh:

  Numeric or `NULL`. Retain only variants with `P_VALUE < pval_thresh`.
  Default `1e-5`. Set `NULL` to skip filtering.

- mhc_region:

  Data frame with columns `chr`, `start`, `end` defining a genomic
  region to exclude (typically the MHC locus). Default `NULL` (no
  exclusion). Use `data.frame(chr = "6", start = 25e6, end = 34e6)` for
  the standard MHC window used in this project.

## Value

A data frame formatted as a TwoSampleMR exposure dataset, including
`chr.exposure`, `pos.exposure`, and `eaf.exposure` columns required by
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Details

OneK1K data is available from <https://onek1k.org/>. The expected input
format is the per-cell-type or combined eQTL table with columns: `CHR`,
`POS`, `RSID`, `A1`, `A2`, `A2_FREQ_ONEK1K`, `SPEARMANS_RHO`, `P_VALUE`,
`GENE`, and optionally `CELL_ID`.

`A2` is the effect allele (the allele whose frequency is reported in
`A2_FREQ_ONEK1K`) and `A1` is the other allele. Spearman's rho
(`SPEARMANS_RHO`) is used as the effect size (beta). Standard errors are
derived from the p-value:
`se = |rho| / qnorm(P_VALUE / 2, lower.tail = FALSE)`. Rows with
non-finite or zero SE are dropped. The exposure phenotype label is
`GENE___<cell_type>`.

If the file contains a `CELL_ID` column (combined multi-cell-type file),
rows are filtered to those matching `onek1k_cell_type`.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)

Other sceqtl-format:
[`format_sceqtl_1m_scbloodnl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_1m_scbloodnl.md),
[`format_sceqtl_dice()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dice.md),
[`format_sceqtl_dynamic_cseqtl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dynamic_cseqtl.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mapping_df <- data.frame(
  cell_type         = "cd4nc",
  path_to_eqtl_file = "esnp_table.tsv.gz"
)
exposure <- format_single_cell_onek1k(
  onek1k_mapping   = mapping_df,
  onek1k_cell_type = "cd4nc",
  sample_size      = 463528L
)
} # }
```
