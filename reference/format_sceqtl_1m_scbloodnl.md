# Format single-cell RNA eQTL data from 1M-scBloodNL for MR analysis

Reads and formats eQTL summary statistics from the 1 Million Immune
Cells single-cell blood cohort (1M-scBloodNL; Oelen et al. 2022) for a
single cell-type file.

## Usage

``` r
format_sceqtl_1m_scbloodnl(
  file,
  cell_type = NULL,
  cis_only = TRUE,
  pval_thresh = 1e-05,
  mhc_region = NULL
)
```

## Arguments

- file:

  Character. Path to a 1M-scBloodNL per-cell-type eQTL file (e.g.,
  `"CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz"`).

- cell_type:

  Character or `NULL`. Cell type label used in the phenotype string
  (`GENE___cell_type`). If `NULL`, derived from the filename by
  stripping `_expression_eQTLsFDR-ProbeLevel.txt.gz`.

- cis_only:

  Logical. If `TRUE` (default), retain only rows where
  `CisTrans == "cis"`.

- pval_thresh:

  Numeric or `NULL`. Retain only variants with `PValue < pval_thresh`.
  Default `1e-5`. Set `NULL` to skip filtering.

- mhc_region:

  Data frame with columns `chr`, `start`, `end` defining a genomic
  region to exclude (typically the MHC locus). Default `NULL` (no
  exclusion). Use `data.frame(chr = "6", start = 25e6, end = 34e6)` for
  the standard MHC window used in this project.

## Value

A data frame formatted as a TwoSampleMR exposure dataset, including
`chr.exposure` and `pos.exposure` columns required by
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).
Note: EAF is not available in 1M-scBloodNL files; `eaf.exposure` will be
absent and `exposure_n` must be supplied explicitly to
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Details

Beta and SE are derived from the Z-score and total sample size:
`se = 1 / sqrt(N_total)`, `beta = OverallZScore * se`, where `N_total`
is the sum of per-dataset sample counts from the `DatasetsNrSamples`
column (semicolon-separated integers). The other allele is derived by
splitting `SNPType` on `/` and taking the allele not equal to
`AlleleAssessed`. Rows where the other allele cannot be unambiguously
determined are dropped. The exposure phenotype label is
`HGNCName___<cell_type>`.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)

Other sceqtl-format:
[`format_sceqtl_dice()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dice.md),
[`format_sceqtl_dynamic_cseqtl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dynamic_cseqtl.md),
[`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)

## Examples

``` r
if (FALSE) { # \dontrun{
exposure <- format_sceqtl_1m_scbloodnl(
  file = "CD4T_expression_eQTLsFDR-ProbeLevel.txt.gz"
)
} # }
```
