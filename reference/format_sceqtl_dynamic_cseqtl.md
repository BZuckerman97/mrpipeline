# Format single-cell RNA eQTL data from dynamic_cseqtl for MR analysis

Reads and formats eQTL summary statistics from the dynamic cis
single-cell eQTL (dynamic_cseqtl) dataset for a single cell-type file.

## Usage

``` r
format_sceqtl_dynamic_cseqtl(
  file,
  cell_type = NULL,
  pval_thresh = 1e-05,
  mhc_region = NULL
)
```

## Arguments

- file:

  Character. Path to a dynamic_cseqtl MR file (e.g.,
  `"CD4T_500kb_combined.MR.tsv.gz"`).

- cell_type:

  Character or `NULL`. Cell type label used in the phenotype string
  (`gene___cell_type`). If `NULL`, derived from the filename by
  stripping `_500kb_combined.MR.tsv.gz`.

- pval_thresh:

  Numeric or `NULL`. Retain only variants with `pval < pval_thresh`.
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
Note: EAF is not available in dynamic_cseqtl files; `eaf.exposure` will
be absent and `exposure_n` must be supplied explicitly to
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Details

Input files are expected to be pre-formatted tab-separated files with
columns: `gene`, `SNP`, `beta`, `se`, `effect_allele`, `other_allele`,
`pval`, `chr`, `pos`. The naming convention for files is
`<cell_type>_500kb_combined.MR.tsv.gz`. The exposure phenotype label is
`gene___<cell_type>`.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)

Other sceqtl-format:
[`format_sceqtl_1m_scbloodnl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_1m_scbloodnl.md),
[`format_sceqtl_dice()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dice.md),
[`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)

## Examples

``` r
if (FALSE) { # \dontrun{
exposure <- format_sceqtl_dynamic_cseqtl(
  file = "CD4T_500kb_combined.MR.tsv.gz"
)
} # }
```
