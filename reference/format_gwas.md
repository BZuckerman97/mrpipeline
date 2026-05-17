# Format GWAS summary statistics to a standard schema

Normalises column names from any GWAS summary statistics file to the
canonical schema expected by
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).
Handles datasets where rsIDs are absent by looking them up from a PLINK
bim file, and supports extracting chromosome and position from compound
marker ID columns (e.g. SCALLOP `MarkerName` format: `"CHR:POS:A1_A2"`).

## Usage

``` r
format_gwas(
  path,
  phenotype_id,
  type = c("outcome", "exposure"),
  col_map = NULL,
  bim_path = NULL,
  marker_col = NULL,
  marker_sep = ":",
  log10_pval = FALSE,
  flip_beta = FALSE,
  n = NULL
)
```

## Arguments

- path:

  Character file path (`.tsv`, `.tsv.gz`, `.txt.gz`, etc. —
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  auto-detects compression) or a pre-loaded data frame.

- phenotype_id:

  Character. Trait / phenotype identifier (e.g. `"IL-18"`, `"CAD"`).
  Stored in the `phenotype` output column.

- type:

  `"outcome"` (default) or `"exposure"`. `"outcome"` returns the
  normalised data frame ready for
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)'s
  `outcome` argument. `"exposure"` additionally calls
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html)
  and returns TwoSampleMR-formatted data with `exposure.` column
  suffixes.

- col_map:

  Named list of extra column-name aliases, e.g.
  `list(pval = "PVALUE", n = "SampleSize")`. Only needed for column
  names not already covered by the built-in alias table (see *Column
  normalisation* section). User entries are checked before the built-in
  list.

- bim_path:

  Character. Path to a PLINK bfile prefix (without `.bim`) used to
  recover rsIDs when the data lacks them. Required whenever the rsids
  column is absent.

- marker_col:

  Character. Name of a compound marker ID column in
  `"CHR<sep>POS<sep>..."` format (e.g. `"MarkerName"` for SCALLOP
  files). When supplied, `chr` and `pos` are parsed from this column.

- marker_sep:

  Character. Field separator used in `marker_col`. Default `":"`.

- log10_pval:

  Logical. If `TRUE`, the p-value column is in `-log10` scale and is
  back-transformed via `10^-x`. Default `FALSE`.

- flip_beta:

  Logical. If `TRUE`, multiplies `beta` by `-1` — use when the source
  file encodes the inverse direction of the intended exposure (e.g.
  modelling NLRP3 activation rather than suppression). Default `FALSE`.

- n:

  Integer. Explicit sample size. Added as the `n` column only when no
  sample-size column is already present in the data.

## Value

- `type = "outcome"`: a data frame with columns `rsids`, `chr`, `pos`,
  `beta`, `se`, `eaf`, `pval`, `n`, `effect_allele`, `other_allele`,
  `phenotype` (plus any extra columns from the source file).

- `type = "exposure"`: a TwoSampleMR-formatted data frame (output of
  [`TwoSampleMR::format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.html))
  with `exposure.`-suffixed columns, suitable for
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)'s
  `exposure` argument.

## Column normalisation

The function renames source columns to a fixed canonical schema by
checking a built-in table of known aliases for each target column:

|  |  |
|----|----|
| Canonical | Built-in aliases recognised automatically |
| `rsids` | `rsid`, `rs_id`, `rsID`, `SNP` |
| `chr` | `chromosome`, `Chr`, `CHROM`, `#CHROM`, `CHR` |
| `pos` | `base_pair_location`, `PosB37`, `PosB38`, `BP`, `POS`, `position`, `GENPOS` |
| `beta` | `Beta`, `Effect`, `BETA` |
| `se` | `standard_error`, `StdErr`, `SE`, `sebeta` |
| `eaf` | `effect_allele_frequency`, `Freq1`, `EAFrq`, `A1FREQ`, `af_alt`, `EAF` |
| `pval` | `p_value`, `P-value`, `P`, `Pval`, `p.value` |
| `n` | `N`, `TotalSampleSize`, `n_total` |
| `effect_allele` | `Allele1`, `EA`, `A1`, `ALLELE1`, `effectAllele`, `ALT` |
| `other_allele` | `Allele2`, `OA`, `A2`, `ALLELE0`, `otherAllele`, `REF` |

Supply `col_map` **only** when your dataset uses a column name that does
not appear in the table above — for example, if your p-value column is
called `"PVALUE"`, add `col_map = list(pval = "PVALUE")`. Inspect
[`names()`](https://rdrr.io/r/base/names.html) of your loaded data to
check. User-supplied aliases are checked before the built-in list, so
they take precedence in the event of ambiguity.

## rsID lookup from bim file

When `rsids` is absent (or all `NA`) after column normalisation, and
`bim_path` is supplied, the function inner-joins the data to the PLINK
bim file by chromosome and position to recover rsIDs. Rows without a bim
match are dropped — they are absent from the reference panel and cannot
be used in LD-based analyses. A message reports how many SNPs were
retained.

## Marker column parsing

Set `marker_col` to the name of a compound marker ID column whose values
have the form `"CHR<sep>POS<sep>..."` (e.g. SCALLOP `"MarkerName"`).
`chr` and `pos` are extracted from the first two fields. This step runs
before the rsID lookup so that the extracted coordinates are available
for the bim join.

## Examples

``` r
if (FALSE) { # \dontrun{
# Outcome GWAS whose columns are already in the built-in alias table
cad <- format_gwas(
  path         = "genomics_data/outcome_GWAS/CAD/cad_gwas.tsv.gz",
  phenotype_id = "CAD"
)

# SCALLOP outcome: no rsIDs, chr+pos embedded in MarkerName column
scallop_il6 <- format_gwas(
  path         = "genomics_data/outcome_GWAS/SCALLOP/CVD1_IL6.tsv.gz",
  phenotype_id = "IL6",
  marker_col   = "MarkerName",
  bim_path     = "LD_ref/g1000_eur"
)

# Dataset with a non-standard p-value column not in the alias table
ebi_il18 <- format_gwas(
  path         = "genomics_data/outcome_GWAS/EBI/GCST90428399.tsv.gz",
  phenotype_id = "IL-18",
  col_map      = list(pval = "PVALUE")
)

# Exposure GWAS — flip beta to model NLRP3 activation not suppression
exposure <- format_gwas(
  path         = "NLRP3/Output/NLRP3_CRP_IVs_300kb.tsv",
  phenotype_id = "NLRP3",
  type         = "exposure",
  flip_beta    = TRUE
)
} # }
```
