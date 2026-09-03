# Format single-cell RNA eQTL data from DICE for MR analysis

Reads and formats eQTL summary statistics from the Database of Immune
Cell Expression, eQTLs and Epigenomics (DICE; Schmiedel et al. 2018) for
a single cell-type VCF file.

## Usage

``` r
format_sceqtl_dice(
  file,
  cell_type = NULL,
  pval_thresh = 1e-05,
  mhc_region = NULL
)
```

## Arguments

- file:

  Character. Path to a DICE eQTL VCF file (e.g.,
  `"t_cell_cd4_naive.vcf"`). Plain and gzip-compressed files are both
  supported.

- cell_type:

  Character or `NULL`. Cell type label used in the phenotype string
  (`GENE___cell_type`). If `NULL`, derived from the filename by
  stripping the `.vcf` (or `.vcf.gz`) extension.

- pval_thresh:

  Numeric or `NULL`. Retain only variants with `Pvalue < pval_thresh`.
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
Note: EAF is not available in DICE VCF files; `eaf.exposure` will be
absent and `exposure_n` must be supplied explicitly to
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Details

The VCF INFO field is parsed with regex to extract `Gene`, `GeneSymbol`,
`Pvalue`, and `Beta`. SE is derived as
`se = |Beta| / sqrt(qchisq(Pvalue, df = 1, lower.tail = FALSE))`. Only
biallelic SNPs (single-nucleotide, A/C/G/T) are retained. Palindromic
SNPs (A/T and C/G pairs) are dropped because DICE provides no allele
frequency data to resolve strand ambiguity. `ALT` is treated as the
effect allele and `REF` as the other allele. The exposure phenotype
label is `GeneSymbol___<cell_type>`, falling back to `Gene` when
`GeneSymbol` is absent.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)

Other sceqtl-format:
[`format_sceqtl_1m_scbloodnl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_1m_scbloodnl.md),
[`format_sceqtl_dynamic_cseqtl()`](https://github.com/BZuckerman97/mrpipeline/reference/format_sceqtl_dynamic_cseqtl.md),
[`format_single_cell_onek1k()`](https://github.com/BZuckerman97/mrpipeline/reference/format_single_cell_onek1k.md)

## Examples

``` r
if (FALSE) { # \dontrun{
exposure <- format_sceqtl_dice(file = "t_cell_cd4_naive.vcf")
} # }
```
