# mrpipeline Developer Guide

## Package Architecture

## S3 Classes

### `mr_result`

Every
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
call returns an `mr_result` object — even when no instruments survive
filtering. The `status` field indicates success or failure:

- `"success"` — analysis completed normally
- `"no_instruments"` — no instruments survived
  filtering/clumping/exclusion
- `"no_harmonised_variants"` — harmonisation removed all variants

The `status_reason` field provides a human-readable explanation (e.g.
`"No significant instruments in cis region for 'PCSK9'"`).

[`print.mr_result()`](https://github.com/BZuckerman97/mrpipeline/reference/print.mr_result.md)
and
[`summary.mr_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.mr_result.md)
display status information for non-success results.

### `coloc_result`

Every
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
call returns a `coloc_result` object — even when analysis cannot proceed
(e.g. no SNPs in region). The object contains:

- `coloc_abf` — output of
  [`coloc::coloc.abf()`](https://rdrr.io/pkg/coloc/man/coloc.abf.html),
  or `NULL`
- `coloc_susie` — output of
  [`coloc::coloc.susie()`](https://rdrr.io/pkg/coloc/man/coloc.susie.html),
  or `NULL`
- `coloc_signals` — output of
  [`coloc::coloc.signals()`](https://rdrr.io/pkg/coloc/man/coloc.signals.html),
  or `NULL`
- `coloc_prop_test` — output of `colocPropTest::coloc.prop.test()`, or
  `NULL`
- `n_snps` — integer, number of SNPs used in the analysis
- `harmonised_data` — data frame of harmonised data
- `methods_skipped` — named character vector (method name → reason
  skipped)
- `params` — list of all input parameters
- `status` — one of `"success"`, `"no_snps_in_region"`,
  `"too_few_snps"`, `"no_harmonised_variants"`
- `status_reason` — human-readable explanation when
  `status != "success"`

[`print.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/print.coloc_result.md)
shows a one-line summary: N SNPs, ABF PP.H4, and SuSiE max PP.H4 with
credible set count (if available).

[`summary.coloc_result()`](https://github.com/BZuckerman97/mrpipeline/reference/summary.coloc_result.md)
shows full details: all posterior probabilities (H0–H4), SuSiE credible
set pairs, signals hits, and skipped methods.

### Plot methods

Both S3 classes have
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods defined
in `R/plot.R`:

- `plot.mr_result(x, type)` — wraps TwoSampleMR plotting functions:
  - `"scatter"` (default):
    [`TwoSampleMR::mr_scatter_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_scatter_plot.html)
  - `"forest"`:
    [`TwoSampleMR::mr_forest_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_forest_plot.html)
  - `"funnel"`:
    [`TwoSampleMR::mr_funnel_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_funnel_plot.html)
- `plot.coloc_result(x, type)` — custom ggplot2 plots:
  - `"pp_bar"` (default): bar chart of ABF posterior probabilities
    (H0–H4)
  - `"regional"`: side-by-side regional association plots (-log10(p) vs
    position) for exposure and outcome

Both methods return `NULL` invisibly for non-success results with an
informative message.

## Exported Utility Functions

### `get_gene_coords()`

Defined in `R/get_gene_coords.R`. Queries Ensembl via biomaRt for gene
coordinates. Key design decisions:

- Uses
  [`biomaRt::useEnsembl()`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.html)
  with `host = "grch37.ensembl.org"` for GRCh37 or the default host for
  GRCh38
- Filters to standard chromosomes (1–22, X, Y)
- Deduplicates by gene + chromosome (widest range), preferring autosomes
- Warns about genes not found in Ensembl

## Internal Helper Functions

All helpers live in `R/helpers.R`, are tagged `@keywords internal`, and
are not exported.

### `harmonise_and_filter()`

### `compute_ld_matrix()`

### `clump_instruments()`

### `align_to_ld_matrix()`

### `eaf_to_maf()`, `resolve_sample_size()`

## Code Conventions

### Messages, Warnings, and Errors

Always use the `cli` package. Never use
[`message()`](https://rdrr.io/r/base/message.html),
[`warning()`](https://rdrr.io/r/base/warning.html), or
[`stop()`](https://rdrr.io/r/base/stop.html) directly.

``` r
cli::cli_inform("Processing {protein}...")      # informational (goes to stderr)
cli::cli_warn("Only {n} SNP(s) available.")     # warning
cli::cli_abort("bfile is required for coloc.")  # error (stops execution)
```

The `verbose` parameter is a **planned future feature**. When
implemented, all
[`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
calls should be gated behind it.

### Pipe

Use the native `|>` pipe exclusively. Never use `%>%` from magrittr.

### String Operations

Prefer `stringr` functions over base R equivalents:

``` r
# Good
stringr::str_remove(x, "_.*")
stringr::str_detect(x, "^chr")
stringr::str_replace(x, "chr", "")

# Avoid
gsub("_.*", "", x)
grepl("^chr", x)
sub("chr", "", x)
```

### Namespace

Prefer `pkg::fn()` over `@importFrom`. Only add `@importFrom` entries
when a function is called very frequently in a hot path.

### Documentation

- All exported functions: `@param`, `@return`, `@examples`, `@seealso`,
  `@family`
- Internal helpers: `@keywords internal` only; no `@export`
- When modifying an exported function, update its roxygen docs AND this
  vignette in the same commit

## Adding New MR Methods

## Adding New Coloc Methods

## Test Data

The package ships several bundled datasets (defined in `R/data.R`) and a
minimal LD reference panel. All SNP-level datasets overlap with the LD
panel, so integration tests can run without external downloads.

### Bundled datasets

| Dataset                | Rows | Format                     | Description                                                                                                                                                    |
|------------------------|------|----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `cd40_sumstats`        | 210  | Raw UKB-PPP                | CD40 summary statistics (chr6 MHC + chr20 cis), input for [`format_pqtl_ukbppp()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_ukbppp.md) |
| `cd40_exposure`        | 50   | TwoSampleMR exposure       | CD40 chr20 cis region, pre-formatted. All SNPs in LD panel                                                                                                     |
| `sjogren_sumstats`     | 341  | Raw outcome (case-control) | Sjögren’s Disease GWAS, all chromosomes                                                                                                                        |
| `sjogren_outcome`      | 50   | Standardised outcome       | SjD chr20 CD40 region (1 real + 49 synthetic). All SNPs in LD panel                                                                                            |
| `cd40_decode_gwas`     | 10   | deCODE GWAS                | Dummy deCODE format data, input for [`format_pqtl_decode()`](https://github.com/BZuckerman97/mrpipeline/reference/format_pqtl_decode.md)                       |
| `cd40_decode_variants` | 10   | deCODE variants            | Allele frequencies for `cd40_decode_gwas`                                                                                                                      |

### LD reference panel

`inst/extdata/ld_ref.{bed,bim,fam}` — synthetic PLINK binary fileset
with 50 biallelic SNPs from the CD40 chr20 cis region (positions
44646911–44858502) and 100 individuals. Allele frequencies are drawn
from the real `cd40_sumstats` data, so LD structure is random but allele
frequencies are realistic. Total size ~4 KB.

Access the bfile prefix in tests:

``` r
bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
```

Guard integration tests that need plink:

``` r
skip_if_not(file.exists(paste0(bfile, ".bed")), "LD reference panel not available")
```

## PR Checklist

Before opening any pull request:

1.  `air format .` — auto-format R files (CLI tool)
2.  `lintr::lint_package()` — fix any lint warnings
3.  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
    — confirm site builds without errors
4.  `devtools::check()` — must produce 0 errors, 0 warnings
5.  `devtools::test()` — all tests must pass

## Contributing
