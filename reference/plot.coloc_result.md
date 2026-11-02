# Plot a coloc_result object

Creates visualisations for colocalization results. Requires the
`ggplot2` package.

## Usage

``` r
# S3 method for class 'coloc_result'
plot(x, type = c("pp_bar", "regional", "locuszoom"), ...)
```

## Arguments

- x:

  A `coloc_result` object.

- type:

  Character. Plot type: `"pp_bar"` (default) for a bar chart of ABF
  posterior probabilities (H0–H4), or `"regional"` for exposure and
  outcome regional association plots, coloured by LD when r2 is
  available and with an optional gene track (see
  [`plot_coloc_regional()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_regional.md)).
  `"locuszoom"` is deprecated in favour of `"regional"` with `ens_db`.

- ...:

  For `type = "regional"`, forwarded to
  [`plot_coloc_regional()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_regional.md)
  (`ens_db`, `bfile`, `plink_bin` and `index_snp`). For
  `type = "locuszoom"`, forwarded to
  [`plot_coloc_locuszoom()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_locuszoom.md).
  Ignored for `"pp_bar"`.

## Value

For `type = "pp_bar"`/`"regional"`, a ggplot object (or `NULL` invisibly
if the result status is not `"success"` or the required data is
unavailable). For `type = "locuszoom"`, `NULL` invisibly always – see
[`plot_coloc_locuszoom()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_locuszoom.md),
which draws directly to the current graphics device rather than
returning a re-renderable object.

## Examples

``` r
if (FALSE) { # \dontrun{
bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
result <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
  bfile = bfile
)
plot(result, type = "pp_bar")
plot(result, type = "regional")

# Colour by LD with another SNP (recomputed, so it needs the panel)
plot(result, type = "regional", index_snp = "rs1883832", bfile = bfile)

# Add a gene track. A character ens_db must be attached, not just installed.
library(EnsDb.Hsapiens.v75)
plot(result, type = "regional", ens_db = "EnsDb.Hsapiens.v75")
} # }
```
