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
  posterior probabilities (H0–H4), `"regional"` for side-by-side
  regional association plots, or `"locuszoom"` for LD-coloured regional
  plots via the `locuszoomr` package (see
  [`plot_coloc_locuszoom()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_locuszoom.md)
  for its additional required arguments, passed through `...`).

- ...:

  For `type = "locuszoom"`, forwarded to
  [`plot_coloc_locuszoom()`](https://github.com/BZuckerman97/mrpipeline/reference/plot_coloc_locuszoom.md)
  (`ens_db`, `bfile`, and optionally `plink_bin`/ `index_snp`). Ignored
  for other types.

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
result <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
plot(result, type = "pp_bar")
plot(result, type = "regional")

# locuszoomr resolves a character ens_db via the search path, so the
# annotation package must be attached, not just installed.
library(EnsDb.Hsapiens.v75)
plot(
  result,
  type = "locuszoom",
  ens_db = "EnsDb.Hsapiens.v75",
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
} # }
```
