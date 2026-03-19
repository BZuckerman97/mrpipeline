# Plot a coloc_result object

Creates visualisations for colocalization results. Requires the
`ggplot2` package.

## Usage

``` r
# S3 method for class 'coloc_result'
plot(x, type = c("pp_bar", "regional"), ...)
```

## Arguments

- x:

  A `coloc_result` object.

- type:

  Character. Plot type: `"pp_bar"` (default) for a bar chart of ABF
  posterior probabilities (H0–H4), or `"regional"` for side-by-side
  regional association plots.

- ...:

  Ignored.

## Value

A ggplot object. Returns `NULL` invisibly if the result status is not
`"success"` or if the required data is unavailable.

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
} # }
```
