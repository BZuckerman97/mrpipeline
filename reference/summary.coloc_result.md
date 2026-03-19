# Summarise a coloc_result object

Displays full details: region info, all posterior probabilities (H0-H4),
SuSiE credible set pairs, signals hits, and skipped methods.

## Usage

``` r
# S3 method for class 'coloc_result'
summary(object, ...)
```

## Arguments

- object:

  A `coloc_result` object.

- ...:

  Ignored.

## Value

Invisibly returns `object`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
summary(result)
} # }
```
