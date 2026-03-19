# Print a coloc_result object

Displays a one-line summary: number of SNPs, ABF PP.H4 if available,
SuSiE max PP.H4 and credible set count if available. Shows status if not
success.

## Usage

``` r
# S3 method for class 'coloc_result'
print(x, ...)
```

## Arguments

- x:

  A `coloc_result` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_coloc(
  exposure = cd40_exposure,
  outcome = sjogren_outcome,
  gene_chr = "20", gene_start = 44746911, gene_end = 44758502,
  bfile = system.file("extdata", "ld_ref", package = "mrpipeline")
)
print(result)
} # }
```
