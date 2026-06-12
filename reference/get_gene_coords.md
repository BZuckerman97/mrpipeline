# Look up gene coordinates from Ensembl via biomaRt

Queries Ensembl for genomic coordinates of one or more HGNC gene
symbols. Results are filtered to standard chromosomes (1–22, X, Y),
deduplicated (widest range per gene, preferring autosomes), and returned
as a tibble.

## Usage

``` r
get_gene_coords(genes, build = c("grch38", "grch37"))
```

## Arguments

- genes:

  Character vector of HGNC gene symbols.

- build:

  Genome build: `"grch38"` (default) or `"grch37"`.

## Value

A tibble with columns `hgnc_symbol`, `chromosome`, `start`, `end`. Genes
not found in Ensembl are dropped with a warning.

## Examples

``` r
get_gene_coords("CD40")
#> Ensembl site unresponsive, trying asia mirror
#> # A tibble: 1 × 4
#>   hgnc_symbol chromosome    start      end
#>   <chr>       <chr>         <int>    <int>
#> 1 CD40        20         46118240 46130378
get_gene_coords(c("CD40", "APOE"), build = "grch37")
#> Error in req_perform(request): HTTP 504 Gateway Timeout.
```
