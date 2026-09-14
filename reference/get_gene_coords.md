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
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Ensembl site unresponsive, trying asia mirror
#> Error in req_perform(html_request): Failed to perform HTTP request.
#> Caused by error in `curl::curl_fetch_memory()`:
#> ! Timeout was reached [jun2026.archive.ensembl.org]:
#> Operation timed out after 9903 milliseconds with 0 bytes received
get_gene_coords(c("CD40", "APOE"), build = "grch37")
#> Error in req_perform(request): HTTP 504 Gateway Timeout.
```
