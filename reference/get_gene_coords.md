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
#> Ensembl site unresponsive, trying useast mirror
#> Error in .chooseEnsemblMirror(mirror = mirror, http_config = http_config): Unable to query any Ensembl site
get_gene_coords(c("CD40", "APOE"), build = "grch37")
#> Error in req_perform(request): Failed to perform HTTP request.
#> Caused by error in `curl::curl_fetch_memory()`:
#> ! Unsupported protocol [grch37.ensembl.org]:
#> Received HTTP/0.9 when not allowed
```
