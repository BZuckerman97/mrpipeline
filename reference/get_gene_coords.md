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
#> Error in .processResults(postRes, mart = mart, hostURLsep = sep, fullXmlQuery = fullXmlQuery,     quote = quote, numAttributes = length(attributes)): Query ERROR: caught BioMart::Exception::Database: Could not connect to mysql database ensembl_mart_115: DBI connect('database=ensembl_mart_115;host=127.0.0.1;port=5316','ensro',...) failed: Can't connect to MySQL server on '127.0.0.1' (111) at /nfs/public/ro/ensweb/live/mart/www_115/biomart-perl/lib/BioMart/Configuration/DBLocation.pm line 98.
get_gene_coords(c("CD40", "APOE"), build = "grch37")
#> Error in req_perform(request): HTTP 504 Gateway Timeout.
```
