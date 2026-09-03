# Skip VCF meta-information lines and return data as a data.frame

Counts lines beginning with `##`, reads the `#CHROM` header, then uses
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
to read the remaining rows with the correct column names.

## Usage

``` r
read_vcf_data(file)
```

## Arguments

- file:

  Path to a VCF file (plain or gzip-compressed).

## Value

A data frame of VCF records.
