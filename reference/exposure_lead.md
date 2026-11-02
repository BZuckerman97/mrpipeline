# Lead exposure SNP

The SNP with the lowest `pval.exposure`. Ties go to the largest
\|beta/se\|: GWAS p-values are often floored (e.g. at 1e-200), which can
tie many SNPs at a strong locus, and row order says nothing about which
of them leads.

## Usage

``` r
exposure_lead(dat)
```

## Arguments

- dat:

  Harmonised data frame.

## Value

The rsID, or `NULL` if every `pval.exposure` is missing.
