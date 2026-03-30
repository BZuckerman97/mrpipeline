# Resolve PLINK resource option from R option or environment variable

Checks `getOption("mrpipeline.plink_{param}")` first, then falls back to
the environment variable `MRPIPELINE_PLINK_{PARAM}`. Returns `NULL` if
neither is set (PLINK auto-detects).

## Usage

``` r
plink_option(param)
```

## Arguments

- param:

  Either `"threads"` or `"memory"`.

## Value

Integer or `NULL`.
