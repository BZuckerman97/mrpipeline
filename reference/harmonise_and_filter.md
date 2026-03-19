# Harmonise exposure and outcome data, filter, and deduplicate

Wraps
[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html),
filters to `mr_keep == TRUE`, and removes duplicate SNPs (keeping the
first occurrence).

## Usage

``` r
harmonise_and_filter(exposure, outcome)
```

## Arguments

- exposure:

  Data frame of formatted exposure data.

- outcome:

  Data frame of formatted outcome data.

## Value

A data frame of harmonised data, filtered and deduplicated.
