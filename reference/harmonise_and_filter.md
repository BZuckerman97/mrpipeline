# Harmonise exposure and outcome data, filter, and deduplicate

Wraps
[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html),
runs the allele orientation check
([`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md))
on the raw harmonised output, filters to `mr_keep == TRUE`, and removes
duplicate SNPs (keeping the first occurrence).

## Usage

``` r
harmonise_and_filter(
  exposure,
  outcome,
  allele_check = c("error", "warn", "none"),
  action = 2,
  check = TRUE,
  verbose = FALSE
)
```

## Arguments

- exposure:

  Data frame of formatted exposure data.

- outcome:

  Data frame of formatted outcome data.

- allele_check:

  One of `"error"` (default), `"warn"` or `"none"`. Passed to
  [`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md).

- action:

  `1`, `2` (default) or `3`. Passed to
  [`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html);
  see
  [`validate_harmonise_action()`](https://github.com/BZuckerman97/mrpipeline/reference/validate_harmonise_action.md)
  for what each level does.

- check:

  Logical. If `FALSE`, skip the allele orientation check entirely
  (nothing is recorded). Used by
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
  which has already checked its instruments together with a sample of
  shared SNPs via
  [`check_allele_orientation_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation_gwas.md).
  Default `TRUE`.

- verbose:

  Logical. Passed to
  [`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md).
  Default `FALSE`. A check that cannot reach a verdict warns regardless.

## Value

A data frame of harmonised data, filtered and deduplicated.

## Details

The allele orientation check runs *before* the `mr_keep` filter so that
SNPs dropped later (e.g. for a missing `beta`/`se`) still contribute
their allele frequencies to the verdict, and it runs on every call –
including a zero-row harmonisation – so the record returned by
[`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)
always describes the most recent exposure/outcome pair rather than a
stale one.
