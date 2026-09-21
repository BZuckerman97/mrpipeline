# Allele orientation check on instruments plus a sample of shared SNPs

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)'s
instrument set is often too small for
[`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md)
to reach a verdict (a cis-MR may have three instruments; the check needs
ten informative non-palindromic SNPs). But
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
holds the *full* exposure and outcome GWAS before it narrows the outcome
to instrument rsIDs, so this helper harmonises the instruments together
with up to `n_sample` additional SNPs shared by the two datasets and
runs the verdict on that larger set. The extra SNPs are taken at evenly
spaced positions through the sorted shared rsIDs – deterministic, so
results are reproducible and the caller's RNG state is untouched, and
effectively random with respect to genomic position.

## Usage

``` r
check_allele_orientation_gwas(
  exposure,
  outcome,
  instrument_snps,
  allele_check = c("error", "warn", "none"),
  n_sample = 1000L,
  action = 2,
  verbose = FALSE,
  call = rlang::caller_env()
)
```

## Arguments

- exposure:

  Data frame of TwoSampleMR-formatted exposure data (the full dataset
  passed to
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
  not just the instruments).

- outcome:

  Data frame in
  [`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
  outcome format (`rsids`, `effect_allele`, `other_allele`, `beta`,
  `se`, `eaf`, ...).

- instrument_snps:

  Character vector of instrument rsIDs to always include in the check
  set.

- allele_check:

  One of `"error"` (default), `"warn"` or `"none"`.

- n_sample:

  Integer. Maximum number of non-instrument shared SNPs to add. Default
  `1000L`.

- action:

  `1`, `2` (default) or `3`. Passed to
  [`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html),
  so that the check describes the same harmonisation the analysis will
  use. The verdict itself is unaffected: `action` gates only the
  frequency-based second flip applied to palindromic variants, which the
  check excludes anyway.

- verbose:

  Logical. Passed to
  [`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md).
  A check that cannot reach a verdict warns regardless.

- call:

  Environment reported in the condition. Default
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html).

## Value

The diagnostic record, invisibly (see
[`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)).

## Details

Cost scales with the total rows of the two inputs, not with the number
of instruments: the rsID intersection and the radix sort of the shared
rsIDs are the only steps that touch every row (about 6 s for a 42M-row
exposure against a 12M-row outcome, i.e. two genome-wide GWAS; well
under a second for a 200k-SNP exposure). Formatting and harmonising the
~1000 sampled SNPs takes a few milliseconds. The shared rsIDs are sorted
with `method = "radix"`: base
[`sort()`](https://rdrr.io/r/base/sort.html) uses shell sort with
per-comparison locale collation for character vectors, which is ~25x
slower on 10M rsIDs and made the check the dominant cost of a cis-MR
(issue \#40). The check needs the *full* frames to sample from, so it
cannot be made cheaper by narrowing the inputs first, and it runs on
every call – it is not cached across
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
calls that share a pair.

The check is skipped (with a `"skipped"` record, and a warning that
orientation is unverified) when the exposure lacks `SNP`/`eaf.exposure`
or the outcome lacks the
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
outcome columns, or when nothing overlaps.
