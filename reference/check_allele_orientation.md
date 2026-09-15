# Detect swapped effect/other alleles from harmonised allele frequencies

Some GWAS files label their allele columns `A1`/`A2` meaning REF/ALT,
with `BETA` and the frequency column oriented to `A2` (EPACTS/RAREMETAL
style), whereas
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
– like PLINK, regenie and METAL – reads `A1` as the effect allele.
Feeding such a file in swaps effect and other allele for every variant
and silently inverts every beta (GitHub issue \#18). This check detects
that at harmonisation time, the only point where two datasets coexist,
so no external reference panel is needed.

## Usage

``` r
check_allele_orientation(
  harmonised,
  allele_check = c("error", "warn", "none"),
  threshold = 0.7,
  min_n = 10L,
  n_sampled = NA_integer_,
  verbose = FALSE,
  call = rlang::caller_env()
)
```

## Arguments

- harmonised:

  Raw output of
  [`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
  (before any `mr_keep` filtering). Needs columns `SNP`, `eaf.exposure`,
  `eaf.outcome`, `effect_allele.exposure`, `other_allele.exposure`,
  `palindromic` and `remove`; otherwise the check is skipped.

- allele_check:

  One of `"error"` (default), `"warn"` or `"none"`. Controls what
  happens on failure; the diagnostic record is stored in every mode.

- threshold:

  Proportion of informative variants that must be complementary for the
  check to fail. Default `0.70`.

- min_n:

  Minimum number of informative variants required to reach a verdict.
  Default `10L`.

- n_sampled:

  Integer. Number of non-instrument SNPs that
  [`check_allele_orientation_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation_gwas.md)
  added to the harmonised set, for the record only. `NA` (default) when
  the check ran on a harmonisation that was not sampled.

- verbose:

  Logical. If `TRUE`, report a passing or skipped verdict via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html).
  Default `FALSE`.

- call:

  Environment. The calling frame reported in the condition. Default
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html).

## Value

The diagnostic record (see
[`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)
for its structure), invisibly. The same record is stored so that
[`last_allele_check()`](https://github.com/BZuckerman97/mrpipeline/reference/last_allele_check.md)
can return it – this matters on the error path, where the return value
is otherwise lost.

## How the verdict is reached

[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
aligns non-palindromic variants purely by allele letters, negating
`beta.outcome` and replacing `eaf.outcome` with `1 - eaf.outcome`
whenever the outcome's effect allele is the exposure's other allele.
When one dataset's allele labels are swapped, that alignment is applied
to *every* variant, so after harmonisation `eaf.outcome` ends up
describing the other allele: `eaf.outcome ~ 1 - eaf.exposure` across the
set. Genuine cohort differences (e.g. ancestry) produce *scatter* around
`eaf.exposure`; this bug produces systematic *complementarity*.

The statistic is therefore: among *informative* variants – both EAFs
present, `palindromic == FALSE`, `remove == FALSE`, one row per SNP –
the proportion whose `eaf.exposure` is closer to `1 - eaf.outcome` than
to `eaf.outcome` (ties, e.g. `eaf.outcome == 0.5`, count as *not*
complementary). The check fails when that proportion exceeds `threshold`
and at least `min_n` informative variants were available; with fewer it
is `"skipped"`. Variants with EAF near 0.5 are equally likely to fall
either side, so they can only dilute the proportion towards 0.5 – they
cannot cause a spurious failure, only mask a real one, which the 0.70
threshold tolerates.

## Why palindromic variants are excluded

For A/T and C/G variants the strand cannot be resolved from letters, so
`harmonise_data()` resolves it *from the allele frequencies*: after the
letter-based swap it flips again if `eaf.exposure` and `eaf.outcome` sit
on opposite sides of 0.5. When the effect allele and the frequency are
both mis-assigned those two flips cancel, so a palindromic variant's
`eaf.outcome` always looks consistent and its harmonised beta can be
*identical* to the correct value even though every non-palindromic beta
in the same set is inverted. Including them would only dilute the
statistic; and once the check fails, every palindromic strand call in
that pair is unreliable regardless of what its beta looks like.

## What the check cannot tell you

The comparison is symmetric: a failure means the two datasets disagree,
not which one is wrong. Break the tie with an independent frequency
reference (`plink --freq` on the LD panel; see `ref_frq` in
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)),
a variant with a well-established effect direction, or provenance (a
dataset that has harmonised cleanly against others is not the suspect).
