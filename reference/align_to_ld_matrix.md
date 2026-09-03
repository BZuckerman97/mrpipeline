# Align harmonised data to an LD matrix, correcting allele orientation

Subsets the harmonised data frame and LD matrix to their shared SNPs,
reorders both to match, and determines – per SNP – whether
`beta.exposure`/`beta.outcome` must be sign-flipped before being used
alongside `ld_matrix$ld` in
[`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html)/[`susieR::susie_rss()`](https://rdrr.io/pkg/susieR/man/susie_rss.html).

## Usage

``` r
align_to_ld_matrix(harmonised_data, ld_matrix, verbose = TRUE)
```

## Arguments

- harmonised_data:

  Data frame with `SNP`, `effect_allele.exposure`,
  `other_allele.exposure`, and (if present) `eaf.exposure`/
  `eaf.outcome` columns – i.e. the output of
  [`harmonise_and_filter()`](https://github.com/BZuckerman97/mrpipeline/reference/harmonise_and_filter.md).

- ld_matrix:

  A named list as returned by
  [`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md):
  `ld` and `alleles`.

- verbose:

  Logical. If `TRUE`, emit a
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
  summarising how many SNPs matched/flipped/dropped. Default `TRUE`.

## Value

A named list with elements:

- `data`: subset and reordered data frame (unmatched SNPs dropped;
  beta/allele columns untouched)

- `ld_matrix`: subset and reordered matrix (never re-signed itself)

- `ld_sign`: numeric vector (`+1`/`-1`), same length and order as
  `data`'s rows, to multiply into beta vectors when building coloc/
  SuSiE dataset objects

## Details

[`compute_ld_matrix()`](https://github.com/BZuckerman97/mrpipeline/reference/compute_ld_matrix.md)'s
`ld` matrix is signed relative to the reference panel's own, arbitrary
A1 allele (`ld_matrix$alleles$ld_a1`), which need not match
`harmonised_data$effect_allele.exposure`. `coloc.abf()` is unaffected by
this (it never uses `LD`, and each SNP's Bayes factor depends only on
`beta^2`, so a per-SNP sign flip changes nothing), but `runsusie()` fits
a joint model across all SNPs from `LD` and `beta`/`varbeta` together –
an inconsistent sign for even a subset of SNPs produces an internally
incoherent fit, which is what `susie_rss()`'s `check_prior` safety check
(the "estimated prior variance is unreasonably large" error) exists to
catch.

Rather than mutating `harmonised_data$beta.exposure`/`beta.outcome`
directly (those columns, and the accompanying allele-label columns, are
also used for reporting/instrument export elsewhere and should keep
reflecting the true exposure-outcome harmonisation), this returns a
`ld_sign` vector of `+1`/`-1` to be multiplied into the beta vectors
*only* when constructing the `dataset_exp`/`dataset_out` lists passed to
[`coloc::runsusie()`](https://rdrr.io/pkg/coloc/man/runsusie.html) (see
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)).

SNPs whose alleles don't match the reference panel's A1/A2 at all (e.g.
indels, multi-allelic mismatches) are dropped from both the returned
data and LD matrix. Palindromic SNPs (A/T, C/G) with an EAF in the
ambiguous zone (0.42-0.58, in either exposure or outcome) are also
dropped, since allele-identity matching alone cannot distinguish "same
strand" from "opposite strand, coincidentally same two letters" for
these – this mirrors the equivalent safety filter in the reference
single-cell coloc pipeline (`harmonise_for_coloc()` in
`single_cell_MR_IMID/SSZ_scMR_scripts/Scripts/scMR_onek1k_1M_coloc.R`)
that motivated this fix.
