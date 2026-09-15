# Retrieve the most recent allele orientation check

Every call to
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
or
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
runs an allele orientation check after harmonisation (see
[`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md))
and stores its full diagnostic record, whether it passed, failed or was
skipped. This function returns that record, following the
[`httr2::last_response()`](https://httr2.r-lib.org/reference/last_response.html)
pattern. It is most useful after an `mrpipeline_allele_check_error` –
when the analysis function's return value is lost – to inspect every
variant behind the verdict, but the record is stored on a pass too, so a
borderline result can still be examined.

## Usage

``` r
last_allele_check()
```

## Value

The record list, or `NULL` (with a message) if no check has run yet in
this session.

## Details

The check detects a GWAS file whose `A1`/`A2` columns mean REF/ALT
rather than effect/other, which silently inverts every beta. See the
*What does A1 mean?* section of
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
for the two conventions and how to tell them apart.

## Record structure

A named list:

- `status`: `"pass"`, `"fail"` or `"skipped"` (fewer than `min_n`
  informative SNPs, or no allele frequencies).

- `n`: number of informative SNPs – non-palindromic, alleles compatible,
  both `eaf.exposure` and `eaf.outcome` present.

- `n_complementary`: how many of those have `eaf.exposure` closer to
  `1 - eaf.outcome` than to `eaf.outcome`.

- `prop`: `n_complementary / n`.

- `threshold`, `min_n`: the criteria in force (fail when
  `prop > threshold` and `n >= min_n`).

- `allele_check`: the mode in force (`"error"`, `"warn"` or `"none"`).

- `exposure`, `outcome`: phenotype names of the harmonised pair (the
  `exposure_id`/`outcome_id` passed to
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
  or
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)),
  plus `id.exposure`, `id.outcome`, TwoSampleMR's internal identifiers.

- `n_sampled`: for
  [`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md),
  the number of non-instrument SNPs shared by the two datasets that were
  added to the check set (up to 1000); `NA` for
  [`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md),
  which checks every SNP in its window.

- `variants`: a data frame with one row per SNP that had both allele
  frequencies – palindromic and allele-incompatible SNPs included but
  flagged – with columns `SNP`, `effect_allele`, `other_allele`,
  `eaf.exposure`, `eaf.outcome`, `eaf.outcome_flipped`
  (`1 - eaf.outcome`), `palindromic`, `remove`, `informative`, `score`
  (positive when complementary; larger is more clear-cut) and
  `complementary` (`NA` for non-informative rows). Sort by `score`
  descending to see the worst offenders.

## See also

[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
(`allele_check` argument),
[`format_gwas()`](https://github.com/BZuckerman97/mrpipeline/reference/format_gwas.md)
(section *What does A1 mean?*).

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mr(exposure, "CD40", outcome, "AMD", instruments = ivs)

chk <- last_allele_check()
chk$status
chk$variants[order(-chk$variants$score), ]
} # }
```
