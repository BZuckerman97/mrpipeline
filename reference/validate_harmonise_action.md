# Validate a TwoSampleMR harmonisation action level

[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
accepts `action` as a vector, applying a different level per outcome.
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
handle exactly one outcome, so a vector here is a mistake worth catching
rather than silently recycling.

## Usage

``` r
validate_harmonise_action(action)
```

## Arguments

- action:

  Value to validate.

## Value

`action`, unchanged. It is deliberately not coerced to integer:
[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html)
embeds `action` as a column in its output, so coercing `2` to `2L` would
make `mrpipeline`'s harmonised frame differ from a plain
`harmonise_data()` call on the same data by the storage mode of that
column alone.

## Details

The three levels:

|  |  |
|----|----|
| `action` | Behaviour |
| 1 | Assume all alleles are on the forward strand: no frequency-based flip |
| 2 | Infer the positive strand, resolving palindromes from allele frequencies (default) |
| 3 | As 2, but set `mr_keep = FALSE` for every palindromic, ambiguous or incompatible SNP |

Only the palindrome handling differs. The letter-based alignment of
non-palindromic variants – negating `beta.outcome` and replacing
`eaf.outcome` with `1 - eaf.outcome` when the outcome's effect allele is
the exposure's other allele – happens at every level, which is why
[`check_allele_orientation()`](https://github.com/BZuckerman97/mrpipeline/reference/check_allele_orientation.md)'s
verdict does not depend on `action`.
