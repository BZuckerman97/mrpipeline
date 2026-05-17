# mrpipeline User Guide

## Overview

## Gene Coordinate Lookup

Use
[`get_gene_coords()`](https://github.com/BZuckerman97/mrpipeline/reference/get_gene_coords.md)
to programmatically retrieve gene coordinates from Ensembl via biomaRt.
This is useful for defining cis regions in
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md)
without hard-coding coordinates.

``` r

# GRCh38 (default)
coords <- get_gene_coords(c("CD40", "APOE"))
coords
#> # A tibble: 2 × 4
#>   hgnc_symbol chromosome     start       end
#>   <chr>       <chr>          <int>     <int>
#> 1 APOE        19          44905791  44909393
#> 2 CD40        20          44746911  44758502

# GRCh37
coords_37 <- get_gene_coords("CD40", build = "grch37")
```

The returned tibble can feed directly into
[`run_mr()`](https://github.com/BZuckerman97/mrpipeline/reference/run_mr.md)
and
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md):

``` r

cd40 <- get_gene_coords("CD40", build = "grch38")

mr_res <- run_mr(
  exposure = exposure_data,
  exposure_id = "CD40",
  outcome = outcome_data,
  outcome_id = "SjD",
  instrument_region = list(
    chromosome = cd40$chromosome,
    start = cd40$start,
    end = cd40$end
  )
)
```

## Controlling PLINK Resource Usage

When running many parallel jobs (e.g. on an HPC cluster), PLINK’s
default behaviour of auto-detecting available threads and memory can
cause problems — each worker may try to reserve half the node’s RAM. Use
`plink_threads` and `plink_memory` to cap resource usage per call.

You can set these per-call, via R options, or via environment variables:

``` r

# Per-call (highest priority)
result <- run_mr(
 exposure = exposure_data, exposure_id = "CD40",
 outcome = outcome_data, outcome_id = "SjD",
 bfile = "/path/to/ld_ref",
 plink_threads = 1,
 plink_memory = 2000
)

# Via R options (apply to all subsequent calls in the session)
options(
 mrpipeline.plink_threads = 1,
 mrpipeline.plink_memory = 2000
)

# Via environment variables (e.g. in .Renviron or a job script)
# MRPIPELINE_PLINK_THREADS=1
# MRPIPELINE_PLINK_MEMORY=2000
```

Priority order: explicit argument \> R option \> environment variable \>
PLINK auto-detect (`NULL`). The same parameters are available in
[`run_coloc()`](https://github.com/BZuckerman97/mrpipeline/reference/run_coloc.md).

## Setting Up an LD Reference Panel

## Formatting Exposure Data

## Running MR Analyses

### Cis-MR (quick start with API clumping)

### Cis-MR with local LD reference (recommended)

### Sensitivity analyses

### LD-corrected IVW

### Genome-wide MR

### Excluding genomic regions

Use the `exclude_regions` argument to remove instruments falling in
specific genomic regions. This is commonly used to exclude the MHC
region on chromosome 6, which can introduce spurious associations due to
complex LD structure.

Supply a data frame with columns `chr`, `start`, and `end`:

``` r

# Exclude MHC region (GRCh37 coordinates: chr6:28,477,797–33,448,354)
mhc_grch37 <- data.frame(chr = "6", start = 28477797, end = 33448354)

result <- run_mr(
  exposure = exposure_data,
  exposure_id = "PCSK9",
  outcome = outcome_data,
  outcome_id = "CHD",
  exclude_regions = mhc_grch37
)

# GRCh38 coordinates: chr6:28,510,120–33,480,577
mhc_grch38 <- data.frame(chr = "6", start = 28510120, end = 33480577)
```

You can exclude multiple regions by stacking rows. For example, when
studying age-related macular degeneration (AMD), you might exclude both
the CFH and ARMS2/HTRA1 loci:

``` r

amd_exclusions <- data.frame(
  chr = c("1", "10"),
  start = c(196621008, 124214077),
  end = c(196716634, 124274424)
)

result <- run_mr(
  exposure = exposure_data,
  exposure_id = "CFH",
  outcome = outcome_data,
  outcome_id = "AMD",
  exclude_regions = amd_exclusions
)
```

### Manual instrument sets

You can bypass automatic instrument selection by supplying a character
vector of rsIDs to the `instruments` argument:

``` r

result <- run_mr(
  exposure = exposure_data,
  exposure_id = "PCSK9",
  outcome = outcome_data,
  outcome_id = "CHD",
  instruments = c("rs11591147", "rs2479409", "rs11583680")
)
```

By default, instruments missing from the exposure data produce a warning
but the analysis continues with the available SNPs. Set
`instruments_strict = TRUE` to error instead:

``` r

result <- run_mr(
  exposure = exposure_data,
  exposure_id = "PCSK9",
  outcome = outcome_data,
  outcome_id = "CHD",
  instruments = c("rs11591147", "rs2479409"),
  instruments_strict = TRUE
)
```

## Interpreting MR Results

Use [`print()`](https://rdrr.io/r/base/print.html) for a one-line
summary and [`summary()`](https://rdrr.io/r/base/summary.html) for full
details:

``` r

mr_res
summary(mr_res)
```

### Plotting MR results

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) produces
diagnostic plots using TwoSampleMR plotting functions (requires
`ggplot2`):

``` r

plot(mr_res, type = "scatter") # scatter plot (default)
plot(mr_res, type = "forest")  # forest plot
plot(mr_res, type = "funnel")  # funnel plot
```

## Colocalization

### Quick colocalization (ABF only)

The simplest colocalization test uses Approximate Bayes Factors (ABF).
Supply pre-formatted exposure data, raw outcome data, the gene region,
and an LD reference panel:

``` r

result <- run_coloc(
  exposure = formatted_exposure,
  outcome = outcome_data,
  gene_chr = 20,
  gene_start = 44746911,
  gene_end = 44758502,
  exposure_n = 35559,
  outcome_n = 400000,
  bfile = "/path/to/ld_reference",
  methods = "abf"
)
print(result)
summary(result)

# Access posterior probabilities directly
result$coloc_abf$summary
```

### Full colocalization (SuSiE + signals)

For multi-signal colocalization, include `"susie"` and `"signals"` in
the `methods` argument. SuSiE fine-maps each trait independently, then
`coloc.susie()` and `coloc.signals()` test colocalization across all
pairs of credible sets:

``` r

result <- run_coloc(
  exposure = formatted_exposure,
  outcome = outcome_data,
  gene_chr = 20,
  gene_start = 44746911,
  gene_end = 44758502,
  exposure_n = 35559,
  outcome_n = 400000,
  bfile = "/path/to/ld_reference",
  methods = c("abf", "susie", "signals")
)
summary(result)
```

### Case-control outcomes

For case-control outcomes, set `outcome_type = "cc"` and provide the
proportion of cases via `outcome_s`:

``` r

result <- run_coloc(
  exposure = formatted_exposure,
  outcome = outcome_data,
  gene_chr = 20,
  gene_start = 44746911,
  gene_end = 44758502,
  exposure_n = 35559,
  outcome_n = 400000,
  outcome_type = "cc",
  outcome_s = 0.3,
  bfile = "/path/to/ld_reference",
  methods = c("abf", "susie", "signals")
)
```

### Plotting coloc results

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
`coloc_result` objects supports two plot types (requires `ggplot2`):

``` r

plot(coloc_res, type = "pp_bar")    # bar chart of ABF posterior probabilities (default)
plot(coloc_res, type = "regional")  # regional association plots
```

## Two-Stage Batch Workflow
