# mrpipeline — Package Developer Guide for Claude Code

## Package Overview

`mrpipeline` is an R package for Mendelian randomisation (MR) and colocalization analysis, focused on proteomic GWAS data.

**Core exported functions:**

- `run_mr()` — Cis-MR, genome-wide MR, or manual-instrument MR with sensitivity methods
- `run_coloc()` — Colocalization (coloc.abf, SuSiE, coloc.signals, colocPropTest)
- `format_pqtl_decode()` — Format deCODE proteomics GWAS to TwoSampleMR exposure format
- `format_pqtl_ukbppp()` — Format UKB-PPP pQTL data to TwoSampleMR exposure format
- `format_single_cell_onek1k()` — Format OneK1K single-cell eQTL data

**Internal helpers** (not exported, in `R/helpers.R`):
`harmonise_and_filter()`, `compute_ld_matrix()`, `clump_instruments()`, `align_to_ld_matrix()`, `eaf_to_maf()`, `resolve_sample_size()`

**S3 classes:** `mr_result` (from `run_mr()`), `coloc_result` (from `run_coloc()`)

## Build and Check Commands

```r
devtools::load_all()      # load package for interactive development
devtools::document()      # regenerate NAMESPACE and man/ from roxygen2
devtools::check()         # R CMD Check (must pass: 0 errors, 0 warnings)
devtools::test()          # run testthat tests
pkgdown::build_site()     # build documentation site (must succeed)
```

```bash
air format .              # auto-format R code (run before every commit)
```

```r
lintr::lint_package()     # style/lint checks (run before every commit)
```

## Code Conventions

**Messages / warnings / errors:** Use `cli` package exclusively.
```r
cli::cli_inform("Loading {protein} data...")   # informational
cli::cli_warn("Only {n} SNPs retained.")       # warning
cli::cli_abort("bfile is required for LD-corrected MR.")  # error
```

**Pipe:** Native `|>` only. Never `%>%`.

**String operations:** Use `stringr` functions.
```r
stringr::str_remove(x, "_.*")      # not gsub("_.*", "", x)
stringr::str_detect(x, "^chr")     # not grepl("^chr", x)
stringr::str_starts(x, "rs")       # not startsWith(x, "rs")
```

**Namespace:** Prefer `pkg::fn()` over `@importFrom`.

**Internal helpers:** Tag with `@keywords internal`; do NOT use `@export`.

**Documentation:** Use roxygen2 for all functions. Vignettes (`mrpipeline-user-guide.Rmd`, `mrpipeline-developer-guide.Rmd`) are **living documents** — update them in the same commit as any API change:
- When adding/changing/removing a function parameter: update roxygen, user guide (usage examples), developer guide (architecture/internals)
- When adding/removing internal helpers: update the developer guide's helpers section
- When modifying S3 class structure: update the developer guide's S3 classes section

**`verbose` argument:** Planned future feature. When added, gate all `cli::cli_inform()` calls behind it.

## PR Checklist

Before opening any pull request:

1. `air format .` — auto-format R files
2. `lintr::lint_package()` — fix any lint warnings
3. `pkgdown::build_site()` — confirm site builds without errors
4. `devtools::check()` — must produce 0 errors, 0 warnings
5. `devtools::test()` — all tests must pass

## Workflow

Feature branches are created off `dev` and merged back to `dev` via PR. Each feature should have a corresponding GitHub issue.

Branch naming: `phase-N/short-description` (e.g. `phase-1/shared-helpers`).

## Architecture Pointers

- Full architectural context: `vignettes/mrpipeline-developer-guide.Rmd`
- End-to-end usage examples: `vignettes/mrpipeline-user-guide.Rmd`
- Test data: `R/data.R` (lazy-loaded datasets), `inst/extdata/` (LD reference panel)
- Tests requiring plink/LD reference: use `testthat::skip_if_not(file.exists(bfile))` for CI safety
