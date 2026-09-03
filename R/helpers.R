# Internal helper functions shared by run_mr() and run_coloc()
# All functions in this file are @keywords internal and NOT exported.

#' Resolve PLINK resource option from R option or environment variable
#'
#' Checks `getOption("mrpipeline.plink_{param}")` first, then falls back to
#' the environment variable `MRPIPELINE_PLINK_{PARAM}`. Returns `NULL` if
#' neither is set (PLINK auto-detects).
#'
#' @param param Either `"threads"` or `"memory"`.
#' @return Integer or `NULL`.
#' @keywords internal
plink_option <- function(param) {
  opt <- getOption(paste0("mrpipeline.plink_", param))
  if (!is.null(opt)) {
    return(as.integer(opt))
  }
  env <- Sys.getenv(paste0("MRPIPELINE_PLINK_", toupper(param)), unset = "")
  if (nzchar(env)) as.integer(env) else NULL
}

#' Harmonise exposure and outcome data, filter, and deduplicate
#'
#' Wraps [TwoSampleMR::harmonise_data()], filters to `mr_keep == TRUE`,
#' and removes duplicate SNPs (keeping the first occurrence).
#'
#' @param exposure Data frame of formatted exposure data.
#' @param outcome Data frame of formatted outcome data.
#'
#' @return A data frame of harmonised data, filtered and deduplicated.
#'
#' @importFrom rlang .data
#' @keywords internal
harmonise_and_filter <- function(exposure, outcome) {
  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure,
    outcome_dat = outcome
  )

  # Guard: return 0-row frame if harmonisation produced no usable output
  # (e.g. no SNP overlap, or all SNPs removed as palindromic). This lets
  # run_mr() hit its nrow == 0 early-return rather than erroring on a missing
  # `mr_keep` column.
  if (nrow(harmonised) == 0 || !"mr_keep" %in% names(harmonised)) {
    return(harmonised[0, , drop = FALSE])
  }

  harmonised <- harmonised |>
    dplyr::filter(.data$mr_keep == TRUE) |>
    dplyr::filter(!duplicated(.data$SNP))

  harmonised
}

#' Compute LD correlation matrix from a local reference panel
#'
#' Calls [ieugwasr::ld_matrix()], which (for a local `bfile`) runs
#' `plink --r square --keep-allele-order` and names rows/columns
#' `rsid_A1_A2`, where A1/A2 are the reference panel's own `.bim` allele
#' coding. This is a **signed** correlation matrix (`--r`, not `--r2`) and
#' the sign of every entry is anchored to that specific, arbitrary A1
#' allele -- which need not match the effect allele used anywhere else in
#' the pipeline. That allele identity is only recoverable from the
#' `rsid_A1_A2` name suffix, so it is parsed out here and returned
#' alongside the matrix (rather than discarded) so [align_to_ld_matrix()]
#' can re-orient betas/z-scores to this matrix's sign convention before
#' they are used together in [coloc::runsusie()]/`susieR::susie_rss()`.
#' See `run_coloc()`'s "LD matrix orientation" section for the full
#' rationale.
#'
#' @param snps Character vector of rsIDs.
#' @param bfile Path to PLINK bfile prefix (without .bed/.bim/.fam).
#' @param plink_bin Path to PLINK binary. If `NULL`, auto-detected via
#'   [genetics.binaRies::get_plink_binary()].
#' @param plink_threads Number of threads for PLINK. `NULL` = auto-detect.
#' @param plink_memory Memory (MB) for PLINK. `NULL` = auto-detect.
#'
#' @return A named list:
#'   - `ld`: square signed correlation matrix with rsID-only row/column
#'     names
#'   - `alleles`: data frame with columns `SNP`, `ld_a1`, `ld_a2` giving the
#'     reference panel's allele coding for each SNP in `ld` -- the allele
#'     that `ld`'s sign is anchored to (`ld_a1`) and its counterpart
#'
#' @keywords internal
compute_ld_matrix <- function(
  snps,
  bfile,
  plink_bin = NULL,
  plink_threads = NULL,
  plink_memory = NULL
) {
  if (is.null(plink_bin)) {
    plink_bin <- genetics.binaRies::get_plink_binary()
  }

  ld <- ieugwasr::ld_matrix(
    variants = snps,
    bfile = bfile,
    plink_bin = plink_bin,
    threads = plink_threads,
    memory = plink_memory
  )

  full_names <- rownames(ld)
  rsid <- stringr::str_remove(full_names, "_.*")
  parts <- stringr::str_match(full_names, "^[^_]+_([^_]+)_([^_]+)$")

  alleles <- data.frame(
    SNP = rsid,
    ld_a1 = parts[, 2],
    ld_a2 = parts[, 3],
    stringsAsFactors = FALSE
  )

  rownames(ld) <- rsid
  colnames(ld) <- rsid

  list(ld = ld, alleles = alleles)
}

#' LD (r2) of every SNP to a specified index SNP, via local reference panel
#'
#' Reuses [compute_ld_matrix()] -- no separate LD mechanism, no network
#' calls, and consistent sign/allele handling with the rest of the package.
#' [compute_ld_matrix()]'s matrix is signed (`--r`, not `--r2`) but squaring
#' removes the sign, which is all a colocalization regional plot's LD-colour
#' scale needs.
#'
#' @param snps Character vector of rsIDs.
#' @param index_snp rsID of the SNP to compute LD against. Must be one of
#'   `snps`.
#' @param bfile Path to PLINK bfile prefix (without .bed/.bim/.fam).
#' @param plink_bin Path to PLINK binary. If `NULL`, auto-detected via
#'   [genetics.binaRies::get_plink_binary()].
#' @param plink_threads Number of threads for PLINK. `NULL` = auto-detect.
#' @param plink_memory Memory (MB) for PLINK. `NULL` = auto-detect.
#'
#' @return A data frame with columns `SNP` and `r2`.
#'
#' @keywords internal
compute_ld_to_index <- function(
  snps,
  index_snp,
  bfile,
  plink_bin = NULL,
  plink_threads = NULL,
  plink_memory = NULL
) {
  m <- compute_ld_matrix(snps, bfile, plink_bin, plink_threads, plink_memory)
  r <- m$ld[, index_snp]
  data.frame(SNP = names(r), r2 = as.numeric(r)^2, stringsAsFactors = FALSE)
}

#' Clump instruments using LD reference
#'
#' Computes pseudo p-values (to ensure correct ranking by PLINK) and
#' performs LD clumping via [ieugwasr::ld_clump()].
#'
#' Pseudo p-values are generated by:
#' 1. Computing z-scores: `beta / se`
#' 2. Computing precise -log10(p) via [Rmpfr::pnorm()]
#' 3. Ranking SNPs by -log10(p) descending
#' 4. Assigning evenly spaced pseudo p-values from 1e-100 to 0.9
#'
#' @param dat Data frame with columns `rsid`, `pval`, `id`, and optionally
#'   `beta` and `se` (used for pseudo p-value computation). If `beta` and
#'   `se` are absent, `pval` is used directly.
#' @param rsq_thresh R-squared clumping threshold.
#' @param clump_kb Clumping window in kb. Default 10000.
#' @param bfile Path to PLINK bfile prefix. If `NULL`, uses API clumping.
#' @param plink_bin Path to PLINK binary. If `NULL`, auto-detected.
#' @param pop Population for API clumping. Default `"EUR"`.
#'
#' @return Data frame of clumped variants (output of `ieugwasr::ld_clump()`).
#'
#' @keywords internal
clump_instruments <- function(
  dat,
  rsq_thresh,
  clump_kb = 10000,
  bfile = NULL,
  plink_bin = NULL,
  pop = "EUR",
  plink_threads = NULL,
  plink_memory = NULL
) {
  # Compute pseudo p-values if beta and se are available
  if (all(c("beta", "se") %in% colnames(dat))) {
    dat <- dat |>
      dplyr::mutate(
        z = .data$beta / .data$se,
        log10p = sapply(
          -log10(
            2 *
              Rmpfr::pnorm(
                Rmpfr::mpfr(abs(.data$z), precBits = 100),
                lower.tail = FALSE
              )
          ),
          as.numeric
        )
      ) |>
      dplyr::arrange(dplyr::desc(.data$log10p)) |>
      dplyr::mutate(
        pval = seq(from = 1e-100, to = 0.9, length.out = dplyr::n())
      )
  }

  clump_args <- list(
    dat = dat[, c("rsid", "pval", "id")],
    clump_kb = clump_kb,
    clump_r2 = rsq_thresh
  )

  if (!is.null(bfile)) {
    if (is.null(plink_bin)) {
      plink_bin <- genetics.binaRies::get_plink_binary()
    }
    clump_args$bfile <- path.expand(bfile)
    clump_args$plink_bin <- plink_bin
    clump_args$threads <- plink_threads
    clump_args$memory <- plink_memory
  } else {
    if (!isTRUE(getOption("mrpipeline.api_message_shown"))) {
      cli::cli_inform(
        c(
          "i" = "Using API for LD clumping (population: {pop}).",
          "i" = "For large datasets, a local {.arg bfile} is recommended."
        )
      )
      options(mrpipeline.api_message_shown = TRUE)
    }
    clump_args$pop <- pop
  }

  do.call(ieugwasr::ld_clump, clump_args)
}

#' Align harmonised data to an LD matrix, correcting allele orientation
#'
#' Subsets the harmonised data frame and LD matrix to their shared SNPs,
#' reorders both to match, and determines -- per SNP -- whether
#' `beta.exposure`/`beta.outcome` must be sign-flipped before being used
#' alongside `ld_matrix$ld` in [coloc::runsusie()]/`susieR::susie_rss()`.
#'
#' [compute_ld_matrix()]'s `ld` matrix is signed relative to the reference
#' panel's own, arbitrary A1 allele (`ld_matrix$alleles$ld_a1`), which need
#' not match `harmonised_data$effect_allele.exposure`. `coloc.abf()` is
#' unaffected by this (it never uses `LD`, and each SNP's Bayes factor
#' depends only on `beta^2`, so a per-SNP sign flip changes nothing), but
#' `runsusie()` fits a joint model across all SNPs from `LD` and
#' `beta`/`varbeta` together -- an inconsistent sign for even a subset of
#' SNPs produces an internally incoherent fit, which is what
#' `susie_rss()`'s `check_prior` safety check (the "estimated prior
#' variance is unreasonably large" error) exists to catch.
#'
#' Rather than mutating `harmonised_data$beta.exposure`/`beta.outcome`
#' directly (those columns, and the accompanying allele-label columns, are
#' also used for reporting/instrument export elsewhere and should keep
#' reflecting the true exposure-outcome harmonisation), this returns a
#' `ld_sign` vector of `+1`/`-1` to be multiplied into the beta vectors
#' *only* when constructing the `dataset_exp`/`dataset_out` lists passed to
#' `coloc::runsusie()` (see `run_coloc()`).
#'
#' SNPs whose alleles don't match the reference panel's A1/A2 at all (e.g.
#' indels, multi-allelic mismatches) are dropped from both the returned
#' data and LD matrix. Palindromic SNPs (A/T, C/G) with an EAF in the
#' ambiguous zone (0.42-0.58, in either exposure or outcome) are also
#' dropped, since allele-identity matching alone cannot distinguish "same
#' strand" from "opposite strand, coincidentally same two letters" for
#' these -- this mirrors the equivalent safety filter in the reference
#' single-cell coloc pipeline (`harmonise_for_coloc()` in
#' `single_cell_MR_IMID/SSZ_scMR_scripts/Scripts/scMR_onek1k_1M_coloc.R`)
#' that motivated this fix.
#'
#' @param harmonised_data Data frame with `SNP`, `effect_allele.exposure`,
#'   `other_allele.exposure`, and (if present) `eaf.exposure`/
#'   `eaf.outcome` columns -- i.e. the output of [harmonise_and_filter()].
#' @param ld_matrix A named list as returned by [compute_ld_matrix()]:
#'   `ld` and `alleles`.
#' @param verbose Logical. If `TRUE`, emit a [cli::cli_inform()] summarising
#'   how many SNPs matched/flipped/dropped. Default `TRUE`.
#'
#' @return A named list with elements:
#'   - `data`: subset and reordered data frame (unmatched SNPs dropped;
#'     beta/allele columns untouched)
#'   - `ld_matrix`: subset and reordered matrix (never re-signed itself)
#'   - `ld_sign`: numeric vector (`+1`/`-1`), same length and order as
#'     `data`'s rows, to multiply into beta vectors when building coloc/
#'     SuSiE dataset objects
#'
#' @keywords internal
align_to_ld_matrix <- function(harmonised_data, ld_matrix, verbose = TRUE) {
  ld <- ld_matrix$ld
  alleles <- ld_matrix$alleles

  shared <- intersect(harmonised_data$SNP, rownames(ld))

  if (length(shared) == 0) {
    cli::cli_abort("No SNPs in common between harmonised data and LD matrix.")
  }

  data_out <- harmonised_data[match(shared, harmonised_data$SNP), ]
  ld_out <- ld[shared, shared, drop = FALSE]
  al <- alleles[match(shared, alleles$SNP), ]

  ea <- toupper(data_out$effect_allele.exposure)
  oa <- toupper(data_out$other_allele.exposure)
  # object_usage_linter false positives below: ld_a1/ld_a2 are used inside
  # dplyr::case_when()'s formula RHS/LHS, and n_match/n_flip/n_drop via cli's
  # glue-style "{var}" interpolation -- neither is traceable by lintr's
  # static analysis.
  ld_a1 <- toupper(al$ld_a1) # nolint: object_usage_linter.
  ld_a2 <- toupper(al$ld_a2) # nolint: object_usage_linter.

  orientation <- dplyr::case_when(
    ea == ld_a1 & oa == ld_a2 ~ "match",
    ea == ld_a2 & oa == ld_a1 ~ "flip",
    TRUE ~ "drop"
  )

  palindromic <- (ea == "A" & oa == "T") |
    (ea == "T" & oa == "A") |
    (ea == "C" & oa == "G") |
    (ea == "G" & oa == "C")

  eaf_exp <- data_out$eaf.exposure
  eaf_out <- data_out$eaf.outcome
  ambiguous_eaf <- (!is.na(eaf_exp) & eaf_exp >= 0.42 & eaf_exp <= 0.58) |
    (!is.na(eaf_out) & eaf_out >= 0.42 & eaf_out <= 0.58)

  orientation[palindromic & ambiguous_eaf] <- "drop"

  n_match <- sum(orientation == "match") # nolint: object_usage_linter.
  n_flip <- sum(orientation == "flip") # nolint: object_usage_linter.
  n_drop <- sum(orientation == "drop") # nolint: object_usage_linter.

  if (verbose) {
    cli::cli_inform(paste0(
      "LD alignment: {n_match} matched, {n_flip} flipped, {n_drop} dropped ",
      "(allele mismatch or ambiguous palindromic SNP vs. reference panel)."
    ))
  }

  keep <- orientation != "drop"
  ld_sign <- ifelse(orientation[keep] == "flip", -1, 1)

  data_out <- data_out[keep, ]
  ld_out <- ld_out[keep, keep, drop = FALSE]

  list(data = data_out, ld_matrix = ld_out, ld_sign = ld_sign)
}

#' Convert effect allele frequency to minor allele frequency
#'
#' @param eaf Numeric vector of effect allele frequencies.
#'
#' @return Numeric vector of minor allele frequencies.
#'
#' @keywords internal
eaf_to_maf <- function(eaf) {
  ifelse(eaf < 0.5, eaf, 1 - eaf)
}

#' Resolve sample size from multiple sources
#'
#' Attempts to determine sample size from (in order of priority):
#' 1. An explicitly provided value
#' 2. The median of a data column
#' 3. `NULL` (caller decides whether to error or warn)
#'
#' @param explicit_n Explicit sample size (numeric scalar or `NULL`).
#' @param data_column Numeric vector (e.g. `samplesize.exposure` column),
#'   or `NULL`.
#' @param label Character label for messages (e.g. `"exposure"`).
#'
#' @return Integer sample size, or `NULL` if unavailable.
#'
#' @keywords internal
resolve_sample_size <- function(
  explicit_n = NULL,
  data_column = NULL,
  label = "dataset"
) {
  if (!is.null(explicit_n)) {
    return(as.integer(explicit_n))
  }

  if (!is.null(data_column) && length(data_column) > 0) {
    n <- stats::median(data_column, na.rm = TRUE)
    if (!is.na(n)) {
      cli::cli_inform(
        "Using median sample size from {label} data: {.val {as.integer(n)}}."
      )
      return(as.integer(n))
    }
  }

  NULL
}
