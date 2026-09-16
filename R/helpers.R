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
#' Wraps [TwoSampleMR::harmonise_data()], runs the allele orientation check
#' ([check_allele_orientation()]) on the raw harmonised output, filters to
#' `mr_keep == TRUE`, and removes duplicate SNPs (keeping the first
#' occurrence).
#'
#' The allele orientation check runs *before* the `mr_keep` filter so that
#' SNPs dropped later (e.g. for a missing `beta`/`se`) still contribute their
#' allele frequencies to the verdict, and it runs on every call -- including
#' a zero-row harmonisation -- so the record returned by
#' [last_allele_check()] always describes the most recent exposure/outcome
#' pair rather than a stale one.
#'
#' @param exposure Data frame of formatted exposure data.
#' @param outcome Data frame of formatted outcome data.
#' @param allele_check One of `"error"` (default), `"warn"` or `"none"`.
#'   Passed to [check_allele_orientation()].
#' @param action `1`, `2` (default) or `3`. Passed to
#'   [TwoSampleMR::harmonise_data()]; see [validate_harmonise_action()] for
#'   what each level does.
#' @param check Logical. If `FALSE`, skip the allele orientation check
#'   entirely (nothing is recorded). Used by [run_mr()], which has already
#'   checked its instruments together with a sample of shared SNPs via
#'   [check_allele_orientation_gwas()]. Default `TRUE`.
#' @param verbose Logical. Passed to [check_allele_orientation()]. Default
#'   `FALSE`. A check that cannot reach a verdict warns regardless.
#'
#' @return A named list with elements:
#'   - `data`: the harmonised data, filtered to `mr_keep == TRUE` and
#'     deduplicated -- what the analysis runs on
#'   - `raw`: the complete, unfiltered [TwoSampleMR::harmonise_data()] output,
#'     every row and every column, including the `mr_keep`, `palindromic`,
#'     `ambiguous` and `remove` flags that explain *why* a variant was
#'     dropped (GitHub issue #17)
#'
#' @importFrom rlang .data
#' @keywords internal
harmonise_and_filter <- function(
  exposure,
  outcome,
  allele_check = c("error", "warn", "none"),
  action = 2,
  check = TRUE,
  verbose = FALSE
) {
  allele_check <- rlang::arg_match(allele_check)
  action <- validate_harmonise_action(action)

  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure,
    outcome_dat = outcome,
    action = action
  )

  if (check) {
    check_allele_orientation(
      harmonised,
      allele_check = allele_check,
      verbose = verbose,
      call = rlang::caller_env()
    )
  }

  # Guard: return 0-row frame if harmonisation produced no usable output
  # (e.g. no SNP overlap, or all SNPs removed as palindromic). This lets
  # run_mr() hit its nrow == 0 early-return rather than erroring on a missing
  # `mr_keep` column. `raw` is still handed back, because a harmonisation that
  # produced nothing usable is precisely when the caller wants to see why.
  if (nrow(harmonised) == 0 || !"mr_keep" %in% names(harmonised)) {
    return(list(data = harmonised[0, , drop = FALSE], raw = harmonised))
  }

  filtered <- harmonised |>
    dplyr::filter(.data$mr_keep == TRUE) |>
    dplyr::filter(!duplicated(.data$SNP))

  list(data = filtered, raw = harmonised)
}

#' Detect swapped effect/other alleles from harmonised allele frequencies
#'
#' Some GWAS files label their allele columns `A1`/`A2` meaning REF/ALT, with
#' `BETA` and the frequency column oriented to `A2` (EPACTS/RAREMETAL style),
#' whereas [format_gwas()] -- like PLINK, regenie and METAL -- reads `A1` as
#' the effect allele. Feeding such a file in swaps effect and other allele
#' for every variant and silently inverts every beta (GitHub issue #18). This
#' check detects that at harmonisation time, the only point where two
#' datasets coexist, so no external reference panel is needed.
#'
#' @section How the verdict is reached:
#' [TwoSampleMR::harmonise_data()] aligns non-palindromic variants purely by
#' allele letters, negating `beta.outcome` and replacing `eaf.outcome` with
#' `1 - eaf.outcome` whenever the outcome's effect allele is the exposure's
#' other allele. When one dataset's allele labels are swapped, that alignment
#' is applied to *every* variant, so after harmonisation `eaf.outcome` ends
#' up describing the other allele: `eaf.outcome ~ 1 - eaf.exposure` across
#' the set. Genuine cohort differences (e.g. ancestry) produce *scatter*
#' around `eaf.exposure`; this bug produces systematic *complementarity*.
#'
#' The statistic is therefore: among *informative* variants -- both EAFs
#' present, `palindromic == FALSE`, `remove == FALSE`, one row per SNP --
#' the proportion whose `eaf.exposure` is closer to `1 - eaf.outcome` than
#' to `eaf.outcome` (ties, e.g. `eaf.outcome == 0.5`, count as *not*
#' complementary). The check fails when that proportion exceeds `threshold`
#' and at least `min_n` informative variants were available; with fewer it
#' is `"skipped"` -- and a warning says so, because a pair that could not be
#' checked is not a pair that passed, and the two must not look alike in a
#' log (GitHub issue #21). Variants with EAF near 0.5 are equally likely to fall
#' either side, so they can only dilute the proportion towards 0.5 -- they
#' cannot cause a spurious failure, only mask a real one, which the 0.70
#' threshold tolerates.
#'
#' @section Why palindromic variants are excluded:
#' For A/T and C/G variants the strand cannot be resolved from letters, so
#' `harmonise_data()` resolves it *from the allele frequencies*: after the
#' letter-based swap it flips again if `eaf.exposure` and `eaf.outcome` sit
#' on opposite sides of 0.5. When the effect allele and the frequency are
#' both mis-assigned those two flips cancel, so a palindromic variant's
#' `eaf.outcome` always looks consistent and its harmonised beta can be
#' *identical* to the correct value even though every non-palindromic beta
#' in the same set is inverted. Including them would only dilute the
#' statistic; and once the check fails, every palindromic strand call in that
#' pair is unreliable regardless of what its beta looks like.
#'
#' @section What the check cannot tell you:
#' The comparison is symmetric: a failure means the two datasets disagree,
#' not which one is wrong. Break the tie with an independent frequency
#' reference (`plink --freq` on the LD panel; see `ref_frq` in
#' [format_gwas()]), a variant with a well-established effect direction, or
#' provenance (a dataset that has harmonised cleanly against others is not
#' the suspect).
#'
#' @param harmonised Raw output of [TwoSampleMR::harmonise_data()] (before
#'   any `mr_keep` filtering). Needs columns `SNP`, `eaf.exposure`,
#'   `eaf.outcome`, `effect_allele.exposure`, `other_allele.exposure`,
#'   `palindromic` and `remove`; otherwise the check is skipped.
#' @param allele_check One of `"error"` (default), `"warn"` or `"none"`.
#'   Controls what happens on failure; the diagnostic record is stored in
#'   every mode. `"none"` also suppresses the warning emitted when the check
#'   cannot reach a verdict.
#' @param threshold Proportion of informative variants that must be
#'   complementary for the check to fail. Default `0.70`.
#' @param min_n Minimum number of informative variants required to reach a
#'   verdict. Default `10L`.
#' @param n_sampled Integer. Number of non-instrument SNPs that
#'   [check_allele_orientation_gwas()] added to the harmonised set, for the
#'   record only. `NA` (default) when the check ran on a harmonisation that
#'   was not sampled.
#' @param verbose Logical. If `TRUE`, report a passing verdict via
#'   [cli::cli_inform()]. A verdict that could *not* be reached is always
#'   reported, as a warning, regardless of `verbose` -- unless
#'   `allele_check = "none"`, which silences it along with the check itself.
#'   Default `FALSE`.
#' @param call Environment. The calling frame reported in the condition.
#'   Default [rlang::caller_env()].
#'
#' @return The diagnostic record (see [last_allele_check()] for its
#'   structure), invisibly. The same record is stored so that
#'   [last_allele_check()] can return it -- this matters on the error path,
#'   where the return value is otherwise lost.
#'
#' @keywords internal
check_allele_orientation <- function(
  harmonised,
  allele_check = c("error", "warn", "none"),
  threshold = 0.70,
  min_n = 10L,
  n_sampled = NA_integer_,
  verbose = FALSE,
  call = rlang::caller_env()
) {
  allele_check <- rlang::arg_match(allele_check)

  col_or_na <- function(col) {
    if (col %in% names(harmonised)) unique(harmonised[[col]]) else NA_character_
  }
  # TwoSampleMR's id.* columns are random hashes for format_gwas() output;
  # the exposure/outcome name columns carry the phenotype ids users know.
  name_exp <- col_or_na("exposure") # nolint: object_usage_linter.
  name_out <- col_or_na("outcome") # nolint: object_usage_linter.

  record <- function(status, variants, n = 0L, n_comp = 0L, prop = NA_real_) {
    .mrpipeline_env$last_allele_check <- list(
      status = status,
      n = as.integer(n),
      n_complementary = as.integer(n_comp),
      prop = prop,
      threshold = threshold,
      min_n = as.integer(min_n),
      allele_check = allele_check,
      exposure = name_exp,
      outcome = name_out,
      id.exposure = col_or_na("id.exposure"),
      id.outcome = col_or_na("id.outcome"),
      n_sampled = as.integer(n_sampled),
      variants = variants
    )
    invisible(.mrpipeline_env$last_allele_check)
  }

  empty_variants <- data.frame(
    SNP = character(),
    effect_allele = character(),
    other_allele = character(),
    eaf.exposure = numeric(),
    eaf.outcome = numeric(),
    eaf.outcome_flipped = numeric(),
    palindromic = logical(),
    remove = logical(),
    informative = logical(),
    score = numeric(),
    complementary = logical(),
    stringsAsFactors = FALSE
  )

  # A skip is not a pass, and the two used to look alike in a log: an
  # informational "check skipped" scrolling past between ordinary progress
  # messages reads as silence, i.e. as no problem found (issue #21). Warn
  # instead, at the level a reader scans for, unless the user has switched
  # the check off entirely.
  unverified <- function(reason) {
    if (allele_check == "none") {
      return(invisible(NULL))
    }
    pair <- if (is.na(name_exp[1]) || is.na(name_out[1])) {
      "Allele orientation could not be checked: "
    } else {
      "Allele orientation could not be checked for {.val {name_exp}} vs {.val {name_out}}: "
    }
    cli::cli_warn(
      c(
        paste0(pair, reason, "."),
        "!" = paste0(
          "Orientation is unverified for this pair. A mis-assigned effect ",
          "allele inverts every estimate and leaves no other trace, so an ",
          "unchecked pair is not the same as a clean one."
        ),
        "i" = paste0(
          "Verify it by anchoring on a variant whose effect direction for ",
          "this trait is established beyond doubt -- no frequencies needed. ",
          "See {.code ?format_gwas}, section {.emph What does A1 mean?} ",
          "(subsection {.emph If the check cannot run})."
        ),
        "i" = paste0(
          "Silence this with {.code allele_check = \"none\"} once you have ",
          "confirmed orientation another way."
        )
      ),
      class = "mrpipeline_allele_check_unverified",
      call = call
    )
  }

  needed <- c(
    "SNP",
    "eaf.exposure",
    "eaf.outcome",
    "effect_allele.exposure",
    "other_allele.exposure",
    "palindromic",
    "remove"
  )

  if (nrow(harmonised) == 0 || !all(needed %in% names(harmonised))) {
    unverified(
      "no harmonised variants carry allele frequencies in both datasets"
    )
    return(record("skipped", empty_variants))
  }

  variants <- harmonised |>
    dplyr::filter(!is.na(.data$eaf.exposure), !is.na(.data$eaf.outcome)) |>
    dplyr::filter(!duplicated(.data$SNP)) |>
    dplyr::transmute(
      SNP = .data$SNP,
      effect_allele = .data$effect_allele.exposure,
      other_allele = .data$other_allele.exposure,
      eaf.exposure = as.numeric(.data$eaf.exposure),
      eaf.outcome = as.numeric(.data$eaf.outcome),
      eaf.outcome_flipped = 1 - .data$eaf.outcome,
      palindromic = !is.na(.data$palindromic) & .data$palindromic,
      remove = !is.na(.data$remove) & .data$remove,
      informative = !.data$palindromic & !.data$remove,
      # Positive => eaf.exposure is closer to the complement of eaf.outcome
      # than to eaf.outcome itself. Only meaningful for informative rows.
      score = abs(.data$eaf.exposure - .data$eaf.outcome) -
        abs(.data$eaf.exposure - .data$eaf.outcome_flipped),
      complementary = dplyr::if_else(.data$informative, .data$score > 0, NA)
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  n <- sum(variants$informative)
  n_comp <- sum(variants$complementary, na.rm = TRUE)
  prop <- if (n > 0) n_comp / n else NA_real_

  if (n < min_n) {
    unverified(
      "only {n} informative non-palindromic SNP{?s} with allele frequencies in both datasets (need {min_n})"
    )
    return(record("skipped", variants, n, n_comp, prop))
  }

  status <- if (prop > threshold) "fail" else "pass"
  rec <- record(status, variants, n, n_comp, prop)

  if (status == "pass") {
    if (verbose) {
      cli::cli_inform(
        "Allele orientation check passed: {n_comp}/{n} non-palindromic SNPs complementary."
      )
    }
    return(invisible(rec))
  }
  if (allele_check == "none") {
    return(invisible(rec))
  }

  offenders <- variants |>
    dplyr::filter(!is.na(.data$complementary), .data$complementary) |>
    dplyr::arrange(dplyr::desc(.data$score))
  offenders <- offenders[seq_len(min(5L, nrow(offenders))), , drop = FALSE]
  offender_lines <- sprintf(
    "%s (%s): eaf.exposure = %.3f, eaf.outcome = %.3f, 1 - eaf.outcome = %.3f",
    offenders$SNP,
    offenders$effect_allele,
    offenders$eaf.exposure,
    offenders$eaf.outcome,
    offenders$eaf.outcome_flipped
  )
  pct <- round(100 * prop) # nolint: object_usage_linter.

  msg <- c(
    paste0(
      "Possible effect/other allele mis-assignment between ",
      "{.val {name_exp}} and {.val {name_out}}: {n_comp}/{n} ({pct}%) ",
      "non-palindromic SNPs have {.field eaf.exposure} closer to ",
      "{.code 1 - eaf.outcome} than to {.field eaf.outcome}."
    ),
    "i" = paste0(
      "This usually means one of the two datasets labels A1/A2 as REF/ALT ",
      "rather than effect/other (see {.code ?format_gwas}, section ",
      "{.emph What does A1 mean?}), so effect and other alleles are swapped ",
      "and every beta is sign-inverted."
    ),
    "i" = "Worst offenders (effect allele in brackets):",
    rlang::set_names(offender_lines, rep(" ", length(offender_lines))),
    "!" = paste0(
      "Palindromic SNPs in this pair were strand-aligned from these ",
      "mis-assigned frequencies, so their alignment is unreliable even where ",
      "the beta looks unchanged."
    ),
    "!" = paste0(
      "The check cannot tell WHICH dataset is wrong. Compare each dataset's ",
      "eaf to an independent panel ({.code plink --freq} on the LD ",
      "reference; see {.arg ref_frq} in {.fn format_gwas}), check a variant ",
      "with a well-established effect direction, or suspect the dataset that ",
      "has not harmonised cleanly elsewhere."
    ),
    ">" = paste0(
      "Fix: re-run {.fn format_gwas} on the offending dataset with ",
      "{.code col_map = list(effect_allele = \"A2\", other_allele = \"A1\")} ",
      "and re-run the analysis."
    ),
    "i" = paste0(
      "Inspect the full record with {.fn last_allele_check}. Downgrade with ",
      "{.code allele_check = \"warn\"} or disable with ",
      "{.code allele_check = \"none\"}."
    )
  )

  if (allele_check == "error") {
    cli::cli_abort(msg, class = "mrpipeline_allele_check_error", call = call)
  }
  cli::cli_warn(msg, class = "mrpipeline_allele_check_warning", call = call)
  invisible(rec)
}

#' Allele orientation check on instruments plus a sample of shared SNPs
#'
#' [run_mr()]'s instrument set is often too small for
#' [check_allele_orientation()] to reach a verdict (a cis-MR may have three
#' instruments; the check needs ten informative non-palindromic SNPs). But
#' `run_mr()` holds the *full* exposure and outcome GWAS before it narrows
#' the outcome to instrument rsIDs, so this helper harmonises the
#' instruments together with up to `n_sample` additional SNPs shared by the
#' two datasets and runs the verdict on that larger set. The extra SNPs are
#' taken at evenly spaced positions through the sorted shared rsIDs --
#' deterministic, so results are reproducible and the caller's RNG state is
#' untouched, and effectively random with respect to genomic position.
#'
#' Cost is dominated by the rsID intersection (about 0.5 s for a 200k-SNP
#' exposure against a 10M-row outcome); formatting and harmonising ~1000
#' SNPs takes a few milliseconds.
#'
#' The check is skipped (with a `"skipped"` record, and a warning that
#' orientation is unverified) when the exposure lacks `SNP`/`eaf.exposure` or
#' the outcome lacks the [format_gwas()] outcome columns, or when nothing
#' overlaps.
#'
#' @param exposure Data frame of TwoSampleMR-formatted exposure data (the
#'   full dataset passed to `run_mr()`, not just the instruments).
#' @param outcome Data frame in [format_gwas()] outcome format (`rsids`,
#'   `effect_allele`, `other_allele`, `beta`, `se`, `eaf`, ...).
#' @param instrument_snps Character vector of instrument rsIDs to always
#'   include in the check set.
#' @param allele_check One of `"error"` (default), `"warn"` or `"none"`.
#' @param n_sample Integer. Maximum number of non-instrument shared SNPs to
#'   add. Default `1000L`.
#' @param action `1`, `2` (default) or `3`. Passed to
#'   [TwoSampleMR::harmonise_data()], so that the check describes the same
#'   harmonisation the analysis will use. The verdict itself is unaffected:
#'   `action` gates only the frequency-based second flip applied to
#'   palindromic variants, which the check excludes anyway.
#' @param verbose Logical. Passed to [check_allele_orientation()]. A check
#'   that cannot reach a verdict warns regardless.
#' @param call Environment reported in the condition. Default
#'   [rlang::caller_env()].
#'
#' @return The diagnostic record, invisibly (see [last_allele_check()]).
#'
#' @keywords internal
check_allele_orientation_gwas <- function(
  exposure,
  outcome,
  instrument_snps,
  allele_check = c("error", "warn", "none"),
  n_sample = 1000L,
  action = 2,
  verbose = FALSE,
  call = rlang::caller_env()
) {
  allele_check <- rlang::arg_match(allele_check)
  action <- validate_harmonise_action(action)

  skipped <- function() {
    check_allele_orientation(
      data.frame(),
      allele_check = allele_check,
      verbose = verbose,
      call = call
    )
  }

  needed_exp <- c("SNP", "eaf.exposure")
  needed_out <- c("rsids", "effect_allele", "other_allele", "beta", "se", "eaf")
  if (
    !all(needed_exp %in% names(exposure)) ||
      !all(needed_out %in% names(outcome))
  ) {
    return(skipped())
  }

  shared <- setdiff(intersect(exposure$SNP, outcome$rsids), instrument_snps)
  shared <- shared[!is.na(shared)]
  if (length(shared) > n_sample) {
    shared <- sort(shared)[round(seq(1, length(shared), length.out = n_sample))]
  }
  check_snps <- unique(c(instrument_snps, shared))

  exposure_sub <- exposure[exposure$SNP %in% check_snps, , drop = FALSE]
  outcome_sub <- outcome[outcome$rsids %in% check_snps, , drop = FALSE]
  if (nrow(exposure_sub) == 0 || nrow(outcome_sub) == 0) {
    return(skipped())
  }

  # Sampled SNPs are only used for this check, so format_data()'s warnings
  # about rows missing beta/se (excluded from "the MR tests") are noise here.
  outcome_data <- suppressWarnings(suppressMessages(TwoSampleMR::format_data(
    outcome_sub,
    type = "outcome",
    phenotype_col = "phenotype",
    header = TRUE,
    snp_col = "rsids",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    beta_col = "beta",
    se_col = "se",
    samplesize_col = "n",
    pval_col = "pval",
    pos_col = "pos",
    chr_col = "chr",
    log_pval = FALSE
  )))

  harmonised <- suppressMessages(TwoSampleMR::harmonise_data(
    exposure_dat = exposure_sub,
    outcome_dat = outcome_data,
    action = action
  ))

  check_allele_orientation(
    harmonised,
    allele_check = allele_check,
    n_sampled = length(shared),
    verbose = verbose,
    call = call
  )
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

  # ieugwasr builds these names from a read.table() of the .bim, which turns
  # a lone "T" allele into logical TRUE (a single-SNP matrix, or every A1 =
  # T). Left as "TRUE" it matches no exposure allele and the SNP is dropped
  # at alignment.
  alleles <- data.frame(
    SNP = rsid,
    ld_a1 = stringr::str_replace(parts[, 2], "^TRUE$", "T"),
    ld_a2 = stringr::str_replace(parts[, 3], "^TRUE$", "T"),
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

#' Read a delimited file, checking compression support first
#'
#' Thin wrapper around [data.table::fread()] that fails early, with an
#' actionable message, when `path` is compressed but the `R.utils` package
#' (which `fread()` needs to decompress `.gz`/`.bz2` files) is unavailable.
#' `R.utils` is only a `Suggests` of `data.table`, so declaring it in
#' `mrpipeline`'s `Imports` is what actually guarantees the gzipped inputs
#' that GWAS summary statistics are routinely distributed as (issue #20).
#' This guard is a backstop for a broken library rather than an expected
#' failure mode, and keeps the error at the top of the call stack instead of
#' deep inside `fread()`.
#'
#' @param path Character scalar file path.
#' @param ... Passed to [data.table::fread()].
#'
#' @return The value of [data.table::fread()].
#'
#' @keywords internal
fread_file <- function(path, ...) {
  check_gz_support(path)
  data.table::fread(path, ...)
}

#' Check that compressed files can be read
#'
#' @param path Character scalar file path. Non-character input is ignored.
#' @param has_rutils Logical, or `NULL` (default) to look `R.utils` up with
#'   `requireNamespace()`. Only consulted when `path` is compressed. Exposed
#'   so tests can exercise the failure branch without uninstalling the
#'   package.
#'
#' @return `invisible(TRUE)` if `path` is readable, otherwise aborts.
#'
#' @keywords internal
check_gz_support <- function(path, has_rutils = NULL) {
  compressed <- is.character(path) &&
    length(path) == 1L &&
    stringr::str_detect(path, "\\.(gz|bz2)$")
  if (!compressed) {
    return(invisible(TRUE))
  }
  if (is.null(has_rutils)) {
    # decompressFile() is the entry point data.table::fread() reaches for to
    # open .gz/.bz2. Naming it explicitly also makes the R.utils dependency
    # visible to R CMD check, which does not count requireNamespace() as use
    # of a declared Import. The `&&` short-circuits, so `R.utils::` is only
    # evaluated once the namespace is known to be available.
    has_rutils <- requireNamespace("R.utils", quietly = TRUE) &&
      is.function(R.utils::decompressFile)
  }
  if (!has_rutils) {
    cli::cli_abort(
      c(
        "Cannot read compressed file {.path {path}}.",
        "x" = "{.pkg data.table} needs {.pkg R.utils} to decompress
               {.field .gz} / {.field .bz2} files, and it is not installed.",
        "i" = "Install it with {.run install.packages(\"R.utils\")}, or
               supply an uncompressed file."
      )
    )
  }
  invisible(TRUE)
}

#' Validate a TwoSampleMR harmonisation action level
#'
#' [TwoSampleMR::harmonise_data()] accepts `action` as a vector, applying a
#' different level per outcome. `run_mr()` and `run_coloc()` handle exactly
#' one outcome, so a vector here is a mistake worth catching rather than
#' silently recycling.
#'
#' The three levels:
#'
#' | `action` | Behaviour |
#' |---|---|
#' | 1 | Assume all alleles are on the forward strand: no frequency-based flip |
#' | 2 | Infer the positive strand, resolving palindromes from allele frequencies (default) |
#' | 3 | As 2, but set `mr_keep = FALSE` for every palindromic, ambiguous or incompatible SNP |
#'
#' Only the palindrome handling differs. The letter-based alignment of
#' non-palindromic variants -- negating `beta.outcome` and replacing
#' `eaf.outcome` with `1 - eaf.outcome` when the outcome's effect allele is
#' the exposure's other allele -- happens at every level, which is why
#' [check_allele_orientation()]'s verdict does not depend on `action`.
#'
#' @param action Value to validate.
#'
#' @return `action`, unchanged. It is deliberately not coerced to integer:
#'   [TwoSampleMR::harmonise_data()] embeds `action` as a column in its
#'   output, so coercing `2` to `2L` would make `mrpipeline`'s harmonised
#'   frame differ from a plain `harmonise_data()` call on the same data by
#'   the storage mode of that column alone.
#'
#' @keywords internal
validate_harmonise_action <- function(action) {
  ok <- is.numeric(action) &&
    length(action) == 1L &&
    !is.na(action) &&
    action %in% 1:3
  if (!ok) {
    cli::cli_abort(
      c(
        "{.arg action} must be a single value of {.val {1:3}}.",
        "x" = "Got {.val {action}}.",
        "i" = paste0(
          "{.val {1}} assumes the forward strand, {.val {2}} resolves ",
          "palindromes from allele frequencies, {.val {3}} drops every ",
          "palindromic, ambiguous or incompatible SNP."
        )
      ),
      call = rlang::caller_env()
    )
  }
  action
}

#' Summarise what happened during harmonisation
#'
#' Counts, from the unfiltered [TwoSampleMR::harmonise_data()] output, how
#' many variants were carried forward and why the rest were not. Used by
#' [summary.mr_result()] and [summary.coloc_result()].
#'
#' A variant can carry more than one flag -- an ambiguous variant is by
#' definition palindromic -- so the reason counts are *not* a partition of
#' `n_dropped` and must not be presented as one. Which flags actually cost a
#' variant its place depends on the harmonisation action: `remove` at every
#' level, `ambiguous` from level 2, `palindromic` only at level 3 (see
#' [validate_harmonise_action()]). `n_incomplete` covers the separate case of
#' a variant dropped by `harmonise_data()` for missing beta/se rather than
#' for any allele problem.
#'
#' @param raw The `raw` element of [harmonise_and_filter()]'s return value.
#'
#' @return A named list of integers: `n_candidates`, `n_kept`, `n_dropped`,
#'   `n_duplicate`, `n_palindromic`, `n_ambiguous`, `n_incompatible` and
#'   `n_incomplete`. All zero when `raw` is empty or lacks the flag columns.
#'
#' @keywords internal
harmonisation_summary <- function(raw) {
  zeros <- list(
    n_candidates = 0L,
    n_kept = 0L,
    n_dropped = 0L,
    n_duplicate = 0L,
    n_palindromic = 0L,
    n_ambiguous = 0L,
    n_incompatible = 0L,
    n_incomplete = 0L
  )
  if (is.null(raw) || !is.data.frame(raw) || nrow(raw) == 0) {
    return(zeros)
  }

  flag <- function(col) {
    if (col %in% names(raw)) {
      !is.na(raw[[col]]) & raw[[col]]
    } else {
      rep(FALSE, nrow(raw))
    }
  }
  mr_keep <- flag("mr_keep")
  palindromic <- flag("palindromic")
  ambiguous <- flag("ambiguous")
  remove <- flag("remove")

  duplicate <- if ("SNP" %in% names(raw)) {
    mr_keep & duplicated(raw$SNP)
  } else {
    rep(FALSE, nrow(raw))
  }

  list(
    n_candidates = nrow(raw),
    n_kept = sum(mr_keep & !duplicate),
    n_dropped = sum(!mr_keep),
    n_duplicate = sum(duplicate),
    n_palindromic = sum(palindromic),
    n_ambiguous = sum(ambiguous),
    n_incompatible = sum(remove),
    # Dropped by harmonise_data() for missing beta/se rather than for any
    # allele problem -- invisible in the three flags above.
    n_incomplete = sum(!mr_keep & !remove & !ambiguous & !palindromic)
  ) |>
    lapply(as.integer)
}

#' Print a harmonisation breakdown
#'
#' Shared by [summary.mr_result()] and [summary.coloc_result()]. Prints
#' nothing when there is no harmonisation record to describe (older result
#' objects, or a run that never reached harmonisation).
#'
#' @param harmonisation The `harmonisation` field of an `mr_result` or
#'   `coloc_result`.
#'
#' @return `invisible(NULL)`, called for its output.
#'
#' @keywords internal
print_harmonisation_summary <- function(harmonisation) {
  h <- harmonisation_summary(harmonisation)
  if (h$n_candidates == 0L) {
    return(invisible(NULL))
  }

  cli::cli_h2("Harmonisation")
  cli::cli_bullets(c(
    "*" = "{h$n_candidates} candidate SNP{?s} -> {h$n_kept} kept, {h$n_dropped} dropped"
  ))

  # Reason counts overlap (every ambiguous variant is palindromic), so they
  # are listed rather than summed, and only non-zero ones are shown.
  reasons <- c(
    "palindromic" = h$n_palindromic,
    "ambiguous" = h$n_ambiguous,
    "incompatible alleles" = h$n_incompatible,
    "incomplete beta/se" = h$n_incomplete
  )
  reasons <- reasons[reasons > 0]
  if (length(reasons) > 0) {
    lines <- paste(reasons, names(reasons)) # nolint: object_usage_linter.
    cli::cli_bullets(c("*" = "Flagged: {lines}"))
  }
  if (h$n_duplicate > 0) {
    cli::cli_bullets(c(
      "*" = "{h$n_duplicate} duplicate SNP row{?s} dropped after filtering"
    ))
  }
  invisible(NULL)
}

#' Warn that a requested method has no LD-corrected form
#'
#' Called by [run_mr()] for every `$results`-producing method that runs while
#' `ld_correct = TRUE` but has no correlated implementation (registry
#' `ld_correctable = FALSE`), so `ld_correct` is never silently ignored: it is
#' either applied, or visibly not applied. Diagnostics (`steiger`,
#' `pleiotropy`, `heterogeneity`, `loo`) do not warn -- they produce no
#' estimate row.
#'
#' @param method The shortcut or raw method name, as the user passed it.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
warn_no_ld_correction <- function(method) {
  cli::cli_warn(
    "{.val {method}} has no LD-corrected form; running uncorrected."
  )
  invisible(NULL)
}
