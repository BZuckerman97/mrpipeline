#' Perform Mendelian randomisation analysis
#'
#' Runs MR with automatic instrument selection (cis-MR, genome-wide, or manual)
#' and optional sensitivity analyses. Returns an `mr_result` S3 object.
#'
#' @section Instrument selection modes:
#' Exactly one of three modes is used, determined by the combination of
#' `instruments` and `instrument_region`:
#'
#' - **Cis-MR** (`instrument_region` provided, `instruments = NULL`): filters
#'   `exposure` to the cis region defined by `instrument_region` +/- `window`,
#'   applies `pval_thresh`, then LD-clumps.
#' - **Genome-wide** (`instrument_region = NULL`, `instruments = NULL`): filters
#'   `exposure` by `pval_thresh` only, then LD-clumps.
#' - **Manual** (`instruments` provided): uses the supplied rsIDs directly.
#'   `instruments_strict` controls whether missing IDs are an error or warning.
#'
#' @section Method dispatch:
#' Every method is looked up in the package's method registry -- see
#' [mr_methods()] and the *Available methods* section below for what each
#' shortcut runs, whether it can be LD-corrected, and how many instruments
#' it needs. Dispatch depends on the number of instruments after clumping:
#' - 1 SNP: Wald ratio only; every multi-SNP method is skipped
#' - 2+ SNPs: `ivw_random`, `ivw_fixed`, `conmix`, `heterogeneity` and any
#'   raw TwoSampleMR method are attempted; `egger`, `weighted_median`,
#'   `presso`, `pleiotropy` and `loo` require >= 3 SNPs
#'
#' Raw TwoSampleMR methods (`mr_*` names from
#' `TwoSampleMR::mr_method_list()$obj` that have no shortcut, e.g.
#' `"mr_raps"`) are dispatched via `TwoSampleMR::mr()` under TwoSampleMR's
#' own label; errors are caught and reported in `$methods_skipped`. Names
#' that are the engine behind a shortcut (`mr_ivw`, `mr_ivw_fe`,
#' `mr_egger_regression`, `mr_weighted_median`) are refused with a pointer
#' to the shortcut, so the same estimator cannot enter unlabelled and
#' uncorrected.
#'
#' ConMix reports `se = NA`: `MendelianRandomization::mr_conmix()` returns a
#' confidence interval that may be asymmetric or multi-modal rather than a
#' standard error, so the `or_lci95`/`or_uci95` columns are `NA` for it too.
#'
#' When `"egger"` is in `methods` and there are >= 3 instruments,
#' `TwoSampleMR::mr_pleiotropy_test()` (the Egger intercept test) is always
#' run automatically and its result stored in `$pleiotropy`. You do not need
#' to add `"pleiotropy"` to `methods` separately. The `"pleiotropy"` shortcut
#' remains available for running the intercept test without Egger.
#'
#' @section LD correction:
#' `ld_correct = TRUE` computes a signed LD matrix for the instruments from
#' `bfile` and re-orients it to the exposure's effect alleles. Instruments
#' absent from the reference panel, or with ambiguous palindromic alleles,
#' are dropped at that step, so an LD-corrected run can have fewer
#' instruments than the same call uncorrected.
#'
#' Methods with a correlated form -- `ivw_random`, `ivw_fixed` and `egger`
#' -- are then fitted by generalised least squares through
#' `MendelianRandomization` with `correl = TRUE`: the weight matrix is
#' `diag(se_y) %*% R %*% diag(se_y)` in place of `diag(se_y^2)`, so two
#' instruments in LD are no longer counted as two independent looks at the
#' causal effect. The IVW estimator is pinned explicitly (`model = "random"`
#' for `ivw_random`, `"fixed"` for `ivw_fixed`) so that each shortcut means
#' the same thing at every instrument count; `MendelianRandomization`'s own
#' default would switch to fixed effects below 4 instruments.
#'
#' Every other `$results` method has no correlated form and runs on the
#' uncorrected data: a warning names each such method, and its row carries
#' `ld_corrected = FALSE`. The diagnostics come from the correlated fits
#' too: `$heterogeneity` holds the generalised Cochran Q for correlated
#' instruments (`Q = r' O^-1 r` on the GLS residuals, from
#' `mr_ivw()@Heter.Stat` and `mr_egger()@Heter.Stat`), `$pleiotropy` the
#' correlated Egger intercept (`mr_egger()@Intercept`), and `$loo` a
#' per-SNP correlated random-effects refit (a block-inverse update, so it
#' stays O(n^3)). Each of those frames carries an `ld_corrected` column on
#' both arms. `steiger` is the one thing left as-is: the Steiger direction
#' test compares per-SNP r^2 values and involves no weight matrix.
#' `ld_correct` is never silently ignored -- it is either applied, or
#' visibly not applied. To compare corrected and uncorrected estimates, call
#' `run_mr()` twice and pass both results to [forest_plot()] as a named
#' list: one `mr_result` is always one instrument set under one weight
#' matrix.
#'
#' If the GLS weight matrix is near-singular (reciprocal condition number
#' below `1e-10`) `run_mr()` warns that every LD-corrected estimate is
#' unstable. The usual causes are identical or near-identical instruments
#' (r^2 ~ 1, from absent clumping or a manual set with a duplicated
#' variant) and more instruments than reference-panel individuals, which
#' makes the sample correlation matrix singular; clump more stringently or
#' drop the duplicate.
#'
#' Random effects are multiplicative: the standard error is inflated by
#' `max(RSE, 1)`, never deflated, so when the instruments are under-dispersed
#' the random- and fixed-effect results coincide exactly.
#'
#' @eval rd_method_table()
#'
#' @param exposure Data frame of formatted exposure data (output of
#'   [TwoSampleMR::format_data()] or `format_pqtl_*()` functions).
#' @param exposure_id Character. Identifier for the exposure (e.g. protein
#'   name).
#' @param outcome Data frame of outcome summary statistics with standardised
#'   columns: `rsids`, `chr`, `pos`, `beta`, `se`, `eaf`, `pval`, `n`,
#'   `effect_allele`, `other_allele`. Formatted internally via
#'   [TwoSampleMR::format_data()].
#' @param outcome_id Character. Identifier for the outcome (e.g. disease name).
#' @param instrument_region List with elements `chromosome`, `start`, `end`
#'   defining the cis region. `NULL` for genome-wide or manual mode.
#' @param window Integer. Window (in bp) to extend either side of
#'   `instrument_region`. Default `100000L`.
#' @param pval_thresh Numeric. P-value threshold for instrument selection.
#'   Default `5e-8`.
#' @param rsq_thresh Numeric. R-squared clumping threshold. Default `0.001`.
#' @param bfile Character. Path to PLINK bfile prefix for local LD operations.
#'   Required when `ld_correct = TRUE`.
#' @param plink_bin Character. Path to PLINK binary. Auto-detected if `NULL`.
#' @param pop Character. Population for API-based LD clumping. Default `"EUR"`.
#' @param instruments Character vector of rsIDs for manual instrument mode, or
#'   `NULL`.
#' @param instruments_strict Logical. If `TRUE`, error when manual instruments
#'   are missing from exposure data. If `FALSE`, warn. Default `FALSE`.
#' @param exclude_regions Data frame with columns `chr`, `start`, `end` defining
#'   genomic regions to exclude instruments from, or `NULL`. For example, to
#'   exclude the MHC region: `data.frame(chr = "6", start = 26e6, end = 34e6)`.
#' @param methods Character vector of methods to run. Named shortcuts (see
#'   [mr_methods()] and the *Available methods* section): `"ivw_random"`,
#'   `"ivw_fixed"`, `"egger"`, `"weighted_median"`, `"presso"`, `"conmix"`,
#'   `"steiger"`, `"pleiotropy"`, `"heterogeneity"`, `"loo"`. You may also
#'   pass any raw method name from `TwoSampleMR::mr_method_list()$obj` that
#'   has no shortcut (e.g. `"mr_raps"`, `"mr_weighted_mode"`); raw names are
#'   never LD-corrected. The former shortcuts `"ivw"` and `"ivw_fe"` are
#'   accepted with a warning and mapped to `"ivw_random"` and `"ivw_fixed"`.
#' @param ld_correct Logical. Fit `ivw_random`, `ivw_fixed` and `egger` by
#'   GLS with the instruments' LD matrix (see the *LD correction* section).
#'   Requires `bfile`. Default `FALSE`.
#' @param exposure_n Numeric. Exposure sample size. If `NULL`, inferred from
#'   `samplesize.exposure` column.
#' @param presso_n_dist Integer. Number of distributions for MR-PRESSO. Default
#'   `1000`.
#' @param plink_threads Integer. Number of threads for PLINK. `NULL` (default)
#'   lets PLINK auto-detect. Read from `getOption("mrpipeline.plink_threads")`
#'   or the `MRPIPELINE_PLINK_THREADS` environment variable via
#'   [plink_option()].
#' @param plink_memory Integer. Memory limit in MB for PLINK. `NULL` (default)
#'   lets PLINK auto-detect. Read from `getOption("mrpipeline.plink_memory")`
#'   or the `MRPIPELINE_PLINK_MEMORY` environment variable via
#'   [plink_option()].
#' @param allele_check Character. What to do when the allele orientation
#'   check finds that effect/other alleles look swapped between exposure and
#'   outcome -- the signature of a GWAS file whose `A1`/`A2` mean REF/ALT,
#'   which silently inverts every beta (see [format_gwas()], section *What
#'   does A1 mean?*). `"error"` (default) aborts, `"warn"` warns and
#'   continues, `"none"` runs the analysis regardless. The check harmonises
#'   the instruments together with up to 1000 further SNPs shared by the two
#'   datasets, so it works even for a cis-MR with a handful of instruments;
#'   it is skipped when fewer than 10 informative non-palindromic SNPs carry
#'   both allele frequencies. The full record is available afterwards from
#'   [last_allele_check()] in every mode.
#' @param harmonise_action `1`, `2` (default) or `3`, passed to
#'   [TwoSampleMR::harmonise_data()]. `1` assumes every allele is on the
#'   forward strand; `2` infers the positive strand, resolving palindromic
#'   variants from their allele frequencies; `3` additionally drops every
#'   palindromic, ambiguous or incompatible SNP. Only palindrome handling
#'   differs -- non-palindromic variants are aligned by allele letter at all
#'   three levels. Reach for `3` when the frequencies that level `2` relies
#'   on cannot be trusted: a failed allele orientation check makes every
#'   palindromic strand call in that pair unreliable, and a dataset without
#'   allele frequencies gives level `2` nothing to resolve them with.
#' @param verbose Logical. If `TRUE`, emit informational messages via
#'   [cli::cli_inform()]. Warnings and errors are always emitted regardless.
#'   Default `TRUE`.
#'
#' @return An `mr_result` object. Check `result$status` for `"success"` vs
#'   failure reasons. Each `$results` row states the estimator that produced
#'   it: `method` (the label), `model` (`"random"`, `"fixed"`, or `NA` where
#'   the distinction does not apply) and `ld_corrected` (whether the LD
#'   matrix was used for that fit), alongside `nsnp`, `b`, `se`, `pval` and
#'   the `or`, `or_lci95`, `or_uci95` columns from
#'   [TwoSampleMR::generate_odds_ratios()]. The `$timing` field contains a
#'   named numeric vector of elapsed seconds for each major step.
#'
#' @examples
#' \dontrun{
#' # Cis-MR using bundled CD40/Sjogren's data
#' bfile <- sub("\\.bed$", "", system.file("extdata", "ld_ref.bed", package = "mrpipeline"))
#' result <- run_mr(
#'   exposure = cd40_exposure,
#'   exposure_id = "CD40",
#'   outcome = sjogren_outcome,
#'   outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = bfile,
#'   methods = c("ivw_random", "egger", "weighted_median")
#' )
#' result
#' summary(result)
#'
#' # LD-corrected: both IVW estimators by GLS; weighted median warns and
#' # runs uncorrected. Compare arms by passing both results to forest_plot().
#' corrected <- run_mr(
#'   exposure = cd40_exposure,
#'   exposure_id = "CD40",
#'   outcome = sjogren_outcome,
#'   outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = bfile,
#'   ld_correct = TRUE,
#'   methods = c("ivw_random", "ivw_fixed", "egger", "weighted_median")
#' )
#' corrected$results[, c("method", "model", "ld_corrected", "b", "se")]
#' forest_plot(list("Uncorrected" = result, "LD-corrected" = corrected))
#' }
#'
#' @seealso [mr_methods()] for the table of methods and what each supports.
#' @export
run_mr <- function(
  exposure,
  exposure_id,
  outcome,
  outcome_id,
  instrument_region = NULL,
  window = 100000L,
  pval_thresh = 5e-8,
  rsq_thresh = 0.001,
  bfile = NULL,
  plink_bin = NULL,
  pop = "EUR",
  instruments = NULL,
  instruments_strict = FALSE,
  exclude_regions = NULL,
  methods = c(
    "ivw_random",
    "egger",
    "weighted_median",
    "presso",
    "conmix",
    "steiger"
  ),
  ld_correct = FALSE,
  exposure_n = NULL,
  presso_n_dist = 1000,
  plink_threads = plink_option("threads"),
  plink_memory = plink_option("memory"),
  allele_check = c("error", "warn", "none"),
  harmonise_action = 2,
  verbose = TRUE
) {
  # --- Validate arguments ---------------------------------------------------

  if (ld_correct && is.null(bfile)) {
    cli::cli_abort("{.arg bfile} is required when {.code ld_correct = TRUE}.")
  }

  allele_check <- rlang::arg_match(allele_check)
  harmonise_action <- validate_harmonise_action(harmonise_action)

  # Every method name is resolved against the registry (see mr_methods()):
  # it is the only source of shortcut names, minimum instrument counts and
  # LD-correctability used anywhere in this function.
  registry <- mr_method_registry()
  shortcut_methods <- registry$shortcut[!is.na(registry$shortcut)]

  # The former "ivw"/"ivw_fe" names are normalised here so nothing downstream
  # ever sees them. Warned on every call: a pipeline looping over exposures
  # should keep seeing it.
  aliases <- c(ivw = "ivw_random", ivw_fe = "ivw_fixed")
  is_alias <- methods %in% names(aliases)
  if (any(is_alias)) {
    old <- methods[is_alias] # nolint: object_usage_linter.
    new <- unname(aliases[methods[is_alias]]) # nolint: object_usage_linter.
    cli::cli_warn(c(
      "{cli::qty(length(old))}Method shortcut{?s} {.val {old}} {?is/are} deprecated; using {.val {new}}.",
      "i" = paste0(
        "{.val ivw_random} is multiplicative random effects at every ",
        "instrument count, including with {.code ld_correct = TRUE} ",
        "(previously fixed effects below 4 instruments)."
      )
    ))
    methods[is_alias] <- new
  }
  methods <- unique(methods)

  # Raw TwoSampleMR names that are the engine behind a shortcut are refused:
  # accepting them would let the same estimator in unlabelled, with no
  # `model` and never LD-corrected. Derived from the registry's engine
  # column, so a new shortcut shadows its raw name automatically.
  tsm_all <- TwoSampleMR::mr_method_list()$obj
  raw_engine <- stringr::str_remove(registry$engine, "^TwoSampleMR::")
  shadowed <- stats::setNames(
    registry$shortcut[!is.na(registry$shortcut) & raw_engine %in% tsm_all],
    raw_engine[!is.na(registry$shortcut) & raw_engine %in% tsm_all]
  )
  hit <- methods[methods %in% names(shadowed)]
  if (length(hit) > 0) {
    use <- unname(shadowed[hit]) # nolint: object_usage_linter.
    cli::cli_abort(c(
      "{cli::qty(length(hit))}{.val {hit}} {?is/are} the function{?s} behind the {.val {use}} shortcut{?s}.",
      "i" = paste0(
        "Use the shortcut instead, so that {.arg ld_correct} and the ",
        "{.field model} column apply to it."
      )
    ))
  }

  tsm_available <- setdiff(tsm_all, c("mr_wald_ratio", names(shadowed)))
  unknown_methods <- setdiff(methods, c(shortcut_methods, tsm_available))
  if (length(unknown_methods) > 0) {
    cli::cli_abort(
      c(
        "Unknown method{?s}: {.val {unknown_methods}}.",
        "i" = "See {.fn mr_methods} for the named shortcuts.",
        "i" = "Or pass any name from {.code TwoSampleMR::mr_method_list()}."
      )
    )
  }

  if (!is.null(exclude_regions)) {
    validate_exclude_regions(exclude_regions)
  }

  params <- list(
    exposure_id = exposure_id,
    outcome_id = outcome_id,
    instrument_region = instrument_region,
    window = window,
    pval_thresh = pval_thresh,
    rsq_thresh = rsq_thresh,
    bfile = bfile,
    pop = pop,
    instruments = instruments,
    instruments_strict = instruments_strict,
    exclude_regions = exclude_regions,
    methods = methods,
    ld_correct = ld_correct,
    exposure_n = exposure_n,
    presso_n_dist = presso_n_dist,
    allele_check = allele_check,
    harmonise_action = harmonise_action
  )

  timing <- numeric(0)

  # --- Instrument selection -------------------------------------------------

  t0 <- proc.time()[["elapsed"]]

  if (!is.null(instruments)) {
    # Manual mode
    if (verbose) {
      cli::cli_inform("Using {length(instruments)} manual instrument{?s}.")
    }
    exposure_iv <- exposure[exposure$SNP %in% instruments, ]

    missing <- setdiff(instruments, exposure_iv$SNP)
    if (length(missing) > 0) {
      msg <- "{length(missing)} instrument{?s} not found in exposure data: {.val {missing}}"
      if (instruments_strict) {
        cli::cli_abort(msg)
      } else {
        cli::cli_warn(msg)
      }
    }

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn("No manual instruments found in exposure data.")
      timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0
      return(new_mr_result(
        status = "no_instruments",
        status_reason = "No manual instruments found in exposure data",
        params = params,
        timing = timing
      ))
    }
  } else if (!is.null(instrument_region)) {
    # Cis-MR mode
    if (verbose) {
      cli::cli_inform(paste0(
        "Cis-MR mode: chr{instrument_region$chromosome}:",
        "{instrument_region$start}-{instrument_region$end} (+/- {window}bp)."
      ))
    }

    exposure_iv <- exposure |>
      dplyr::filter(
        as.character(.data$chr.exposure) ==
          as.character(instrument_region$chromosome),
        .data$pos.exposure >= (instrument_region$start - window),
        .data$pos.exposure <= (instrument_region$end + window)
      ) |>
      dplyr::filter(.data$pval.exposure < pval_thresh)

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn(
        "No significant instruments in cis region for {.val {exposure_id}}."
      )
      timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0
      return(new_mr_result(
        status = "no_instruments",
        status_reason = paste0(
          "No significant instruments in cis region for '",
          exposure_id,
          "'"
        ),
        params = params,
        timing = timing
      ))
    }

    # Clump
    clump_dat <- data.frame(
      rsid = exposure_iv$SNP,
      pval = exposure_iv$pval.exposure,
      id = exposure_iv$id.exposure,
      beta = exposure_iv$beta.exposure,
      se = exposure_iv$se.exposure,
      stringsAsFactors = FALSE
    )
    clumped <- clump_instruments(
      dat = clump_dat,
      rsq_thresh = rsq_thresh,
      bfile = bfile,
      plink_bin = plink_bin,
      pop = pop,
      plink_threads = plink_threads,
      plink_memory = plink_memory
    )
    exposure_iv <- exposure_iv[exposure_iv$SNP %in% clumped$rsid, ]

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn(
        "No instruments remaining after clumping for {.val {exposure_id}}."
      )
      timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0
      return(new_mr_result(
        status = "no_instruments",
        status_reason = paste0(
          "No instruments remaining after clumping for '",
          exposure_id,
          "'"
        ),
        params = params,
        timing = timing
      ))
    }
  } else {
    # Genome-wide mode
    if (verbose) {
      cli::cli_inform(
        "Genome-wide mode: selecting instruments at p < {pval_thresh}."
      )
    }

    exposure_iv <- exposure |>
      dplyr::filter(.data$pval.exposure < pval_thresh)

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn(
        "No genome-wide significant instruments for {.val {exposure_id}}."
      )
      timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0
      return(new_mr_result(
        status = "no_instruments",
        status_reason = paste0(
          "No genome-wide significant instruments for '",
          exposure_id,
          "'"
        ),
        params = params,
        timing = timing
      ))
    }

    # Clump
    clump_dat <- data.frame(
      rsid = exposure_iv$SNP,
      pval = exposure_iv$pval.exposure,
      id = exposure_iv$id.exposure,
      beta = exposure_iv$beta.exposure,
      se = exposure_iv$se.exposure,
      stringsAsFactors = FALSE
    )
    clumped <- clump_instruments(
      dat = clump_dat,
      rsq_thresh = rsq_thresh,
      bfile = bfile,
      plink_bin = plink_bin,
      pop = pop,
      plink_threads = plink_threads,
      plink_memory = plink_memory
    )
    exposure_iv <- exposure_iv[exposure_iv$SNP %in% clumped$rsid, ]

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn(
        "No instruments remaining after clumping for {.val {exposure_id}}."
      )
      timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0
      return(new_mr_result(
        status = "no_instruments",
        status_reason = paste0(
          "No instruments remaining after clumping for '",
          exposure_id,
          "'"
        ),
        params = params,
        timing = timing
      ))
    }
  }

  timing[["instrument_selection"]] <- proc.time()[["elapsed"]] - t0

  # --- Region exclusion -----------------------------------------------------

  t0 <- proc.time()[["elapsed"]]

  if (!is.null(exclude_regions) && "chr.exposure" %in% colnames(exposure_iv)) {
    in_excluded <- rep(FALSE, nrow(exposure_iv))
    for (i in seq_len(nrow(exclude_regions))) {
      chr_match <- as.character(exposure_iv$chr.exposure) ==
        as.character(exclude_regions$chr[i])
      pos_in_region <- exposure_iv$pos.exposure >= exclude_regions$start[i] &
        exposure_iv$pos.exposure <= exclude_regions$end[i]
      in_excluded <- in_excluded | (chr_match & pos_in_region)
    }

    if (any(in_excluded)) {
      n_removed <- sum(in_excluded) # nolint: object_usage_linter.
      exposure_iv <- exposure_iv[!in_excluded, ]
      if (verbose) {
        cli::cli_inform(
          "Removed {n_removed} instrument{?s} in excluded region{?s}."
        )
      }

      if (nrow(exposure_iv) == 0) {
        cli::cli_warn(
          "All instruments for {.val {exposure_id}} fall in excluded regions."
        )
        timing[["region_exclusion"]] <- proc.time()[["elapsed"]] - t0
        return(new_mr_result(
          status = "no_instruments",
          status_reason = paste0(
            "All instruments for '",
            exposure_id,
            "' fall in excluded regions"
          ),
          params = params,
          timing = timing
        ))
      }
    }
  }

  timing[["region_exclusion"]] <- proc.time()[["elapsed"]] - t0

  # --- Allele orientation check ---------------------------------------------

  # Runs on the full exposure/outcome (instruments plus a sample of shared
  # SNPs) BEFORE the outcome is narrowed to instrument rsIDs below, so that a
  # cis-MR with only a few instruments still gets a verdict. The instrument
  # harmonisation further down therefore passes check = FALSE.
  t0 <- proc.time()[["elapsed"]]

  check_allele_orientation_gwas(
    exposure,
    outcome,
    instrument_snps = exposure_iv$SNP,
    allele_check = allele_check,
    action = harmonise_action,
    verbose = verbose
  )

  timing[["allele_check"]] <- proc.time()[["elapsed"]] - t0

  # --- Format outcome and harmonise -----------------------------------------

  t0 <- proc.time()[["elapsed"]]

  # Pre-filter outcome to instrument SNPs before passing to TwoSampleMR.
  # format_gwas() returns the full GWAS file (potentially millions of rows).
  # TwoSampleMR::format_data() on a multi-million-row frame can overflow the
  # C stack; restricting to the ~10-50 instrument rsIDs avoids this and
  # speeds up harmonisation without affecting results.
  if ("rsids" %in% names(outcome) && nrow(exposure_iv) > 0) {
    outcome <- outcome[outcome$rsids %in% exposure_iv$SNP, , drop = FALSE]
  }

  if (nrow(outcome) == 0) {
    cli::cli_warn(
      "No overlapping SNPs between instruments and outcome for {.val {exposure_id}} -> {.val {outcome_id}}."
    )
    timing[["harmonisation"]] <- proc.time()[["elapsed"]] - t0
    return(new_mr_result(
      status = "no_harmonised_variants",
      status_reason = paste0(
        "No overlapping SNPs between instruments and outcome for '",
        exposure_id,
        "' -> '",
        outcome_id,
        "'"
      ),
      params = params,
      timing = timing
    ))
  }

  outcome_data <- TwoSampleMR::format_data(
    outcome,
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
  )

  harmonisation <- harmonise_and_filter(
    exposure_iv,
    outcome_data,
    allele_check = allele_check,
    action = harmonise_action,
    check = FALSE
  )
  harmonised <- harmonisation$data
  harmonisation <- harmonisation$raw

  timing[["harmonisation"]] <- proc.time()[["elapsed"]] - t0

  if (nrow(harmonised) == 0) {
    cli::cli_warn(
      "No variants remaining after harmonisation for {.val {exposure_id}}."
    )
    # The unfiltered frame goes back even here -- an empty result is exactly
    # when the caller needs to see which flag cost each variant its place.
    return(new_mr_result(
      harmonisation = harmonisation,
      status = "no_harmonised_variants",
      status_reason = paste0(
        "No variants remaining after harmonisation for '",
        exposure_id,
        "'"
      ),
      params = params,
      timing = timing
    ))
  }

  # --- Resolve sample size --------------------------------------------------

  exp_n <- resolve_sample_size(
    explicit_n = exposure_n,
    data_column = harmonised$samplesize.exposure,
    label = "exposure"
  )

  # --- LD correction --------------------------------------------------------

  t0 <- proc.time()[["elapsed"]]

  ld_mat <- NULL
  ld_input <- NULL
  if (ld_correct) {
    ld_mat <- compute_ld_matrix(
      snps = harmonised$SNP,
      bfile = bfile,
      plink_bin = plink_bin,
      plink_threads = plink_threads,
      plink_memory = plink_memory
    )
    aligned <- align_to_ld_matrix(harmonised, ld_mat)
    harmonised <- aligned$data
    ld_mat <- aligned$ld_matrix
    # One correlated MRInput serves every LD-corrected method below. The
    # matrix, not the `correl` flag, is what carries the correction:
    # mr_ivw()/mr_egger() read it from this object.
    ld_input <- MendelianRandomization::mr_input(
      bx = harmonised$beta.exposure,
      bxse = harmonised$se.exposure,
      by = harmonised$beta.outcome,
      byse = harmonised$se.outcome,
      correlation = ld_mat
    )
    # Every LD-corrected fit solves this matrix; a near-singular one makes
    # all of them unstable, and until now that happened silently.
    ld_rcond <- rcond(gls_weight_matrix(ld_mat, harmonised$se.outcome))
    if (ld_rcond < 1e-10) {
      cli::cli_warn(c(
        paste0(
          "The LD weight matrix for {.val {exposure_id}} is near-singular ",
          "(reciprocal condition number {signif(ld_rcond, 2)}); every ",
          "LD-corrected estimate is unstable."
        ),
        "i" = paste0(
          "Usual causes: identical or near-identical instruments (r2 ~ 1, ",
          "from absent clumping or a manual set with a duplicated variant), ",
          "or more instruments than reference-panel individuals, which makes ",
          "the sample correlation matrix singular."
        ),
        "i" = "Clump more stringently or drop the duplicate."
      ))
    }
  }

  timing[["ld_correction"]] <- proc.time()[["elapsed"]] - t0

  # Alignment can drop every instrument (none in the panel, or all ambiguous
  # palindromes); without this guard dispatch would fall into the multi-SNP
  # branch with zero rows and fail inside the GLS fit.
  if (nrow(harmonised) == 0) {
    cli::cli_warn(
      "No instruments remaining after LD alignment for {.val {exposure_id}}."
    )
    return(new_mr_result(
      harmonisation = harmonisation,
      ld_matrix = ld_mat,
      status = "no_harmonised_variants",
      status_reason = paste0(
        "No instruments remaining after LD alignment for '",
        exposure_id,
        "'"
      ),
      params = params,
      timing = timing
    ))
  }

  # --- F-statistics ---------------------------------------------------------

  f_per_snp <- (harmonised$beta.exposure / harmonised$se.exposure)^2
  f_stats <- list(
    per_snp = f_per_snp,
    mean = mean(f_per_snp),
    min = min(f_per_snp)
  )

  # --- Method dispatch ------------------------------------------------------

  n_snps <- nrow(harmonised)
  results_list <- list()
  methods_skipped <- character()

  # Skip reason when a method's registry minimum exceeds the instrument
  # count; NULL when it can run.
  too_few <- function(entry) {
    if (n_snps < entry$min_instruments) {
      paste0("Requires >= ", entry$min_instruments, " instruments")
    } else {
      NULL
    }
  }

  # Wald ratio for single instrument

  if (n_snps == 1) {
    t0 <- proc.time()[["elapsed"]]
    wald <- TwoSampleMR::mr(harmonised, method_list = "mr_wald_ratio")
    results_list[["wald"]] <- mr_result_row(
      mr_method_entry(label = "Wald ratio"),
      exposure_id,
      outcome_id,
      nsnp = wald$nsnp,
      b = wald$b,
      se = wald$se,
      pval = wald$pval,
      ld_corrected = FALSE
    )
    timing[["mr_wald_ratio"]] <- proc.time()[["elapsed"]] - t0

    # Skip every multi-SNP estimator (shortcuts + any raw TwoSampleMR names).
    # The diagnostics record their own instrument-count reasons below.
    multi <- registry$shortcut[
      !is.na(registry$shortcut) &
        registry$output == "$results" &
        registry$min_instruments > 1
    ]
    generic_tsm_methods <- methods[!methods %in% shortcut_methods]
    for (m in c(intersect(methods, multi), generic_tsm_methods)) {
      methods_skipped[m] <- "Only 1 instrument (Wald ratio used)"
    }
  } else {
    # IVW, multiplicative random effects
    if ("ivw_random" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("ivw_random")
      if (ld_correct) {
        # `model` is pinned: MendelianRandomization's "default" is fixed
        # effects below 4 instruments, which would silently change the
        # estimator this shortcut names (issue #27).
        fit <- MendelianRandomization::mr_ivw(
          ld_input,
          correl = TRUE,
          model = "random"
        )
        results_list[["ivw_random"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit@SNPs,
          b = fit@Estimate,
          se = fit@StdError,
          pval = fit@Pvalue,
          ld_corrected = TRUE
        )
      } else {
        fit <- TwoSampleMR::mr(harmonised, method_list = "mr_ivw")
        results_list[["ivw_random"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit$nsnp,
          b = fit$b,
          se = fit$se,
          pval = fit$pval,
          ld_corrected = FALSE
        )
      }
      timing[["mr_ivw_random"]] <- proc.time()[["elapsed"]] - t0
    }

    # IVW, fixed effects
    if ("ivw_fixed" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("ivw_fixed")
      if (ld_correct) {
        fit <- MendelianRandomization::mr_ivw(
          ld_input,
          correl = TRUE,
          model = "fixed"
        )
        results_list[["ivw_fixed"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit@SNPs,
          b = fit@Estimate,
          se = fit@StdError,
          pval = fit@Pvalue,
          ld_corrected = TRUE
        )
      } else {
        fit <- TwoSampleMR::mr(harmonised, method_list = "mr_ivw_fe")
        results_list[["ivw_fixed"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit$nsnp,
          b = fit$b,
          se = fit$se,
          pval = fit$pval,
          ld_corrected = FALSE
        )
      }
      timing[["mr_ivw_fixed"]] <- proc.time()[["elapsed"]] - t0
    }

    # Egger (requires >= 3 SNPs; always multiplicative random effects on
    # both paths -- mr_egger() has no model argument)
    if ("egger" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("egger")
      reason <- too_few(entry)
      if (!is.null(reason)) {
        methods_skipped["egger"] <- reason
      } else if (ld_correct) {
        fit <- MendelianRandomization::mr_egger(ld_input, correl = TRUE)
        results_list[["egger"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit@SNPs,
          b = fit@Estimate,
          se = fit@StdError.Est,
          pval = fit@Pvalue.Est,
          ld_corrected = TRUE
        )
      } else {
        fit <- TwoSampleMR::mr(harmonised, method_list = "mr_egger_regression")
        results_list[["egger"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit$nsnp,
          b = fit$b,
          se = fit$se,
          pval = fit$pval,
          ld_corrected = FALSE
        )
      }
      timing[["mr_egger"]] <- proc.time()[["elapsed"]] - t0
    }

    # Weighted median (requires >= 3 SNPs; no LD-corrected form)
    if ("weighted_median" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("weighted_median")
      reason <- too_few(entry)
      if (!is.null(reason)) {
        methods_skipped["weighted_median"] <- reason
      } else {
        if (ld_correct) {
          warn_no_ld_correction("weighted_median")
        }
        fit <- TwoSampleMR::mr(harmonised, method_list = "mr_weighted_median")
        results_list[["weighted_median"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = fit$nsnp,
          b = fit$b,
          se = fit$se,
          pval = fit$pval,
          ld_corrected = FALSE
        )
      }
      timing[["mr_weighted_median"]] <- proc.time()[["elapsed"]] - t0
    }

    # MR-PRESSO (requires >= 3 SNPs; no LD-corrected form)
    if ("presso" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("presso")
      reason <- too_few(entry)
      if (!is.null(reason)) {
        methods_skipped["presso"] <- reason
      } else {
        if (ld_correct) {
          warn_no_ld_correction("presso")
        }
        presso_result <- tryCatch(
          {
            TwoSampleMR::run_mr_presso(
              harmonised,
              NbDistribution = presso_n_dist
            )
          },
          error = function(e) {
            cli::cli_warn("MR-PRESSO failed: {conditionMessage(e)}")
            NULL
          }
        )

        if (!is.null(presso_result) && length(presso_result) > 0) {
          # MR-PRESSO returns a list; extract main result
          presso_main <- presso_result[[1]]$`Main MR results`
          # Use "Raw" estimate (row 1)
          if (!is.null(presso_main) && nrow(presso_main) > 0) {
            results_list[["presso"]] <- mr_result_row(
              entry,
              exposure_id,
              outcome_id,
              nsnp = n_snps,
              b = presso_main$`Causal Estimate`[1],
              se = presso_main$Sd[1],
              pval = presso_main$`P-value`[1],
              ld_corrected = FALSE
            )
          }
        }
      }
      timing[["mr_presso"]] <- proc.time()[["elapsed"]] - t0
    }

    # ConMix (no LD-corrected form: mr_conmix() takes no correlation matrix)
    if ("conmix" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      entry <- mr_method_entry("conmix")
      if (ld_correct) {
        warn_no_ld_correction("conmix")
      }
      conmix_result <- tryCatch(
        {
          mr_input <- MendelianRandomization::mr_input(
            bx = harmonised$beta.exposure,
            bxse = harmonised$se.exposure,
            by = harmonised$beta.outcome,
            byse = harmonised$se.outcome
          )
          MendelianRandomization::mr_conmix(mr_input)
        },
        error = function(e) {
          cli::cli_warn("ConMix failed: {conditionMessage(e)}")
          NULL
        }
      )

      if (!is.null(conmix_result)) {
        results_list[["conmix"]] <- mr_result_row(
          entry,
          exposure_id,
          outcome_id,
          nsnp = n_snps,
          b = conmix_result@Estimate,
          se = NA_real_,
          pval = conmix_result@Pvalue,
          ld_corrected = FALSE
        )
      }
      timing[["mr_conmix"]] <- proc.time()[["elapsed"]] - t0
    }

    # Raw TwoSampleMR methods (names with no shortcut), under TwoSampleMR's
    # own label. Never LD-corrected.
    passthrough <- mr_method_entry(engine = "TwoSampleMR::mr")
    generic_tsm <- methods[!methods %in% shortcut_methods]
    for (m in generic_tsm) {
      t0 <- proc.time()[["elapsed"]]
      if (ld_correct) {
        warn_no_ld_correction(m)
      }
      generic_res <- tryCatch(
        {
          TwoSampleMR::mr(harmonised, method_list = m)
        },
        error = function(e) {
          cli::cli_warn("{.val {m}} failed: {conditionMessage(e)}")
          methods_skipped[m] <<- paste("Failed:", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(generic_res) && nrow(generic_res) > 0) {
        results_list[[m]] <- mr_result_row(
          passthrough,
          exposure_id,
          outcome_id,
          nsnp = generic_res$nsnp,
          b = generic_res$b,
          se = generic_res$se,
          pval = generic_res$pval,
          ld_corrected = FALSE,
          label = generic_res$method
        )
      }
      timing[[paste0("mr_generic_", m)]] <- proc.time()[["elapsed"]] - t0
    }
  }

  # Steiger (works with any number of SNPs if sample sizes available)
  steiger_result <- NULL
  if ("steiger" %in% methods) {
    t0 <- proc.time()[["elapsed"]]
    if (is.null(exp_n)) {
      methods_skipped["steiger"] <- "Exposure sample size not available"
    } else {
      steiger_result <- tryCatch(
        {
          TwoSampleMR::steiger_filtering(harmonised)
        },
        error = function(e) {
          cli::cli_warn("Steiger filtering failed: {conditionMessage(e)}")
          methods_skipped["steiger"] <<- paste("Failed:", conditionMessage(e))
          NULL
        }
      )
    }
    timing[["mr_steiger"]] <- proc.time()[["elapsed"]] - t0
  }

  # Pleiotropy test -- Egger intercept (requires >= 3 SNPs).
  # Runs automatically when "egger" is in methods; the "pleiotropy" shortcut
  # also triggers it independently (e.g. without Egger). Under ld_correct the
  # intercept comes from the correlated Egger fit (issue #31).
  pleiotropy_result <- NULL
  if ("pleiotropy" %in% methods || "egger" %in% methods) {
    t0 <- proc.time()[["elapsed"]]
    reason <- too_few(mr_method_entry("pleiotropy"))
    if (!is.null(reason)) {
      methods_skipped["pleiotropy"] <- reason
    } else {
      pleiotropy_result <- tryCatch(
        {
          if (ld_correct) {
            pleiotropy_correlated(ld_input, harmonised)
          } else {
            cbind(
              TwoSampleMR::mr_pleiotropy_test(harmonised),
              ld_corrected = FALSE
            )
          }
        },
        error = function(e) {
          cli::cli_warn(
            "Pleiotropy test failed: {conditionMessage(e)}"
          )
          methods_skipped["pleiotropy"] <<- paste(
            "Failed:",
            conditionMessage(e)
          )
          NULL
        }
      )
    }
    timing[["mr_pleiotropy"]] <- proc.time()[["elapsed"]] - t0
  }

  # Heterogeneity test -- Cochran's Q (requires >= 2 SNPs). Under ld_correct
  # this is the generalised Q from the correlated fits (issue #31).
  heterogeneity_result <- NULL
  if ("heterogeneity" %in% methods) {
    t0 <- proc.time()[["elapsed"]]
    reason <- too_few(mr_method_entry("heterogeneity"))
    if (!is.null(reason)) {
      methods_skipped["heterogeneity"] <- reason
    } else {
      heterogeneity_result <- tryCatch(
        {
          if (ld_correct) {
            heterogeneity_correlated(ld_input, harmonised)
          } else {
            cbind(
              TwoSampleMR::mr_heterogeneity(harmonised),
              ld_corrected = FALSE
            )
          }
        },
        error = function(e) {
          cli::cli_warn(
            "Heterogeneity test failed: {conditionMessage(e)}"
          )
          methods_skipped["heterogeneity"] <<- paste(
            "Failed:",
            conditionMessage(e)
          )
          NULL
        }
      )
    }
    timing[["mr_heterogeneity"]] <- proc.time()[["elapsed"]] - t0
  }

  # Leave-one-out analysis (requires >= 3 SNPs -- with 2, dropping one just
  # reproduces the remaining SNP's Wald ratio, which isn't informative).
  # Under ld_correct each refit uses the reduced weight matrix (issue #31).
  loo_result <- NULL
  if ("loo" %in% methods) {
    t0 <- proc.time()[["elapsed"]]
    reason <- too_few(mr_method_entry("loo"))
    if (!is.null(reason)) {
      methods_skipped["loo"] <- reason
    } else {
      loo_result <- tryCatch(
        {
          if (ld_correct) {
            loo_correlated(harmonised, ld_mat)
          } else {
            cbind(
              TwoSampleMR::mr_leaveoneout(harmonised),
              ld_corrected = FALSE
            )
          }
        },
        error = function(e) {
          cli::cli_warn(
            "Leave-one-out analysis failed: {conditionMessage(e)}"
          )
          methods_skipped["loo"] <<- paste(
            "Failed:",
            conditionMessage(e)
          )
          NULL
        }
      )
    }
    timing[["mr_loo"]] <- proc.time()[["elapsed"]] - t0
  }

  # --- Assemble results -----------------------------------------------------

  if (length(results_list) == 0) {
    results_df <- TwoSampleMR::generate_odds_ratios(empty_mr_results())
  } else {
    # Base rbind() on purpose, not dplyr::bind_rows(): every row comes from
    # mr_result_row() with an identical schema, and rbind() errors the moment
    # a branch breaks that, where bind_rows() would silently fill NA.
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    results_df <- TwoSampleMR::generate_odds_ratios(results_df)
  }

  new_mr_result(
    results = results_df,
    instruments = harmonised,
    harmonisation = harmonisation,
    f_stats = f_stats,
    steiger = steiger_result,
    pleiotropy = pleiotropy_result,
    heterogeneity = heterogeneity_result,
    loo = loo_result,
    methods_skipped = methods_skipped,
    ld_matrix = ld_mat,
    params = params,
    timing = timing
  )
}

#' Validate exclude_regions argument
#'
#' @param exclude_regions Data frame to validate.
#'
#' @return Invisibly returns `TRUE` if valid; otherwise aborts with an error.
#'
#' @keywords internal
validate_exclude_regions <- function(exclude_regions) {
  if (!is.data.frame(exclude_regions)) {
    cli::cli_abort("{.arg exclude_regions} must be a data frame.")
  }

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(exclude_regions))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg exclude_regions} must have columns {.val {required_cols}}; missing {.val {missing_cols}}."
    )
  }

  if (
    !rlang::is_integerish(exclude_regions$start) ||
      !rlang::is_integerish(exclude_regions$end)
  ) {
    cli::cli_abort(
      "{.arg exclude_regions} columns {.val start} and {.val end} must be whole numbers."
    )
  }

  if (any(exclude_regions$start < 0) || any(exclude_regions$end < 0)) {
    cli::cli_abort(
      "{.arg exclude_regions} columns {.val start} and {.val end} must be positive."
    )
  }

  if (any(exclude_regions$start > exclude_regions$end)) {
    cli::cli_abort(
      "{.arg exclude_regions}: {.val start} must be <= {.val end} for all rows."
    )
  }

  invisible(TRUE)
}
