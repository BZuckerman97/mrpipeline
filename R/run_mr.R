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
#' Methods are dispatched based on the number of instruments after clumping:
#' - 1 SNP: Wald ratio only; all other methods skipped
#' - 2+ SNPs: IVW, IVW-FE, ConMix, Steiger, and any raw TwoSampleMR methods
#'   are attempted; Egger, weighted median, and PRESSO require >= 3 SNPs
#' - 3+ SNPs: all methods in `methods` are attempted
#'
#' Generic TwoSampleMR methods (raw `mr_*` names) are dispatched via
#' `TwoSampleMR::mr()` and errors are caught and reported as skipped.
#'
#' When `ld_correct = TRUE`, IVW and Egger use
#' `MendelianRandomization::mr_ivw()` and `MendelianRandomization::mr_egger()`
#' with `correl = TRUE`.
#'
#' When `"egger"` is in `methods` and there are >= 3 instruments,
#' `TwoSampleMR::mr_pleiotropy_test()` (the Egger intercept test) is always
#' run automatically and its result stored in `$pleiotropy`. You do not need
#' to add `"pleiotropy"` to `methods` separately. The `"pleiotropy"` shortcut
#' remains available for running the intercept test without Egger.
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
#' @param methods Character vector of MR methods to run. Named shortcuts:
#'   `"ivw"` (IVW random effects), `"ivw_fe"` (IVW fixed effects),
#'   `"egger"` (MR Egger), `"weighted_median"` (weighted median),
#'   `"presso"` (MR-PRESSO), `"conmix"` (ContMix), `"steiger"` (Steiger
#'   filtering), `"pleiotropy"` (Egger intercept test; result stored in
#'   `$pleiotropy`, not `$results`). You may also pass any method name from
#'   `TwoSampleMR::mr_method_list()$obj` directly (e.g. `"mr_simple_median"`,
#'   `"mr_raps"`). Note: `"ivw_fe"` does not support `ld_correct = TRUE`.
#' @param ld_correct Logical. Use LD-corrected IVW/Egger via the
#'   `MendelianRandomization` package. Requires `bfile`. Default `FALSE`.
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
#' @param verbose Logical. If `TRUE`, emit informational messages via
#'   [cli::cli_inform()]. Warnings and errors are always emitted regardless.
#'   Default `TRUE`.
#'
#' @return An `mr_result` object. Check `result$status` for `"success"` vs
#'   failure reasons. The `$results` data frame includes `or`, `or_lci95`, and
#'   `or_uci95` columns (from [TwoSampleMR::generate_odds_ratios()]) alongside
#'   the raw `b` and `se`. The `$timing` field contains a named numeric vector
#'   of elapsed seconds for each major step.
#'
#' @examples
#' \dontrun{
#' # Cis-MR using bundled CD40/Sjogren's data
#' bfile <- system.file("extdata", "ld_ref", package = "mrpipeline")
#' result <- run_mr(
#'   exposure = cd40_exposure,
#'   exposure_id = "CD40",
#'   outcome = sjogren_outcome,
#'   outcome_id = "SjD",
#'   instrument_region = list(chromosome = "20", start = 44746911, end = 44758502),
#'   bfile = bfile,
#'   methods = c("ivw", "egger", "weighted_median")
#' )
#' result
#' summary(result)
#' }
#'
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
  methods = c("ivw", "egger", "weighted_median", "presso", "conmix", "steiger"),
  ld_correct = FALSE,
  exposure_n = NULL,
  presso_n_dist = 1000,
  plink_threads = plink_option("threads"),
  plink_memory = plink_option("memory"),
  verbose = TRUE
) {
  # --- Validate arguments ---------------------------------------------------

  if (ld_correct && is.null(bfile)) {
    cli::cli_abort("{.arg bfile} is required when {.code ld_correct = TRUE}.")
  }

  shortcut_methods <- c(
    "ivw",
    "ivw_fe",
    "egger",
    "weighted_median",
    "presso",
    "conmix",
    "steiger",
    "pleiotropy"
  )
  tsm_available <- setdiff(TwoSampleMR::mr_method_list()$obj, "mr_wald_ratio")
  unknown_methods <- setdiff(methods, c(shortcut_methods, tsm_available))
  if (length(unknown_methods) > 0) {
    cli::cli_abort(
      c(
        "Unknown method{?s}: {.val {unknown_methods}}.",
        "i" = "Named shortcuts: {.val {shortcut_methods}}.",
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
    presso_n_dist = presso_n_dist
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
      cli::cli_inform(
        "Cis-MR mode: chr{instrument_region$chromosome}:{instrument_region$start}-{instrument_region$end} (+/- {window}bp)."
      )
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
      in_excluded <- in_excluded |
        (as.character(exposure_iv$chr.exposure) ==
          as.character(exclude_regions$chr[i]) &
          exposure_iv$pos.exposure >= exclude_regions$start[i] &
          exposure_iv$pos.exposure <= exclude_regions$end[i])
    }

    if (any(in_excluded)) {
      n_removed <- sum(in_excluded)
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

  harmonised <- harmonise_and_filter(exposure_iv, outcome_data)

  timing[["harmonisation"]] <- proc.time()[["elapsed"]] - t0

  if (nrow(harmonised) == 0) {
    cli::cli_warn(
      "No variants remaining after harmonisation for {.val {exposure_id}}."
    )
    return(new_mr_result(
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
  }

  timing[["ld_correction"]] <- proc.time()[["elapsed"]] - t0

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

  # Wald ratio for single instrument

  if (n_snps == 1) {
    t0 <- proc.time()[["elapsed"]]
    wald <- TwoSampleMR::mr(harmonised, method_list = "mr_wald_ratio")
    results_list[["Wald ratio"]] <- data.frame(
      exposure = exposure_id,
      outcome = outcome_id,
      method = "Wald ratio",
      nsnp = wald$nsnp,
      b = wald$b,
      se = wald$se,
      pval = wald$pval,
      stringsAsFactors = FALSE
    )
    timing[["mr_ivw"]] <- proc.time()[["elapsed"]] - t0

    # Skip all multi-SNP methods (shortcuts + any raw TwoSampleMR names)
    multi_shortcut_names <- c(
      "ivw",
      "ivw_fe",
      "egger",
      "weighted_median",
      "presso",
      "conmix"
    )
    generic_tsm_methods <- setdiff(methods, c(multi_shortcut_names, "steiger"))
    all_multi <- c(
      intersect(methods, multi_shortcut_names),
      generic_tsm_methods
    )
    for (m in all_multi) {
      methods_skipped[m] <- "Only 1 instrument (Wald ratio used)"
    }
  } else {
    # IVW
    if ("ivw" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      if (ld_correct) {
        mr_input <- MendelianRandomization::mr_input(
          bx = harmonised$beta.exposure,
          bxse = harmonised$se.exposure,
          by = harmonised$beta.outcome,
          byse = harmonised$se.outcome,
          correlation = ld_mat
        )
        ivw_res <- MendelianRandomization::mr_ivw(mr_input, correl = TRUE)
        results_list[["IVW (LD-corrected)"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "IVW (LD-corrected)",
          nsnp = ivw_res@SNPs,
          b = ivw_res@Estimate,
          se = ivw_res@StdError,
          pval = ivw_res@Pvalue,
          stringsAsFactors = FALSE
        )
      } else {
        ivw_mr <- TwoSampleMR::mr(harmonised, method_list = "mr_ivw")
        results_list[["IVW"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "Inverse variance weighted",
          nsnp = ivw_mr$nsnp,
          b = ivw_mr$b,
          se = ivw_mr$se,
          pval = ivw_mr$pval,
          stringsAsFactors = FALSE
        )
      }
      timing[["mr_ivw"]] <- proc.time()[["elapsed"]] - t0
    }

    # IVW fixed effects
    if ("ivw_fe" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      if (ld_correct) {
        cli::cli_warn(
          "{.val 'ivw_fe'} does not support ld_correct; running uncorrected."
        )
      }
      ivw_fe_mr <- TwoSampleMR::mr(harmonised, method_list = "mr_ivw_fe")
      results_list[["IVW (fixed effects)"]] <- data.frame(
        exposure = exposure_id,
        outcome = outcome_id,
        method = "IVW (fixed effects)",
        nsnp = ivw_fe_mr$nsnp,
        b = ivw_fe_mr$b,
        se = ivw_fe_mr$se,
        pval = ivw_fe_mr$pval,
        stringsAsFactors = FALSE
      )
      timing[["mr_ivw_fe"]] <- proc.time()[["elapsed"]] - t0
    }

    # Egger (requires >= 3 SNPs)
    if ("egger" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      if (n_snps < 3) {
        methods_skipped["egger"] <- "Requires >= 3 instruments"
      } else if (ld_correct) {
        mr_input <- MendelianRandomization::mr_input(
          bx = harmonised$beta.exposure,
          bxse = harmonised$se.exposure,
          by = harmonised$beta.outcome,
          byse = harmonised$se.outcome,
          correlation = ld_mat
        )
        egger_res <- MendelianRandomization::mr_egger(mr_input, correl = TRUE)
        results_list[["Egger (LD-corrected)"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "MR Egger (LD-corrected)",
          nsnp = egger_res@SNPs,
          b = egger_res@Estimate,
          se = egger_res@StdError.Est,
          pval = egger_res@Pvalue.Est,
          stringsAsFactors = FALSE
        )
      } else {
        egger_mr <- TwoSampleMR::mr(
          harmonised,
          method_list = "mr_egger_regression"
        )
        results_list[["Egger"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "MR Egger",
          nsnp = egger_mr$nsnp,
          b = egger_mr$b,
          se = egger_mr$se,
          pval = egger_mr$pval,
          stringsAsFactors = FALSE
        )
      }
      timing[["mr_egger"]] <- proc.time()[["elapsed"]] - t0
    }

    # Weighted median (requires >= 3 SNPs)
    if ("weighted_median" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      if (n_snps < 3) {
        methods_skipped["weighted_median"] <- "Requires >= 3 instruments"
      } else {
        wm_mr <- TwoSampleMR::mr(harmonised, method_list = "mr_weighted_median")
        results_list[["Weighted median"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "Weighted median",
          nsnp = wm_mr$nsnp,
          b = wm_mr$b,
          se = wm_mr$se,
          pval = wm_mr$pval,
          stringsAsFactors = FALSE
        )
      }
      timing[["mr_weighted_median"]] <- proc.time()[["elapsed"]] - t0
    }

    # MR-PRESSO (requires >= 3 SNPs)
    if ("presso" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
      if (n_snps < 3) {
        methods_skipped["presso"] <- "Requires >= 3 instruments"
      } else {
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
            results_list[["PRESSO"]] <- data.frame(
              exposure = exposure_id,
              outcome = outcome_id,
              method = "MR-PRESSO",
              nsnp = n_snps,
              b = presso_main$`Causal Estimate`[1],
              se = presso_main$Sd[1],
              pval = presso_main$`P-value`[1],
              stringsAsFactors = FALSE
            )
          }
        }
      }
      timing[["mr_presso"]] <- proc.time()[["elapsed"]] - t0
    }

    # ConMix
    if ("conmix" %in% methods) {
      t0 <- proc.time()[["elapsed"]]
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
        results_list[["ConMix"]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = "ConMix",
          nsnp = n_snps,
          b = conmix_result@Estimate,
          se = NA_real_,
          pval = conmix_result@Pvalue,
          stringsAsFactors = FALSE
        )
      }
      timing[["mr_conmix"]] <- proc.time()[["elapsed"]] - t0
    }

    # Generic TwoSampleMR methods (raw mr_* names not handled by shortcuts)
    shortcut_names <- c(
      "ivw",
      "ivw_fe",
      "egger",
      "weighted_median",
      "presso",
      "conmix",
      "steiger"
    )
    generic_tsm <- methods[!methods %in% shortcut_names]
    for (m in generic_tsm) {
      t0 <- proc.time()[["elapsed"]]
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
        results_list[[m]] <- data.frame(
          exposure = exposure_id,
          outcome = outcome_id,
          method = generic_res$method,
          nsnp = generic_res$nsnp,
          b = generic_res$b,
          se = generic_res$se,
          pval = generic_res$pval,
          stringsAsFactors = FALSE
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
  # also triggers it independently (e.g. without Egger).
  pleiotropy_result <- NULL
  if ("pleiotropy" %in% methods || "egger" %in% methods) {
    t0 <- proc.time()[["elapsed"]]
    if (n_snps < 3) {
      methods_skipped["pleiotropy"] <- "Requires >= 3 instruments"
    } else {
      pleiotropy_result <- tryCatch(
        {
          TwoSampleMR::mr_pleiotropy_test(harmonised)
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

  # --- Assemble results -----------------------------------------------------

  if (length(results_list) == 0) {
    results_df <- data.frame(
      exposure = character(),
      outcome = character(),
      method = character(),
      nsnp = integer(),
      b = numeric(),
      se = numeric(),
      pval = numeric(),
      stringsAsFactors = FALSE
    )
    results_df <- TwoSampleMR::generate_odds_ratios(results_df)
  } else {
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    results_df <- TwoSampleMR::generate_odds_ratios(results_df)
  }

  new_mr_result(
    results = results_df,
    instruments = harmonised,
    f_stats = f_stats,
    steiger = steiger_result,
    pleiotropy = pleiotropy_result,
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
