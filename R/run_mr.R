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
#' - 1 SNP: Wald ratio only
#' - 2 SNPs: IVW (+ ConMix/Steiger if requested); Egger/weighted_median/PRESSO
#'   skipped
#' - 3+ SNPs: all methods in `methods` are attempted
#'
#' When `ld_correct = TRUE`, IVW and Egger use
#' `MendelianRandomization::mr_ivw()` and `MendelianRandomization::mr_egger()`
#' with `correl = TRUE`.
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
#' @param mhc_remove Logical. Remove instruments in the MHC region
#'   (chr6:26-34Mb). Default `FALSE`.
#' @param methods Character vector of MR methods to run. Options: `"ivw"`,
#'   `"egger"`, `"weighted_median"`, `"presso"`, `"conmix"`, `"steiger"`.
#' @param ld_correct Logical. Use LD-corrected IVW/Egger via the
#'   `MendelianRandomization` package. Requires `bfile`. Default `FALSE`.
#' @param exposure_n Numeric. Exposure sample size. If `NULL`, inferred from
#'   `samplesize.exposure` column.
#' @param presso_n_dist Integer. Number of distributions for MR-PRESSO. Default
#'   `1000`.
#'
#' @return An `mr_result` object, or `NULL` if no instruments survive filtering.
#'
#' @export
run_mr <- function(exposure,
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
                   mhc_remove = FALSE,
                   methods = c("ivw", "egger", "weighted_median",
                               "presso", "conmix", "steiger"),
                   ld_correct = FALSE,
                   exposure_n = NULL,
                   presso_n_dist = 1000) {
  # --- Validate arguments ---------------------------------------------------

  if (ld_correct && is.null(bfile)) {
    cli::cli_abort("{.arg bfile} is required when {.code ld_correct = TRUE}.")
  }

  methods <- match.arg(
    methods,
    choices = c("ivw", "egger", "weighted_median", "presso", "conmix",
                "steiger"),
    several.ok = TRUE
  )

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
    mhc_remove = mhc_remove,
    methods = methods,
    ld_correct = ld_correct,
    exposure_n = exposure_n,
    presso_n_dist = presso_n_dist
  )

  # --- Instrument selection -------------------------------------------------

  if (!is.null(instruments)) {
    # Manual mode
    cli::cli_inform("Using {length(instruments)} manual instrument{?s}.")
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
      return(NULL)
    }
  } else if (!is.null(instrument_region)) {
    # Cis-MR mode
    cli::cli_inform("Cis-MR mode: chr{instrument_region$chromosome}:{instrument_region$start}-{instrument_region$end} (+/- {window}bp).")

    exposure_iv <- exposure |>
      dplyr::filter(
        as.character(.data$chr.exposure) == as.character(instrument_region$chromosome),
        .data$pos.exposure >= (instrument_region$start - window),
        .data$pos.exposure <= (instrument_region$end + window)
      ) |>
      dplyr::filter(.data$pval.exposure < pval_thresh)

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn("No significant instruments in cis region for {.val {exposure_id}}.")
      return(NULL)
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
      pop = pop
    )
    exposure_iv <- exposure_iv[exposure_iv$SNP %in% clumped$rsid, ]

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn("No instruments remaining after clumping for {.val {exposure_id}}.")
      return(NULL)
    }
  } else {
    # Genome-wide mode
    cli::cli_inform("Genome-wide mode: selecting instruments at p < {pval_thresh}.")

    exposure_iv <- exposure |>
      dplyr::filter(.data$pval.exposure < pval_thresh)

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn("No genome-wide significant instruments for {.val {exposure_id}}.")
      return(NULL)
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
      pop = pop
    )
    exposure_iv <- exposure_iv[exposure_iv$SNP %in% clumped$rsid, ]

    if (nrow(exposure_iv) == 0) {
      cli::cli_warn("No instruments remaining after clumping for {.val {exposure_id}}.")
      return(NULL)
    }
  }

  # --- MHC removal ----------------------------------------------------------

  if (mhc_remove && "chr.exposure" %in% colnames(exposure_iv)) {
    in_mhc <- in_mhc_region(
      chr = exposure_iv$chr.exposure,
      start = exposure_iv$pos.exposure,
      end = exposure_iv$pos.exposure
    )

    if (any(in_mhc)) {
      n_removed <- sum(in_mhc)
      exposure_iv <- exposure_iv[!in_mhc, ]
      cli::cli_inform("Removed {n_removed} instrument{?s} in MHC region.")

      if (nrow(exposure_iv) == 0) {
        cli::cli_warn("All instruments for {.val {exposure_id}} are in the MHC region.")
        return(NULL)
      }
    }
  }

  # --- Format outcome and harmonise -----------------------------------------

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

  if (nrow(harmonised) == 0) {
    cli::cli_warn("No variants remaining after harmonisation for {.val {exposure_id}}.")
    return(NULL)
  }

  # --- Resolve sample size --------------------------------------------------

  exp_n <- resolve_sample_size(
    explicit_n = exposure_n,
    data_column = harmonised$samplesize.exposure,
    label = "exposure"
  )

  # --- LD correction --------------------------------------------------------

  ld_mat <- NULL
  if (ld_correct) {
    ld_mat <- compute_ld_matrix(
      snps = harmonised$SNP,
      bfile = bfile,
      plink_bin = plink_bin
    )
    aligned <- align_to_ld_matrix(harmonised, ld_mat)
    harmonised <- aligned$data
    ld_mat <- aligned$ld_matrix
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

  # Wald ratio for single instrument

  if (n_snps == 1) {
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

    # Skip all multi-SNP methods
    multi_methods <- intersect(
      methods,
      c("ivw", "egger", "weighted_median", "presso", "conmix")
    )
    for (m in multi_methods) {
      methods_skipped[m] <- "Only 1 instrument (Wald ratio used)"
    }
  } else {
    # IVW
    if ("ivw" %in% methods) {
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
    }

    # Egger (requires >= 3 SNPs)
    if ("egger" %in% methods) {
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
        egger_mr <- TwoSampleMR::mr(harmonised, method_list = "mr_egger_regression")
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
    }

    # Weighted median (requires >= 3 SNPs)
    if ("weighted_median" %in% methods) {
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
    }

    # MR-PRESSO (requires >= 3 SNPs)
    if ("presso" %in% methods) {
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
    }

    # ConMix
    if ("conmix" %in% methods) {
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
    }
  }

  # Steiger (works with any number of SNPs if sample sizes available)
  steiger_result <- NULL
  if ("steiger" %in% methods) {
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
  } else {
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
  }

  new_mr_result(
    results = results_df,
    instruments = harmonised,
    f_stats = f_stats,
    steiger = steiger_result,
    methods_skipped = methods_skipped,
    ld_matrix = ld_mat,
    params = params
  )
}
