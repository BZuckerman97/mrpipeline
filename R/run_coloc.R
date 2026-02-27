#' Perform colocalization analysis
#'
#' Runs colocalization between a pre-formatted exposure and a raw outcome
#' dataset within a specified gene region. Supports coloc.abf, SuSiE-based
#' coloc (coloc.susie), coloc.signals, and colocPropTest. Returns a
#' `coloc_result` S3 object.
#'
#' @section Methods:
#' - `"abf"` — Approximate Bayes Factor colocalization via
#'   [coloc::coloc.abf()]
#' - `"susie"` — SuSiE fine-mapping + colocalization via
#'   [coloc::runsusie()] and [coloc::coloc.susie()]
#' - `"signals"` — Multi-signal colocalization via
#'   [coloc::coloc.signals()]. Requires `"susie"` to have run successfully.
#' - `"prop_test"` — Proportionality test via
#'   `colocPropTest::coloc.prop.test()`. Requires `"signals"` to have run
#'   successfully and the `colocPropTest` package to be installed.
#'
#' @param exposure Data frame of formatted exposure data (output of
#'   [TwoSampleMR::format_data()] or `format_pqtl_*()` functions).
#' @param outcome Data frame of outcome summary statistics with standardised
#'   columns: `rsids`, `chr`, `pos`, `beta`, `se`, `eaf`, `pval`, `n`,
#'   `effect_allele`, `other_allele`.
#' @param gene_chr Character or integer. Chromosome of the gene/region.
#' @param gene_start Integer. Start position (bp) of the gene/region.
#' @param gene_end Integer. End position (bp) of the gene/region.
#' @param coloc_window Integer. Window in base pairs to extend around the
#'   gene region. Default `10000L` (10 kb).
#' @param exposure_n Integer. Exposure sample size. If `NULL`, inferred from
#'   `samplesize.exposure` column.
#' @param outcome_n Integer. Outcome sample size. If `NULL`, inferred from
#'   `n` column of `outcome`.
#' @param exposure_type Character. Type of exposure trait: `"quant"` or
#'   `"cc"`. Default `"quant"`.
#' @param outcome_type Character. Type of outcome trait: `"quant"` or
#'   `"cc"`. Default `"quant"`.
#' @param exposure_s Numeric. Proportion of cases for case-control exposure.
#'   Required when `exposure_type = "cc"`.
#' @param outcome_s Numeric. Proportion of cases for case-control outcome.
#'   Required when `outcome_type = "cc"`.
#' @param exposure_sdY Numeric. Standard deviation of the exposure trait
#'   (for quantitative traits). Default `1`.
#' @param outcome_sdY Numeric. Standard deviation of the outcome trait
#'   (for quantitative traits). Default `1`.
#' @param bfile Character. Path to PLINK bfile prefix (without
#'   .bed/.bim/.fam) for LD reference. Required.
#' @param plink_bin Character. Path to PLINK binary. Auto-detected if `NULL`.
#' @param methods Character vector of colocalization methods to run.
#'   Options: `"abf"`, `"susie"`, `"signals"`, `"prop_test"`. Default
#'   `c("abf", "susie", "signals")`.
#' @param p1 Numeric. Prior probability a SNP is associated with trait 1.
#'   Default `1e-4`.
#' @param p2 Numeric. Prior probability a SNP is associated with trait 2.
#'   Default `1e-4`.
#' @param p12 Numeric. Prior probability a SNP is associated with both
#'   traits. Default `1e-5`.
#'
#' @return A `coloc_result` object. Check `result$status` for `"success"` vs
#'   failure reasons.
#'
#' @seealso [new_coloc_result()] for the S3 class structure,
#'   [print.coloc_result()] and [summary.coloc_result()] for display methods.
#'
#' @family colocalization
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_coloc(
#'   exposure = formatted_exposure,
#'   outcome = outcome_data,
#'   gene_chr = 20, gene_start = 44746911, gene_end = 44758502,
#'   bfile = "/path/to/ld_reference",
#'   methods = "abf"
#' )
#' print(result)
#' summary(result)
#' }
run_coloc <- function(
  exposure,
  outcome,
  gene_chr,
  gene_start,
  gene_end,
  coloc_window = 10000L,
  exposure_n = NULL,
  outcome_n = NULL,
  exposure_type = c("quant", "cc"),
  outcome_type = c("quant", "cc"),
  exposure_s = NULL,
  outcome_s = NULL,
  exposure_sdY = 1,
  outcome_sdY = 1,
  bfile,
  plink_bin = NULL,
  methods = c("abf", "susie", "signals"),
  p1 = 1e-4,

  p2 = 1e-4,
  p12 = 1e-5
) {
  # --- Validate arguments ---------------------------------------------------

  if (missing(bfile)) {
    cli::cli_abort("{.arg bfile} is required for colocalization analysis.")
  }

  exposure_type <- rlang::arg_match(exposure_type)
  outcome_type <- rlang::arg_match(outcome_type)
  methods <- match.arg(
    methods,
    choices = c("abf", "susie", "signals", "prop_test"),
    several.ok = TRUE
  )

  if (exposure_type == "cc" && is.null(exposure_s)) {
    cli::cli_abort(
      "{.arg exposure_s} is required when {.code exposure_type = \"cc\"}."
    )
  }
  if (outcome_type == "cc" && is.null(outcome_s)) {
    cli::cli_abort(
      "{.arg outcome_s} is required when {.code outcome_type = \"cc\"}."
    )
  }

  params <- list(
    gene_chr = gene_chr,
    gene_start = gene_start,
    gene_end = gene_end,
    coloc_window = coloc_window,
    exposure_n = exposure_n,
    outcome_n = outcome_n,
    exposure_type = exposure_type,
    outcome_type = outcome_type,
    exposure_s = exposure_s,
    outcome_s = outcome_s,
    exposure_sdY = exposure_sdY,
    outcome_sdY = outcome_sdY,
    bfile = bfile,
    methods = methods,
    p1 = p1,
    p2 = p2,
    p12 = p12
  )

  # --- Filter exposure to window --------------------------------------------

  chr_val <- as.character(gene_chr) |>
    stringr::str_remove("^chr")
  min_pos <- gene_start - coloc_window
  max_pos <- gene_end + coloc_window

  exposure_filt <- exposure |>
    dplyr::filter(
      as.character(.data$chr.exposure) == chr_val,
      .data$pos.exposure >= min_pos,
      .data$pos.exposure <= max_pos
    )

  if (nrow(exposure_filt) == 0) {
    cli::cli_warn(
      "No exposure SNPs in region chr{chr_val}:{min_pos}-{max_pos}."
    )
    return(new_coloc_result(
      status = "no_snps_in_region",
      status_reason = paste0(
        "No exposure SNPs in region chr",
        chr_val,
        ":",
        min_pos,
        "-",
        max_pos
      ),
      params = params
    ))
  }

  # --- Format outcome to window ---------------------------------------------

  outcome_in_window <- outcome |>
    dplyr::filter(
      as.character(.data$chr) == chr_val,
      .data$pos >= min_pos,
      .data$pos <= max_pos
    )

  if (nrow(outcome_in_window) == 0) {
    cli::cli_warn("No outcome SNPs in region chr{chr_val}:{min_pos}-{max_pos}.")
    return(new_coloc_result(
      status = "no_snps_in_region",
      status_reason = paste0(
        "No outcome SNPs in region chr",
        chr_val,
        ":",
        min_pos,
        "-",
        max_pos
      ),
      params = params
    ))
  }

  outcome_data <- TwoSampleMR::format_data(
    outcome_in_window,
    type = "outcome",
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

  # --- Harmonise ------------------------------------------------------------

  harmonised <- harmonise_and_filter(exposure_filt, outcome_data)

  if (nrow(harmonised) < 3) {
    status <- if (nrow(harmonised) == 0) {
      "no_harmonised_variants"
    } else {
      "too_few_snps"
    }
    reason <- if (nrow(harmonised) == 0) {
      "No variants remaining after harmonisation"
    } else {
      paste0(
        "Only ",
        nrow(harmonised),
        " SNP(s) after harmonisation (need >= 3)"
      )
    }

    cli::cli_warn(reason)
    return(new_coloc_result(
      n_snps = nrow(harmonised),
      harmonised_data = harmonised,
      status = status,
      status_reason = reason,
      params = params
    ))
  }

  # --- Resolve sample sizes -------------------------------------------------

  exp_n <- resolve_sample_size(
    explicit_n = exposure_n,
    data_column = harmonised$samplesize.exposure,
    label = "exposure"
  )
  if (is.null(exp_n)) {
    cli::cli_abort(
      "Exposure sample size could not be determined. Provide {.arg exposure_n}."
    )
  }

  out_n <- resolve_sample_size(
    explicit_n = outcome_n,
    data_column = harmonised$samplesize.outcome,
    label = "outcome"
  )
  if (is.null(out_n)) {
    cli::cli_abort(
      "Outcome sample size could not be determined. Provide {.arg outcome_n}."
    )
  }

  # --- LD matrix ------------------------------------------------------------

  ld_mat <- compute_ld_matrix(
    snps = harmonised$SNP,
    bfile = bfile,
    plink_bin = plink_bin
  )
  aligned <- align_to_ld_matrix(harmonised, ld_mat)
  harmonised <- aligned$data
  ld_mat <- aligned$ld_matrix

  if (nrow(harmonised) < 3) {
    cli::cli_warn(
      "Only {nrow(harmonised)} SNP(s) after LD alignment (need >= 3)."
    )
    return(new_coloc_result(
      n_snps = nrow(harmonised),
      harmonised_data = harmonised,
      status = "too_few_snps",
      status_reason = paste0(
        "Only ",
        nrow(harmonised),
        " SNP(s) after LD alignment (need >= 3)"
      ),
      params = params
    ))
  }

  # --- MAF ------------------------------------------------------------------

  maf <- eaf_to_maf(harmonised$eaf.exposure)

  # --- Build coloc datasets -------------------------------------------------

  dataset_exp <- list(
    beta = harmonised$beta.exposure,
    varbeta = harmonised$se.exposure^2,
    N = exp_n,
    type = exposure_type,
    snp = harmonised$SNP,
    MAF = maf,
    LD = ld_mat
  )
  if (exposure_type == "quant") {
    dataset_exp$sdY <- exposure_sdY
  } else {
    dataset_exp$s <- exposure_s
  }

  dataset_out <- list(
    beta = harmonised$beta.outcome,
    varbeta = harmonised$se.outcome^2,
    N = out_n,
    type = outcome_type,
    snp = harmonised$SNP,
    MAF = maf,
    LD = ld_mat
  )
  if (outcome_type == "quant") {
    dataset_out$sdY <- outcome_sdY
  } else {
    dataset_out$s <- outcome_s
  }

  # --- Run methods ----------------------------------------------------------

  n_snps <- nrow(harmonised)
  coloc_abf_res <- NULL
  coloc_susie_res <- NULL
  coloc_signals_res <- NULL
  coloc_prop_test_res <- NULL
  methods_skipped <- character()

  # ABF
  if ("abf" %in% methods) {
    coloc_abf_res <- tryCatch(
      coloc::coloc.abf(
        dataset1 = dataset_exp,
        dataset2 = dataset_out,
        p1 = p1,
        p2 = p2,
        p12 = p12
      ),
      error = function(e) {
        cli::cli_warn("coloc.abf failed: {conditionMessage(e)}")
        methods_skipped["abf"] <<- conditionMessage(e)
        NULL
      }
    )
  }

  # SuSiE
  susie_exp <- NULL
  susie_out <- NULL
  if ("susie" %in% methods) {
    susie_result <- tryCatch(
      {
        s_exp <- coloc::runsusie(dataset_exp)
        s_out <- coloc::runsusie(dataset_out)
        cs_res <- coloc::coloc.susie(s_exp, s_out)
        list(susie_exp = s_exp, susie_out = s_out, coloc_susie = cs_res)
      },
      error = function(e) {
        cli::cli_warn("SuSiE/coloc.susie failed: {conditionMessage(e)}")
        methods_skipped["susie"] <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(susie_result)) {
      susie_exp <- susie_result$susie_exp
      susie_out <- susie_result$susie_out
      coloc_susie_res <- susie_result$coloc_susie
    }
  }

  # Signals (requires SuSiE)
  if ("signals" %in% methods) {
    if (!"susie" %in% methods) {
      cli::cli_warn(
        "{.val signals} requires {.val susie}; skipping because {.val susie} was not requested."
      )
      methods_skipped["signals"] <- "susie was not requested"
    } else if (is.null(susie_exp) || is.null(susie_out)) {
      cli::cli_warn(
        "{.val signals} requires SuSiE output; skipping because SuSiE failed."
      )
      methods_skipped["signals"] <- "SuSiE failed"
    } else {
      coloc_signals_res <- tryCatch(
        coloc::coloc.signals(
          susie_exp,
          susie_out,
          p1 = p1,
          p2 = p2,
          p12 = p12
        ),
        error = function(e) {
          cli::cli_warn("coloc.signals failed: {conditionMessage(e)}")
          methods_skipped["signals"] <<- conditionMessage(e)
          NULL
        }
      )
    }
  }

  # Prop test (requires signals + colocPropTest)
  if ("prop_test" %in% methods) {
    if (!"signals" %in% methods) {
      cli::cli_warn(
        "{.val prop_test} requires {.val signals}; skipping because {.val signals} was not requested."
      )
      methods_skipped["prop_test"] <- "signals was not requested"
    } else if (is.null(coloc_signals_res)) {
      cli::cli_warn(
        "{.val prop_test} requires coloc.signals output; skipping because signals failed."
      )
      methods_skipped["prop_test"] <- "coloc.signals failed or was skipped"
    } else if (
      !is.data.frame(coloc_signals_res$summary) ||
        nrow(coloc_signals_res$summary) == 0
    ) {
      cli::cli_warn(
        "{.val prop_test} skipped: no signal pairs in coloc.signals output."
      )
      methods_skipped["prop_test"] <- "No signal pairs in coloc.signals output"
    } else if (!requireNamespace("colocPropTest", quietly = TRUE)) {
      cli::cli_warn(
        "{.pkg colocPropTest} package not installed; skipping prop_test."
      )
      methods_skipped["prop_test"] <- "colocPropTest package not installed"
    } else {
      coloc_prop_test_res <- tryCatch(
        {
          prop_fn <- utils::getFromNamespace(
            "coloc.prop.test",
            "colocPropTest"
          )
          prop_fn(coloc_signals_res)
        },
        error = function(e) {
          cli::cli_warn("colocPropTest failed: {conditionMessage(e)}")
          methods_skipped["prop_test"] <<- conditionMessage(e)
          NULL
        }
      )
    }
  }

  # --- Return ---------------------------------------------------------------

  new_coloc_result(
    coloc_abf = coloc_abf_res,
    coloc_susie = coloc_susie_res,
    coloc_signals = coloc_signals_res,
    coloc_prop_test = coloc_prop_test_res,
    n_snps = n_snps,
    harmonised_data = harmonised,
    methods_skipped = methods_skipped,
    params = params
  )
}
