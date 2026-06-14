# -- mrpipeline single-cell eQTL targets pipeline ----------------------------
#
# Runs cis-MR and colocalization for OneK1K and 1M-scBloodNL eQTL data
# against one or more outcome GWAS.
#
# Usage (R console):
#   targets::tar_make()                   # run all
#   targets::tar_make(onek1k_mr_summary)  # run a specific target
#   targets::tar_visnetwork()             # inspect DAG
#   targets::tar_read(onek1k_mr_summary)  # retrieve a result
#
# Before running: edit config.R to set your file paths and parameters.

library(targets)
library(tarchetypes)

source("config.R")

tar_option_set(
  packages = c(
    "mrpipeline",
    "TwoSampleMR",
    "data.table",
    "dplyr",
    "ggplot2"
  ),
  format = "qs" # fast serialisation; install.packages("qs") if needed
)


# =============================================================================
# Helper functions
# =============================================================================

# -- Outcome loading ----------------------------------------------------------

load_outcome_gwas <- function(outcome_cfg) {
  data.table::fread(outcome_cfg$file, nThread = parallel::detectCores()) |>
    as.data.frame()
}


# -- MR batch runner ----------------------------------------------------------
# Loops over all gene x cell_type phenotypes in the exposure data frame and
# calls run_mr() for each. Returns a named list of mr_result objects.

run_mr_batch <- function(exposure, outcome, outcome_id, bfile, cfg) {
  phenotypes <- unique(exposure$exposure)
  if (length(phenotypes) == 0L) {
    return(list())
  }

  results <- lapply(phenotypes, function(pheno) {
    tryCatch(
      mrpipeline::run_mr(
        exposure = exposure,
        exposure_id = pheno,
        outcome = outcome,
        outcome_id = outcome_id,
        bfile = bfile,
        pval_thresh = cfg$mr_pval_thresh,
        rsq_thresh = cfg$mr_rsq_thresh,
        window = cfg$mr_window,
        methods = cfg$mr_methods,
        exclude_regions = cfg$mhc_region
      ),
      error = function(e) {
        cli::cli_warn("MR failed for {pheno}: {conditionMessage(e)}")
        NULL
      }
    )
  })
  names(results) <- phenotypes
  results
}


# -- Summarise a list of mr_result objects ------------------------------------
# Returns a tidy data frame with one row per method per phenotype.

summarise_mr_list <- function(mr_list) {
  rows <- lapply(names(mr_list), function(pheno) {
    res <- mr_list[[pheno]]
    if (is.null(res) || res$status != "success" || nrow(res$results) == 0L) {
      return(data.frame(
        exposure_id = pheno,
        method = NA_character_,
        b = NA_real_,
        se = NA_real_,
        pval = NA_real_,
        n_snps = 0L,
        status = if (is.null(res)) "error" else res$status,
        stringsAsFactors = FALSE
      ))
    }
    r <- res$results
    data.frame(
      exposure_id = pheno,
      method = r$method,
      b = r$b,
      se = r$se,
      pval = r$pval,
      n_snps = nrow(res$instruments),
      status = "success",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# -- Filter MR summary to significant hits for coloc -------------------------
# Returns the subset of phenotypes where at least one method passes the
# coloc_mr_pval threshold.

filter_mr_for_coloc <- function(mr_summary, coloc_mr_pval) {
  mr_summary |>
    dplyr::filter(!is.na(.data$pval) & .data$pval < .env$coloc_mr_pval) |>
    dplyr::pull(.data$exposure_id) |>
    unique()
}


# -- Extract and format OneK1K cis data for a single gene + cell type ---------
# Streams from the 21 GB eqtl_table.tsv.gz via awk to avoid loading the whole
# file. Returns a TwoSampleMR exposure data frame ready for run_coloc().
#
# gene_ct: character in the form "GENE___cell_type" (output of run_mr_batch)

format_onek1k_cis_for_coloc <- function(gene_ct, eqtl_file, onek1k_n) {
  parts <- strsplit(gene_ct, "___")[[1]]
  gene_name <- parts[1]
  cell_type <- parts[2]

  # Column positions in eqtl_table.tsv.gz:
  #   1=CELL_ID  5=GENE  7=CHR  (confirmed from header)
  cmd <- paste0(
    "gunzip -c '",
    eqtl_file,
    "' | ",
    "awk -F'\\t' 'NR==1 || ($5==\"",
    gene_name,
    "\" && $1==\"",
    cell_type,
    "\")'"
  )

  dt <- data.table::fread(
    cmd = cmd,
    nThread = parallel::detectCores()
  ) |>
    as.data.frame()

  if (nrow(dt) == 0L) {
    cli::cli_abort("No cis data found for {gene_name} / {cell_type}.")
  }

  dt <- dt |>
    dplyr::mutate(
      se_derived = abs(.data$SPEARMANS_RHO) /
        stats::qnorm(.data$P_VALUE / 2, lower.tail = FALSE),
      phenotype_label = paste(.data$GENE, cell_type, sep = "___")
    ) |>
    dplyr::filter(is.finite(.data$se_derived) & .data$se_derived > 0)

  n_ct <- onek1k_n[cell_type]
  if (!is.na(n_ct)) {
    dt$samplesize_col <- as.integer(n_ct)
    ss_col <- "samplesize_col"
  } else {
    ss_col <- NULL
  }

  TwoSampleMR::format_data(
    dat = dt,
    type = "exposure",
    phenotype_col = "phenotype_label",
    snp_col = "RSID",
    beta_col = "SPEARMANS_RHO",
    se_col = "se_derived",
    eaf_col = "A2_FREQ_ONEK1K",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    pval_col = "P_VALUE",
    chr_col = "CHR",
    pos_col = "POS",
    samplesize_col = ss_col
  )
}


# -- Run coloc for a single gene x cell_type pair ----------------------------

run_onek1k_coloc <- function(
  gene_ct,
  eqtl_file,
  onek1k_n,
  outcome,
  outcome_id,
  bfile,
  cfg
) {
  coloc_exp <- format_onek1k_cis_for_coloc(gene_ct, eqtl_file, onek1k_n)
  gene_chr <- as.character(unique(coloc_exp$chr.exposure)[1])
  gene_start <- min(coloc_exp$pos.exposure, na.rm = TRUE)
  gene_end <- max(coloc_exp$pos.exposure, na.rm = TRUE)

  # Derive exposure_n from samplesize column if present
  exposure_n <- if ("samplesize.exposure" %in% names(coloc_exp)) {
    as.integer(coloc_exp$samplesize.exposure[1])
  } else {
    NULL
  }

  tryCatch(
    mrpipeline::run_coloc(
      exposure = coloc_exp,
      exposure_id = gene_ct,
      outcome = outcome,
      outcome_id = outcome_id,
      gene_chr = gene_chr,
      gene_start = gene_start,
      gene_end = gene_end,
      coloc_window = cfg$coloc_window,
      exposure_n = exposure_n,
      outcome_type = cfg$outcomes$ra_eur$type,
      outcome_n = cfg$outcomes$ra_eur$n,
      outcome_s = cfg$outcomes$ra_eur$s,
      bfile = bfile,
      methods = cfg$coloc_methods
    ),
    error = function(e) {
      cli::cli_warn("Coloc failed for {gene_ct}: {conditionMessage(e)}")
      NULL
    }
  )
}


# -- Summarise a list of coloc_result objects ---------------------------------

summarise_coloc_list <- function(coloc_list) {
  rows <- lapply(names(coloc_list), function(gene_ct) {
    res <- coloc_list[[gene_ct]]
    if (is.null(res) || res$status != "success") {
      return(data.frame(
        exposure_id = gene_ct,
        n_snps = 0L,
        pp_h4_abf = NA_real_,
        pp4_ratio = NA_real_,
        status = if (is.null(res)) "error" else res$status,
        stringsAsFactors = FALSE
      ))
    }
    pp_h4 <- pp_h3 <- NA_real_
    if (!is.null(res$coloc_abf)) {
      s <- res$coloc_abf$summary
      pp_h4 <- s["PP.H4.abf"]
      pp_h3 <- s["PP.H3.abf"]
    }
    denom <- pp_h3 + pp_h4
    data.frame(
      exposure_id = gene_ct,
      n_snps = res$n_snps,
      pp_h4_abf = as.numeric(pp_h4),
      pp4_ratio = if (!is.na(denom) && denom > 1e-6) {
        as.numeric(pp_h4 / denom)
      } else {
        NA_real_
      },
      status = "success",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# -- Plot helpers -------------------------------------------------------------
# Each saves all plots for a named list of results to a directory and returns
# the vector of file paths (so targets treats it as a file target).

save_mr_plots <- function(mr_list, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saved <- character()
  for (pheno in names(mr_list)) {
    res <- mr_list[[pheno]]
    if (is.null(res) || res$status != "success") {
      next
    }
    safe <- gsub("[^A-Za-z0-9_]", "_", pheno)
    for (type in c("scatter", "forest", "funnel")) {
      p <- tryCatch(plot(res, type = type), error = function(e) NULL)
      if (is.null(p)) {
        next
      }
      # scatter returns a list of ggplots; take the first
      if (is.list(p) && !inherits(p, "ggplot")) {
        p <- p[[1]]
      }
      if (!inherits(p, "ggplot")) {
        next
      }
      path <- file.path(outdir, paste0(safe, "_", type, ".png"))
      ggplot2::ggsave(path, p, width = 8, height = 6)
      saved <- c(saved, path)
    }
  }
  saved
}

save_coloc_plots <- function(coloc_list, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saved <- character()
  for (gene_ct in names(coloc_list)) {
    res <- coloc_list[[gene_ct]]
    if (is.null(res) || res$status != "success") {
      next
    }
    safe <- gsub("[^A-Za-z0-9_]", "_", gene_ct)
    for (type in c("pp_bar", "regional")) {
      p <- tryCatch(plot(res, type = type), error = function(e) NULL)
      if (is.null(p) || !inherits(p, "ggplot")) {
        next
      }
      path <- file.path(outdir, paste0(safe, "_", type, ".png"))
      ggplot2::ggsave(path, p, width = 8, height = 6)
      saved <- c(saved, path)
    }
  }
  saved
}


# =============================================================================
# Pipeline
# =============================================================================

# Static branching values: one row per cell type per dataset
onek1k_values <- tibble::tibble(cell_type = cfg$onek1k_cell_types)
scbloodnl_values <- tibble::tibble(cell_type = cfg$scbloodnl_cell_types)

list(
  # -- Outcome data (loaded once, shared across all analyses) ----------------
  tar_target(
    outcome_ra,
    load_outcome_gwas(cfg$outcomes$ra_eur)
  ),

  # ==========================================================================
  # OneK1K
  # ==========================================================================

  # -- Format MR exposures: one target per cell type, reading esnp_table -----
  # esnp_table is small enough to load fully; CELL_ID column is used to filter.
  tarchetypes::tar_map(
    values = onek1k_values,
    names = "cell_type",

    tar_target(
      onek1k_exposure_mr,
      mrpipeline::format_single_cell_onek1k(
        onek1k_mapping = data.frame(
          cell_type = cell_type,
          path_to_eqtl_file = cfg$onek1k_esnp_file
        ),
        onek1k_cell_type = cell_type,
        sample_size = cfg$onek1k_n[cell_type]
      )
    ),

    # -- Run MR for all genes in this cell type --------------------------------
    tar_target(
      onek1k_mr_ra,
      run_mr_batch(
        exposure = onek1k_exposure_mr,
        outcome = outcome_ra,
        outcome_id = cfg$outcomes$ra_eur$id,
        bfile = cfg$bfile,
        cfg = cfg
      )
    ),

    # -- Tidy summary table for this cell type --------------------------------
    tar_target(
      onek1k_mr_ra_summary,
      summarise_mr_list(onek1k_mr_ra)
    ),

    # -- MR plots (scatter / forest / funnel) for successful results ----------
    tar_target(
      onek1k_mr_ra_plots,
      save_mr_plots(
        onek1k_mr_ra,
        outdir = file.path(cfg$output_dir, "plots", "mr", "onek1k", cell_type)
      ),
      format = "file"
    )
  ),

  # -- Aggregate all OneK1K MR summaries across cell types -------------------
  tar_target(
    onek1k_mr_summary_all,
    dplyr::bind_rows(
      list(
        bin = onek1k_mr_ra_summary_bin,
        bmem = onek1k_mr_ra_summary_bmem,
        cd4et = onek1k_mr_ra_summary_cd4et,
        cd4nc = onek1k_mr_ra_summary_cd4nc,
        cd4sox4 = onek1k_mr_ra_summary_cd4sox4,
        cd8et = onek1k_mr_ra_summary_cd8et,
        cd8nc = onek1k_mr_ra_summary_cd8nc,
        cd8s100b = onek1k_mr_ra_summary_cd8s100b,
        dc = onek1k_mr_ra_summary_dc,
        monoc = onek1k_mr_ra_summary_monoc,
        mononc = onek1k_mr_ra_summary_mononc,
        nk = onek1k_mr_ra_summary_nk,
        nkr = onek1k_mr_ra_summary_nkr,
        plasma = onek1k_mr_ra_summary_plasma
      ),
      .id = "cell_type"
    )
  ),

  # -- Identify gene x cell_type pairs to take to coloc ----------------------
  # (any MR method p < coloc_mr_pval, excluding failed/no-instrument results)
  tar_target(
    onek1k_coloc_candidates,
    filter_mr_for_coloc(onek1k_mr_summary_all, cfg$coloc_mr_pval)
  ),

  # -- Coloc: extract cis data from 21 GB file and run run_coloc() -----------
  # Runs sequentially over candidates; this is the most expensive step.
  tar_target(
    onek1k_coloc_ra_list,
    {
      out <- lapply(onek1k_coloc_candidates, function(gene_ct) {
        run_onek1k_coloc(
          gene_ct = gene_ct,
          eqtl_file = cfg$onek1k_eqtl_file,
          onek1k_n = cfg$onek1k_n,
          outcome = outcome_ra,
          outcome_id = cfg$outcomes$ra_eur$id,
          bfile = cfg$bfile,
          cfg = cfg
        )
      })
      names(out) <- onek1k_coloc_candidates
      out
    }
  ),

  # -- Coloc summary table ---------------------------------------------------
  tar_target(
    onek1k_coloc_summary,
    summarise_coloc_list(onek1k_coloc_ra_list)
  ),

  # -- Coloc plots (pp_bar + regional) for all results ----------------------
  tar_target(
    onek1k_coloc_ra_plots,
    save_coloc_plots(
      onek1k_coloc_ra_list,
      outdir = file.path(cfg$output_dir, "plots", "coloc", "onek1k")
    ),
    format = "file"
  ),

  # ==========================================================================
  # 1M-scBloodNL
  # NOTE: coloc is skipped -- the public FDR-level files do not contain full
  # cis summary statistics required for valid colocalization.
  # ==========================================================================

  tarchetypes::tar_map(
    values = scbloodnl_values,
    names = "cell_type",

    tar_target(
      scbloodnl_exposure_mr,
      mrpipeline::format_sceqtl_1m_scbloodnl(
        file = file.path(
          cfg$scbloodnl_dir,
          cfg$scbloodnl_condition,
          paste0(cell_type, "_expression_eQTLsFDR-ProbeLevel.txt.gz")
        ),
        cell_type = cell_type,
        cis_only = TRUE
      )
    ),

    tar_target(
      scbloodnl_mr_ra,
      run_mr_batch(
        exposure = scbloodnl_exposure_mr,
        outcome = outcome_ra,
        outcome_id = cfg$outcomes$ra_eur$id,
        bfile = cfg$bfile,
        cfg = cfg
      )
    ),

    tar_target(
      scbloodnl_mr_ra_summary,
      summarise_mr_list(scbloodnl_mr_ra)
    ),

    tar_target(
      scbloodnl_mr_ra_plots,
      save_mr_plots(
        scbloodnl_mr_ra,
        outdir = file.path(
          cfg$output_dir,
          "plots",
          "mr",
          "scbloodnl",
          cell_type
        )
      ),
      format = "file"
    )
  ),

  # -- Aggregate all 1M-scBloodNL MR summaries -------------------------------
  tar_target(
    scbloodnl_mr_summary_all,
    dplyr::bind_rows(
      list(
        B = scbloodnl_mr_ra_summary_B,
        CD4T = scbloodnl_mr_ra_summary_CD4T,
        CD8T = scbloodnl_mr_ra_summary_CD8T,
        DC = scbloodnl_mr_ra_summary_DC,
        NK = scbloodnl_mr_ra_summary_NK,
        monocyte = scbloodnl_mr_ra_summary_monocyte,
        megakaryocyte = scbloodnl_mr_ra_summary_megakaryocyte
      ),
      .id = "cell_type"
    )
  ),

  # ==========================================================================
  # Combined results across both datasets
  # ==========================================================================

  tar_target(
    combined_mr_summary,
    dplyr::bind_rows(
      dplyr::mutate(onek1k_mr_summary_all, dataset = "onek1k"),
      dplyr::mutate(scbloodnl_mr_summary_all, dataset = "1m_scbloodnl")
    )
  ),

  # -- Save combined summary to disk -----------------------------------------
  tar_target(
    combined_mr_csv,
    {
      dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(cfg$output_dir, "combined_mr_summary.csv")
      utils::write.csv(combined_mr_summary, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  tar_target(
    onek1k_coloc_csv,
    {
      path <- file.path(cfg$output_dir, "onek1k_coloc_summary.csv")
      utils::write.csv(onek1k_coloc_summary, path, row.names = FALSE)
      path
    },
    format = "file"
  )
)
