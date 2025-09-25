#' Run a Colocalization Pipeline with Optional SuSiE and colocPropTest
#'
#' This function performs colocalization analysis using coloc.abf,
#' optionally integrates SuSiE results for colocalization (using coloc.susie and coloc.signals),
#' and can run colocPropTest on the results of coloc.signals.
#'
#' @param exposure_data A data frame containing pre-formatted exposure GWAS summary statistics.
#'                      Expected to have columns like SNP, beta.exposure, se.exposure,
#'                      eaf.exposure, effect_allele.exposure, other_allele.exposure,
#'                      pval.exposure, chr.exposure, pos.exposure, samplesize.exposure (optional).
#'                      This can be the output of functions like `format_pqtl_decode` or `format_pqtl_ukbppp` or `format_single_cell_onek1k`.
#' @param outcome_gwas A data frame containing the outcome GWAS summary statistics file.
#' @param exposure_name Character string, name for the exposure trait.
#' @param outcome_name Character string, name for the outcome trait.
#' @param exposure_type Character, type of exposure trait ("quant" or "cc").
#' @param outcome_type Character, type of outcome trait ("quant" or "cc").
#' @param exposure_n Integer, sample size for the exposure GWAS. If not provided and
#'                   `exposure_data` contains a `samplesize.exposure` column, that will be used.
#' @param outcome_n Integer, sample size for the outcome GWAS.
#' @param outcome_s Numeric, proportion of cases if outcome_type is "cc". Default NULL.
#' @param exposure_sdY Numeric, standard deviation of the exposure trait if exposure_type is "quant".
#'                     Assumed to be 1 if betas are standardized. Default 1.
#'
#' @param out_snp_col Character, SNP ID column name in outcome data file.
#' @param out_beta_col Character, beta column name in outcome data file.
#' @param out_se_col Character, standard error column name in outcome data file.
#' @param out_eaf_col Character, effect allele frequency column name in outcome data file.
#' @param out_effect_allele_col Character, effect allele column name in outcome data file.
#' @param out_other_allele_col Character, other allele column name in outcome data file.
#' @param out_pval_col Character, p-value column name in outcome data file.
#' @param out_chr_col Character, chromosome column name in outcome data file.
#' @param out_pos_col Character, position column name in outcome data file.
#' @param out_n_col Character, per-SNP sample size column in outcome data file (optional).
#'
#' @param perform_coloc_abf Logical, whether to run `coloc.abf`.
#' @param perform_susie Logical, whether to run SuSiE, `coloc.susie`, and `coloc.signals`.
#' @param perform_coloc_prop_test Logical, whether to run `colocPropTest`.
#'                                 Requires `perform_susie` to be TRUE.
#' @param mhc_remove Logical, whether to remove the MHC region
#' @param gene_chr Character or Integer, chromosome of the gene/region of interest.
#' @param gene_start Integer, start position of the gene/region.
#' @param gene_end Integer, end position of the gene/region.
#' @param window_kb Integer, window in kilobases to extend around the gene/region.
#' @param ld_bfile_path Path to PLINK bfile (prefix, without .bed/.bim/.fam) for LD reference.
#' @param plink_bin_path Path to PLINK executable. Auto-detected if not provided.
#'
#' @param coloc_p1 Numeric, prior probability P(H1).
#' @param coloc_p2 Numeric, prior probability P(H2).
#' @param coloc_p12 Numeric, prior probability P(H4).
#'
#' @param output_dir Path to a directory where results and logs might be saved.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A list containing results from the requested analyses.
#'         Elements may include: `harmonized_data_aligned`, `ld_matrix`,
#'         `coloc_abf_results`, `susie_exposure`, `susie_outcome`,
#'         `coloc_susie_results`, `coloc_signals_results`, `coloc_prop_test_results`.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create dummy exposure data (already formatted)
#' dummy_exposure_data <- data.frame(
#'   SNP = paste0("rs", 1:100),
#'   chr.exposure = 1,
#'   pos.exposure = (1:100) * 1000 + 10000, # Ensure positions are within region
#'   effect_allele.exposure = "A",
#'   other_allele.exposure = "G",
#'   beta.exposure = rnorm(100, 0, 0.1),
#'   se.exposure = runif(100, 0.01, 0.05),
#'   pval.exposure = runif(100),
#'   eaf.exposure = runif(100, 0.01, 0.99),
#'   samplesize.exposure = 10000,
#'   phenotype.exposure = "MyExposure"
#' )
#'
#' # Create dummy outcome GWAS file
#' temp_dir <- tempdir()
#' outcome_file <- file.path(temp_dir, "outcome_gwas.txt")
#' ld_bfile <- file.path(temp_dir, "ld_ref") # Needs actual .bed, .bim, .fam
#'
#' write.table(data.frame(
#'   SNPID = paste0("rs", 1:150), CHR = 1, POS = (1:150) * 1000 + 10000,
#'   EA = "A", NEA = "G", BETA_outcome = rnorm(150, 0, 0.1),
#'   SE_outcome = runif(150, 0.01, 0.05), P_outcome = runif(150),
#'   EAF_outcome = runif(150, 0.01, 0.99), N_outcome = 20000
#' ), outcome_file, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # Create dummy LD reference files (these won't work for actual LD calculation)
#' # file.create(paste0(ld_bfile, c(".bed", ".bim", ".fam")))
#'
#' results <- run_coloc_pipeline(
#'   exposure_data = dummy_exposure_data,
#'   outcome_gwas_path = outcome_file,
#'   exposure_name = "MyExposure",
#'   outcome_name = "MyOutcome",
#'   exposure_type = "quant",
#'   outcome_type = "cc",
#'   exposure_n = 10000, # Can be omitted if samplesize.exposure in exposure_data
#'   outcome_n = 20000,
#'   outcome_s = 0.1,
#'   out_snp_col = "SNPID", out_beta_col = "BETA_outcome", out_se_col = "SE_outcome",
#'   out_eaf_col = "EAF_outcome", out_effect_allele_col = "EA", out_other_allele_col = "NEA",
#'   out_pval_col = "P_outcome", out_chr_col = "CHR", out_pos_col = "POS", out_n_col = "N_outcome",
#'   perform_coloc_abf = TRUE,
#'   perform_susie = FALSE, # Set to TRUE if you have a working LD ref and susieR
#'   perform_coloc_prop_test = FALSE, # Set to TRUE if perform_susie is TRUE
#'   gene_chr = 1,
#'   gene_start = 20000,
#'   gene_end = 80000,
#'   window_kb = 50,
#'   ld_bfile_path = ld_bfile,
#'   output_dir = file.path(temp_dir, "coloc_results")
#' )
#' print(results$coloc_abf_results$summary)
#' }
run_coloc_pipeline <- function(
  exposure_data, outcome_gwas_path,
  exposure_name, outcome_name,
  exposure_type = "quant", outcome_type = "cc",
  exposure_n = NULL, outcome_n, outcome_s = NULL, exposure_sdY = 1,

  out_snp_col, out_beta_col, out_se_col,
  out_eaf_col, out_effect_allele_col, out_other_allele_col,
  out_pval_col, out_chr_col, out_pos_col, out_n_col = NULL,

  perform_coloc_abf = TRUE, perform_susie = TRUE, perform_coloc_prop_test = TRUE,
  mhc_remove = FALSE,
  gene_chr, gene_start, gene_end, window_kb = 200,
  ld_bfile_path,
  plink_bin_path = NULL,

  coloc_p1, coloc_p2, coloc_p12,

  output_dir,
  verbose = TRUE
) {

  # --- 0. Argument Checks and Setup ---
  if (!is.data.frame(exposure_data)) {
    stop("exposure_data must be a pre-formatted data frame.")
  }
  if (perform_coloc_prop_test && !perform_susie) {
    warning("perform_coloc_prop_test requires perform_susie to be TRUE. Setting perform_susie to TRUE.")
    perform_susie <- TRUE
  }
  if (outcome_type == "cc" && is.null(outcome_s)) {
    stop("If outcome_type is 'cc', outcome_s (proportion of cases) must be provided.")
  }
  if (is.null(plink_bin_path)) {
    plink_bin_path <- tryCatch(genetics.binaRies::get_plink_binary(), error = function(e) NULL)
    if (is.null(plink_bin_path) && (perform_susie || file.exists(paste0(ld_bfile_path, ".bed")))) {
        warning("PLINK binary not found or specified. LD matrix calculation might fail.")
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  results_list <- list()

  # --- 1. Prepare Exposure Data (already formatted) ---
  if (verbose) message("Using pre-formatted exposure data...")
  exposure_dat <- exposure_data
  # Ensure phenotype column exists or create it
  if (!"phenotype.exposure" %in% colnames(exposure_dat)) {
      exposure_dat$phenotype.exposure <- exposure_name
  }
  # Determine exposure_n if not provided but samplesize.exposure column exists
  if (is.null(exposure_n) && "samplesize.exposure" %in% colnames(exposure_dat)) {
    exposure_n <- as.integer(stats::median(exposure_dat$samplesize.exposure, na.rm = TRUE))
    if (verbose) message(paste("Using median sample size from exposure_data for exposure_n:", exposure_n))
  } else if (is.null(exposure_n)) {
    stop("exposure_n must be provided if 'samplesize.exposure' column is not in exposure_data.")
  }


  # --- 2. Load and Format Outcome Data ---
  if (verbose) message("Loading and formatting outcome data...")
  outcome_raw <- data.table::fread(outcome_gwas_path)
  outcome_dat_formatted <- TwoSampleMR::format_data(
    outcome_raw, type = "outcome",
    snp_col = out_snp_col, beta_col = out_beta_col, se_col = out_se_col,
    eaf_col = out_eaf_col, effect_allele_col = out_effect_allele_col,
    other_allele_col = out_other_allele_col, pval_col = out_pval_col,
    chr_col = out_chr_col, pos_col = out_pos_col,
    samplesize_col = out_n_col,
    phenotype_col = outcome_name # Creates outcome_name column
  )
  # format_data adds "outcome." prefix, remove it for easier use
  colnames(outcome_dat_formatted) <- sub("^outcome\\.", "", colnames(outcome_dat_formatted))
  outcome_dat_formatted$phenotype <- outcome_name # Ensure phenotype column is correctly named


  # --- 3. Optionally remove MHC, filter by region and Harmonize ---
  # Optional: Filter out MHC region from the already region-filtered data
  if (isTRUE(mhc_remove)) {
    if ("chr.exposure" %in% colnames(exposure_dat)) {
      if (verbose) message("Assessing whether MHC region is present and should be removed...")

      mhc_chr <- "6"
      mhc_start <- 26000000
      mhc_end <- 34000000

      # Check if all SNPs in the region are within MHC
      if (all(exposure_dat$chr.exposure == mhc_chr &
              exposure_dat$pos.exposure >= mhc_start &
              exposure_dat$pos.exposure <= mhc_end)) {
        if (verbose) message("All SNPs for this protein in the specified region are within MHC. Skipping analysis.")
        return(NULL)
      }

      # Filter out any SNPs in MHC region
      initial_snp_count <- nrow(exposure_dat)
      exposure_dat <- exposure_dat |>
        dplyr::filter(!(as.character(chr.exposure) == mhc_chr & pos.exposure >= mhc_start & pos.exposure <= mhc_end))
      if (verbose && nrow(exposure_dat) < initial_snp_count) message("Filtered out SNPs in MHC region (chr6:26-34Mb).")
    } else {
      warning("`mhc_remove` is TRUE, but 'chr.exposure' column not found in exposure data. MHC filter skipped.")
    }
  }

  if (verbose) message("Filtering data by region and harmonizing...")
  current_chr_val <- gsub("chr", "", as.character(gene_chr))
  min_pos <- gene_start - (window_kb * 1000)
  max_pos <- gene_end + (window_kb * 1000)

  # Ensure column names for filtering are correct for pre-formatted exposure_dat
  exposure_filt <- exposure_dat |>
    dplyr::filter(as.character(chr.exposure) == current_chr_val & pos.exposure >= min_pos & pos.exposure <= max_pos)

  outcome_filt <- outcome_dat_formatted |>
    dplyr::filter(chr == current_chr_val & pos >= min_pos & pos <= max_pos)

  if (nrow(exposure_filt) == 0 || nrow(outcome_filt) == 0) {
    warning("No SNPs remaining in one or both datasets after filtering for the region.")
    return(list(error = "No SNPs in region", harmonized_data = NULL))
  }

  # Rename exposure columns to match TwoSampleMR expectations before harmonization
  # TwoSampleMR::harmonise_data expects columns like SNP, beta.exposure, se.exposure etc.
  # The pre-formatted exposure_data should already have these.

  harmonized_data <- TwoSampleMR::harmonise_data(exposure_filt, outcome_filt)
  harmonized_data <- harmonized_data |> dplyr::filter(mr_keep == TRUE)

  if (nrow(harmonized_data) == 0) {
    warning("No SNPs remaining after initial harmonization.")
    return(list(error = "No SNPs post-initial-harmonization", harmonized_data = NULL))
  }

  # --- 4. Prepare LD Matrix & Align with Harmonized Data ---
  ld_matrix <- NULL
  harmonized_data_aligned <- harmonized_data # Keep original harmonized for coloc.abf if LD fails

  if (perform_susie || perform_coloc_abf) {
    if (verbose) message("Preparing LD matrix and aligning data...")
    common_snps_for_ld <- harmonized_data$SNP

    if (length(common_snps_for_ld) > 1 && !is.null(ld_bfile_path) && !is.null(plink_bin_path) && file.exists(paste0(ld_bfile_path, ".bed"))) {
      tryCatch({
        ld_matrix_raw <- ieugwasr::ld_matrix(
          variants = common_snps_for_ld,
          bfile = ld_bfile_path,
          plink_bin = plink_bin_path
        )
        # Clean LD matrix SNP names (rsID_A_G -> rsID)
        # This is important if your harmonized_data$SNP is just rsID
        cleaned_ld_rownames <- gsub("_.*", "", rownames(ld_matrix_raw))
        cleaned_ld_colnames <- gsub("_.*", "", colnames(ld_matrix_raw))
        rownames(ld_matrix_raw) <- cleaned_ld_rownames
        colnames(ld_matrix_raw) <- cleaned_ld_colnames

        # Ensure harmonized_data matches ld_matrix (allele flipping and reordering)
        # This complex block aligns alleles between summary stats and LD matrix
        temp_harmonized_data <- harmonized_data

        # Identify SNPs in harmonized_data that need allele flipping to match LD matrix conventions
        # This assumes LD matrix from ieugwasr might have specific allele orders in its _A_G suffix (though we removed it)
        # The core idea is to match the SNP IDs first.
        snps_in_ld <- intersect(temp_harmonized_data$SNP, rownames(ld_matrix_raw))
        if(length(snps_in_ld) == 0) stop("No SNPs in common between harmonized data and LD matrix after cleaning LD SNP names.")

        temp_harmonized_data <- temp_harmonized_data[temp_harmonized_data$SNP %in% snps_in_ld, ]
        ld_matrix <- ld_matrix_raw[snps_in_ld, snps_in_ld]

        # The following block is for detailed allele alignment if rsID_A_G format was strictly used by LD matrix
        # For now, we assume simple rsID matching after gsub is sufficient for ieugwasr::ld_matrix
        # If more complex allele matching is needed (e.g. if ld_matrix was from plink --r square with --make-just-sq),
        # the provided snippet would be adapted here.
        # The key is that `harmonized_data_aligned` and `ld_matrix` must have SNPs in the same order
        # and alleles oriented consistently.

        # Simplified alignment: ensure order and subset
        harmonized_data_aligned <- temp_harmonized_data[match(rownames(ld_matrix), temp_harmonized_data$SNP), ]
        # Remove any NA rows that might result if a SNP in ld_matrix rownames isn't in harmonized_data_aligned$SNP (should not happen with intersect)
        harmonized_data_aligned <- harmonized_data_aligned[!is.na(harmonized_data_aligned$SNP), ]
        # Re-filter ld_matrix to be absolutely sure
        ld_matrix <- ld_matrix[harmonized_data_aligned$SNP, harmonized_data_aligned$SNP]

        if(nrow(harmonized_data_aligned) == 0 || nrow(ld_matrix) == 0) stop("Data alignment with LD matrix resulted in zero SNPs.")

        results_list$ld_matrix <- ld_matrix
        results_list$harmonized_data_aligned <- harmonized_data_aligned

      }, error = function(e) {
        warning(paste("LD matrix calculation or data alignment failed:", e$message))
        ld_matrix <<- NULL
        harmonized_data_aligned <<- harmonized_data # Fallback to non-LD-aligned data
      })
    } else if (length(common_snps_for_ld) <= 1) {
        warning("Too few SNPs to calculate LD matrix.")
        harmonized_data_aligned <- harmonized_data
    } else {
        warning("LD bfile path/PLINK binary not provided/found or LD file missing. Cannot calculate LD matrix.")
        harmonized_data_aligned <- harmonized_data
    }
  } else {
      harmonized_data_aligned <- harmonized_data # If not using LD
  }

  # Use harmonized_data_aligned for subsequent steps
  current_harmonized_data <- harmonized_data_aligned
  if (nrow(current_harmonized_data) == 0) {
    warning("No SNPs remaining after attempting LD alignment.")
    return(list(error = "No SNPs post-LD-alignment", harmonized_data_aligned = NULL))
  }

  # --- 5. SuSiE Analysis, coloc.signals, and colocPropTest ---
  if (perform_susie) {
    if (verbose) message("Performing SuSiE analysis...")
    if (!is.null(ld_matrix) && nrow(ld_matrix) > 0 && nrow(current_harmonized_data) > 0) {
      if (!"eaf.exposure" %in% colnames(current_harmonized_data)) {
          warning("eaf.exposure not found for SuSiE. SuSiE might be less accurate or fail.")
          current_harmonized_data$eaf.exposure <- 0.1 # Fallback
      }
      current_harmonized_data$MAF <- ifelse(current_harmonized_data$eaf.exposure < 0.5, current_harmonized_data$eaf.exposure, 1 - current_harmonized_data$eaf.exposure)

      dataset_susie_exp <- list(
        beta = current_harmonized_data$beta.exposure,
        varbeta = current_harmonized_data$se.exposure^2,
        N = exposure_n,
        type = exposure_type,
        MAF = current_harmonized_data$MAF,
        LD = ld_matrix,
        snp = current_harmonized_data$SNP,
        sdY = if (exposure_type == "quant") exposure_sdY else NULL
      )
      dataset_susie_out <- list(
        beta = current_harmonized_data$beta.outcome,
        varbeta = current_harmonized_data$se.outcome^2,
        N = outcome_n,
        type = outcome_type,
        s = if (outcome_type == "cc") outcome_s else NULL,
        MAF = current_harmonized_data$MAF, # Assuming same MAF after harmonization
        LD = ld_matrix,
        snp = current_harmonized_data$SNP,
        sdY = if (outcome_type == "quant") 1 else NULL
      )

      tryCatch({
        susie_exp_results <- coloc::runsusie(dataset_susie_exp, verbose = verbose)
        susie_out_results <- coloc::runsusie(dataset_susie_out, verbose = verbose)
        results_list$susie_exposure <- susie_exp_results
        results_list$susie_outcome <- susie_out_results

        if (verbose) message("Performing coloc with SuSiE results (coloc.susie)...")
        coloc_susie_res <- coloc::coloc.susie(susie_exp_results, susie_out_results)
        results_list$coloc_susie_results <- coloc_susie_res

        if (verbose) message("Performing coloc.signals with SuSiE results...")
        # coloc.signals requires p1, p2, p12 similar to coloc.abf
        coloc_signals_res <- coloc::coloc.signals(
            susie_exp_results, susie_out_results,
            p1 = coloc_p1, p2 = coloc_p2, p12 = coloc_p12
        )
        results_list$coloc_signals_results <- coloc_signals_res

        if (perform_coloc_prop_test) {
          if (verbose) message("Performing colocPropTest analysis...")
          if (!is.null(coloc_signals_res) && nrow(coloc_signals_res$summary) > 0) {
            if(requireNamespace("colocPropTest", quietly = TRUE)) {
                # coloc.prop.test expects the output of coloc.signals
                coloc_prop_res <- colocPropTest::coloc.prop.test(coloc_signals_res)
                results_list$coloc_prop_test_results <- coloc_prop_res
            } else {
                warning("colocPropTest package not installed. Skipping this step.")
            }
          } else {
            warning("Skipping colocPropTest: coloc.signals results not available or empty.")
          }
        }

      }, error = function(e) {
        warning(paste("SuSiE, coloc.susie, or coloc.signals analysis failed:", e$message))
      })
    } else {
      warning("Skipping SuSiE-based analyses: LD matrix not available or no overlapping SNPs after alignment.")
    }
  }

  # --- 6. coloc.abf Analysis ---
  if (perform_coloc_abf) {
    if (verbose) message("Performing coloc.abf analysis...")
     if (!"eaf.exposure" %in% colnames(current_harmonized_data) || !"eaf.outcome" %in% colnames(current_harmonized_data)) {
        warning("EAF not found in harmonized data for coloc.abf. Results may be less reliable.")
     }

    dataset1_coloc <- list(
      beta = current_harmonized_data$beta.exposure,
      varbeta = current_harmonized_data$se.exposure^2,
      N = exposure_n,
      type = exposure_type,
      snp = current_harmonized_data$SNP,
      MAF = if ("eaf.exposure" %in% colnames(current_harmonized_data)) current_harmonized_data$eaf.exposure else rep(0.1, nrow(current_harmonized_data)),
      sdY = if (exposure_type == "quant") exposure_sdY else NULL
    )
    dataset2_coloc <- list(
      beta = current_harmonized_data$beta.outcome,
      varbeta = current_harmonized_data$se.outcome^2,
      N = outcome_n,
      type = outcome_type,
      s = if (outcome_type == "cc") outcome_s else NULL,
      snp = current_harmonized_data$SNP,
      MAF = if ("eaf.outcome" %in% colnames(current_harmonized_data)) current_harmonized_data$eaf.outcome else rep(0.1, nrow(current_harmonized_data)),
      sdY = if (outcome_type == "quant") 1 else NULL
    )

    # Add LD to coloc.abf if available and desired (can improve accuracy)
    # Use the aligned LD matrix
    if(!is.null(ld_matrix) && nrow(ld_matrix) > 0 && all(current_harmonized_data$SNP %in% rownames(ld_matrix))) {
        # Ensure dataset SNPs match LD matrix order if LD is used
        ordered_snps <- intersect(rownames(ld_matrix), current_harmonized_data$SNP)
        if (length(ordered_snps) == nrow(current_harmonized_data) && length(ordered_snps) == nrow(ld_matrix)) {
            # This assumes current_harmonized_data is already ordered by LD matrix if ld_matrix is not NULL
            dataset1_coloc$LD <- ld_matrix[ordered_snps, ordered_snps]
            dataset2_coloc$LD <- ld_matrix[ordered_snps, ordered_snps]
        } else {
            warning("SNP mismatch or order issue when adding LD to coloc.abf datasets. Proceeding without LD for coloc.abf.")
        }
    }

    tryCatch({
      coloc_abf_results_run <- coloc::coloc.abf(
        dataset1 = dataset1_coloc,
        dataset2 = dataset2_coloc,
        p1 = coloc_p1, p2 = coloc_p2, p12 = coloc_p12
      )
      results_list$coloc_abf_results <- coloc_abf_results_run
    }, error = function(e) {
      warning(paste("coloc.abf analysis failed:", e$message))
    })
  }

  if (verbose) message("Pipeline finished.")
  return(results_list)
}
