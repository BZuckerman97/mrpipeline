#' Performs MR
#'
#' @importFrom utils timestamp
#'
#' @param exposure_id character. e.g. "syn52361761" -> we should change this to the protein/gene name
#' @param exposure data frame. Summary statistics for exposure, must include exposure_id column ie gene name which we use as above?
#' @param outcome data frame. Summary statistics for outcome, must be formatted as above
#' @param outcome_id character. Name for outcome e.g. 'SjD'
#' @param instrument_region list of the chromosome position, gene start and gene end for each gene of interest this needs to link with the original mapping file
#' @param window integer. Set to cis region if required
#' @param pval_thresh number. 5e-6 by default
#' @param rsq_thresh R square clumping threshold
#' @param bfile path to the LD_folder containing 1000 genome files (g1000_eur)
#'
#' @return List with 2 data frames - MR results and instruments
#' @export
#' @examples

run_mr <- function(exposure,
                   exposure_id,
                   outcome, #' this needs to be the formatted summary statistics for the outcome GWAS
                   outcome_id,
                   # sumstats_info,
                   # downloadLocation,
                   # ref_rsid,
                   instrument_region = list(chromosome = 1L,
                                            start = 1L,
                                            end = 200L),
                   window = 1L,
                   pval_thresh = 5e-6,
                   rsq_thresh = 0.1,
                   bfile) {

# Validate arguments -----------------------------------------------------
# Note: The validation function is currently commented out.
# It should be implemented to check the structure and types of instrument_region.
validate_instrument_region_arg(instrument_region)

  result <- NULL

  print(exposure_id)
  timestamp()
  print("******")

# Filter exposure for IVs ----------------------------------------------------------

  #' Taken out sumstats_info here
  #' To replace with tidyverse language
  #' To move out of the run_mr function and have a specify window
  exposure <- exposure |>
    dplyr::filter(pos.exposure > (instrument_region$start - window) & pos.exposure < (instrument_region$end + window)) |> # Selecting the cis region only (here defined as 200kb before or after the protein-encoding region), uses build 37 positions
    dplyr::filter(pval.exposure < pval_thresh)                                                                                   # Selecting "region-wide" significant cis-pQTLs (here defined as P<5e-6)

  #' This needs to be done before the function starts
  #exposure$CHROM <-
    #  ifelse(exposure$CHROM == 23, "X", exposure$CHROM)                                                                     # Renaming 23rd chromosome "X" for consistency between sum stats


# Get IVs from outcome ----------------------------------------------------

    if (is.null(exposure) || nrow(exposure) == 0) {
      warning(paste0("Skipping ", exposure_id))
      warning("No significant cis pQTLs")
    } else {
      # Only selecting the chromosome of interest to speed up stuff downstream from here
      #' Will need to make sure this is on the format_data() outcome data ie chr.outcome and chr.exposure
      outcome_overlap <- outcome |>
        dplyr::filter(chr %in% exposure$chr.exposure) # |> #' CHROM needs to change to chr.exposure
      #  dplyr::filter(pos %in% exposure$pos.exposure) - TO HIDE THIS UNTIL I CAN TOGGLE BETWEEN BUILD 37 OR BUILD 38 IN THE FORMAT PQTL DATA
      # Only selecting the variants that are overlapping between exposure and outcome sum stats

      if (is.null(outcome_overlap) ||
          nrow(outcome_overlap) == 0) {
        warning(paste0("Skipping ", exposure_id))
        warning("No overlap between outcome and protein exposure")
      } else {


# Reformat outcome df -----------------------------------------------------
        # This step is unnecessary
        # outcome_rsid <- outcome_overlap |>
        #  dplyr::select(chrom, pos, rsids)

        outcome_overlap <- outcome_overlap |>
          dplyr::mutate(phenotype = paste(outcome_id)) #|>
          #dplyr::mutate(id = paste(chrom, pos, alt, ref, sep = ":")) no longer need this as using rsIDs/SNPs not ID

        outcome_overlap <- as.data.frame(outcome_overlap)

        # These column names need to be completed prior to the function
        outcome_data <- TwoSampleMR::format_data(
          outcome_overlap,
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
          log_pval = FALSE
        )

# Harmonise ---------------------------------------------------------------

        harmonised_data_frame <-
          TwoSampleMR::harmonise_data(exposure_dat = exposure, outcome_dat = outcome_data)                                             # This is where the matching happens

        harmonised_data_frame <- harmonised_data_frame |>
          dplyr::arrange(pval.exposure)
# We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
        harmonised_data_frame <- harmonised_data_frame |>
          dplyr::filter(!duplicated(SNP))

        if (nrow(harmonised_data_frame) == 0) {
          warning(paste0("Skipping ", exposure_id))
          warning("No variants remaining after harmonising")
          return(NULL)
        }

# The code below calculates the absolute p-values and generates pseudo p-values to ensure that each SNP is properly ranked by PLINK.
harmonised_data_frame <- harmonised_data_frame %>%
  dplyr::mutate(
    z = beta.exposure / se.exposure,
    log10p = sapply(
      -log10(2 * Rmpfr::pnorm(
        Rmpfr::mpfr(abs(z), precBits = 100),
        lower.tail = FALSE
      )),
      as.numeric
    )
  )

harmonised_data_frame <- harmonised_data_frame %>%
  dplyr::arrange(desc(log10p)) %>%
  dplyr::mutate(pseudo_p = seq(from = 1e-100, to = 0.9, length.out = n()))

# Clump -------------------------------------------------------------------
        print("Clumping")
          clump <-
            ieugwasr::ld_clump(
              dplyr::tibble(
                rsid = harmonised_data_frame$SNP,
                pval = harmonised_data_frame$pseudo_p,
                id = harmonised_data_frame$id.exposure
              ),
              # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
              plink_bin = genetics.binaRies::get_plink_binary(),
              clump_kb = 10000,
              clump_r2 = rsq_thresh,
              bfile = bfile
            )
          harmonised_clumped_final_data_frame <- harmonised_data_frame |>
            dplyr::filter(harmonised_data_frame$SNP %in% clump$rsid)

          # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.

# Perform MR --------------------------------------------------------------

          if (nrow(harmonised_clumped_final_data_frame[harmonised_clumped_final_data_frame$mr_keep,]) == 0) {
            warning(paste0("Skipping ", exposure_id))
            warning("No variants remaining after clumping")
            results <- NULL
          } else {

## Wald ratio --------------------------------------------------------------

            if (nrow(harmonised_clumped_final_data_frame) == 1) {
              # If the genetic instrument includes 1 variant, you use the Wald ratio as your method
              results_mr <-
                TwoSampleMR::mr(harmonised_clumped_final_data_frame, method_list = c("mr_wald_ratio"))
              results <-
                data.frame(
                  exp = exposure_id,
                  outc = paste(outcome_id),
                  pvalthreshold = pval_thresh,
                  rsqthreshold = rsq_thresh,
                  nsnp = results_mr$nsnp,
                  method = results_mr$method,
                  b = results_mr$b,
                  se = results_mr$se,
                  pval = results_mr$pval
                )
            } else if (nrow(harmonised_clumped_final_data_frame) == 2) {


## 2 IVs - IVW not Egger ---------------------------------------------------

              # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
              ld_correlation_matrix <-
                ieugwasr::ld_matrix(
                  harmonised_clumped_final_data_frame$SNP,
                  bfile = bfile,
                  plink_bin = genetics.binaRies::get_plink_binary()
                )

              # To make sure the LD matrix columns match the SNP order precisely
              # Art's code to solve the order issues of SNPid in LD matix
              rownames(ld_correlation_matrix) <- gsub("\\_.*", "", rownames(ld_correlation_matrix))
              harmonised_clumped_final_data_frame_duplicate <- harmonised_clumped_final_data_frame[!(paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_") %in% colnames(ld_correlation_matrix)),]
              colnames(harmonised_clumped_final_data_frame_duplicate)[which(colnames(harmonised_clumped_final_data_frame_duplicate) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
              harmonised_clumped_final_data_frame_duplicate$beta.exposure <- harmonised_clumped_final_data_frame_duplicate$beta.exposure*-1
              harmonised_clumped_final_data_frame_duplicate$beta.outcome <- harmonised_clumped_final_data_frame_duplicate$beta.outcome*-1
              harmonised_clumped_final_data_frame_duplicate$eaf.exposure <- 1-harmonised_clumped_final_data_frame_duplicate$eaf.exposure
              harmonised_clumped_final_data_frame_duplicate$eaf.outcome <- 1-harmonised_clumped_final_data_frame_duplicate$eaf.outcome
              harmonised_clumped_final_data_frame <- harmonised_clumped_final_data_frame[(paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_") %in% colnames(ld_correlation_matrix)),]
              harmonised_clumped_final_data_frame <- rbind(harmonised_clumped_final_data_frame, harmonised_clumped_final_data_frame_duplicate)
              harmonised_clumped_final_data_frame <- harmonised_clumped_final_data_frame[ order(match(harmonised_clumped_final_data_frame$SNP, rownames(ld_correlation_matrix))), ]
              harmonised_clumped_final_data_frame$marker_ld <- paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_")

              harmonised_clumped_correlated_final_data_frame <-
                MendelianRandomization::mr_input( # This variable name is inconsistent, consider renaming
                  bx = harmonised_clumped_final_data_frame$beta.exposure, # Use harmonised_clumped_final_data_frame
                  bxse = harmonised_clumped_final_data_frame$se.exposure,
                  by = harmonised_clumped_final_data_frame$beta.outcome,
                  byse = harmonised_clumped_final_data_frame$se.outcome,
                  correlation = ld_correlation_matrix # Corrected variable name
                )
              output_mr_ivw_corr <-
                MendelianRandomization::mr_ivw(harmonised_clumped_correlated_final_data_frame, correl = TRUE)
              results <-
                data.frame(
                  exp = exposure_id,
                  outc = paste(outcome_id),
                  pvalthreshold = pval_thresh,
                  rsqthreshold = rsq_thresh,
                  nsnp = output_mr_ivw_corr@SNPs,
                  method = "Inverse variance weighted (correlation inc)",
                  b = output_mr_ivw_corr@Estimate,
                  se = output_mr_ivw_corr@StdError,
                  pval = output_mr_ivw_corr@Pvalue
                )
            }
            else { # This block is for nrow > 2
              # Ensure LD matrix is calculated for >2 IVs if not already done
              # (The existing code calculates it in the nrow == 2 block,
              #  it might be better to calculate it once if nrow >= 2)
              # For this fix, we assume ld_correlation_matrix is available or recalculated here if needed.
              # If ld_correlation_matrix from the nrow==2 block is intended to be used,
              # ensure its scope or recalculate. For simplicity, let's assume it needs to be
              # available or recalculated if the logic implies it's different for >2 IVs.
              # The complex SNP/allele ordering logic is also present in the nrow==2 block.
              # This logic should also be applied here if necessary before mr_input.

              # Re-applying the SNP ordering and allele flipping logic for > 2 IVs
              # Calculate ld_correlation_matrix for the >2 IVs case
              ld_correlation_matrix <-
                ieugwasr::ld_matrix(
                  harmonised_clumped_final_data_frame$SNP,
                  bfile = bfile,
                  plink_bin = genetics.binaRies::get_plink_binary()
                )

              # This is duplicated from the nrow == 2 block and should ideally be refactored.
              rownames(ld_correlation_matrix) <- gsub("\\_.*", "", rownames(ld_correlation_matrix))
              harmonised_clumped_final_data_frame_duplicate <- harmonised_clumped_final_data_frame[!(paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_") %in% colnames(ld_correlation_matrix)),]
              colnames(harmonised_clumped_final_data_frame_duplicate)[which(colnames(harmonised_clumped_final_data_frame_duplicate) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
              harmonised_clumped_final_data_frame_duplicate$beta.exposure <- harmonised_clumped_final_data_frame_duplicate$beta.exposure*-1
              harmonised_clumped_final_data_frame_duplicate$beta.outcome <- harmonised_clumped_final_data_frame_duplicate$beta.outcome*-1
              harmonised_clumped_final_data_frame_duplicate$eaf.exposure <- 1-harmonised_clumped_final_data_frame_duplicate$eaf.exposure
              harmonised_clumped_final_data_frame_duplicate$eaf.outcome <- 1-harmonised_clumped_final_data_frame_duplicate$eaf.outcome
              harmonised_clumped_final_data_frame <- harmonised_clumped_final_data_frame[(paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_") %in% colnames(ld_correlation_matrix)),]
              harmonised_clumped_final_data_frame <- rbind(harmonised_clumped_final_data_frame, harmonised_clumped_final_data_frame_duplicate)
              harmonised_clumped_final_data_frame <- harmonised_clumped_final_data_frame[ order(match(harmonised_clumped_final_data_frame$SNP, rownames(ld_correlation_matrix))), ]
              harmonised_clumped_final_data_frame$marker_ld <- paste(harmonised_clumped_final_data_frame$SNP, harmonised_clumped_final_data_frame$effect_allele.exposure, harmonised_clumped_final_data_frame$other_allele.exposure, sep="_")

              harmonised_clumped_correlated_final_data_frame <-
                MendelianRandomization::mr_input(
                  bx = harmonised_clumped_final_data_frame$beta.exposure,
                  bxse = harmonised_clumped_final_data_frame$se.exposure,
                  by = harmonised_clumped_final_data_frame$beta.outcome,
                  byse = harmonised_clumped_final_data_frame$se.outcome,
                  correlation = ld_correlation_matrix
                )
              output_mr_ivw_corr <- # Define for >2 IVs
                MendelianRandomization::mr_ivw(harmonised_clumped_correlated_final_data_frame, correl = TRUE)

              # Format IVW results
              results <- data.frame(
                exp = exposure_id, outc = paste(outcome_id), pvalthreshold = pval_thresh, rsqthreshold = rsq_thresh,
                nsnp = output_mr_ivw_corr@SNPs, method = "Inverse variance weighted (correlation inc)",
                b = output_mr_ivw_corr@Estimate, se = output_mr_ivw_corr@StdError, pval = output_mr_ivw_corr@Pvalue
              )
            }
          }

          if (is.null(results) || nrow(results) == 0) {
            warning(paste0("Skipping ", exposure_id))
            warning("No results returned from MR analysis")
          } else {
            df_sum <-
              data.frame(
                exp = as.character(NA),
                outc = as.character(NA),
                nsnp = NA,
                method = NA,
                b = NA,
                se = NA,
                pval = NA
              )[-1,]
            df_instr <-
              data.frame(
                pos.exposure = NA,
                pos_id = NA,
                effect_allele.exposure = NA,
                other_allele.exposure = NA,
                effect_allele.outcome = NA,
                other_allele.outcome = NA,
                beta.exposure = NA,
                beta.outcome = NA,
                eaf.exposure = NA,
                eaf.outcome = NA,
                remove = NA,
                palindromic = NA,
                ambiguous = NA,
                id.outcome = NA,
                chr.outcome = NA,
                pos.outcome = NA,
                pval.outcome = NA,
                se.outcome = NA,
                outcome = NA,
                mr_keep.outcome = NA,
                pval_origin.outcome = NA,
                chr.exposure = NA,
                samplesize.exposure = NA,
                se.exposure = NA,
                pval.exposure = NA,
                exposure = NA,
                pval = NA,
                mr_keep.exposure = NA,
                pval_origin.exposure = NA,
                id.exposure = NA,
                action = NA,
                mr_keep = NA,
                samplesize.outcome = NA,
                SNP = NA
              )[-1,]

            df_sum <- rbind(df_sum, results)
            df_instr <- rbind(df_instr, harmonised_clumped_final_data_frame)
            result <- list(results = df_sum,
                           instruments = df_instr)
          }
      }
    }

  message(paste(exposure_id, "done"))

  return(result)
}

# Private functions --------------------------------------------------------


validate_instrument_region_arg <- function(instrument_region) {

  if (is.null(instrument_region)) {
    invisible(TRUE)
  }

  # check is list

  # check correct names in list

  # check chromosome is type integer

  # check start and end are both integers

  invisible(TRUE)
}
