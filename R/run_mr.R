#' Performs MR
#'
#' @param exposure_id character. e.g. "syn52361761" -> we should change this to the protein/gene name
#' @param expoure data frame. Summary statistics for exposure, must include exposure_id column ie gene name which we use as above?
#' @param outcome data frame. Summary statistics for outcome, must be formatted as above
#' @param outcome_id character. Name for outcome e.g. 'SjD'
#' @param pval_thresh number. 5e-6 by default
#' @param rsq_thresh R square clumping threshold
#' @param instrument_region list of the chromosome position, gene start and gene end for each gene of interest this needs to link with the original mapping file
#' @return List with 2 data frames - MR results and instruments
run_mr <- function(exposure,
                   exposure_id,
                   outcome,
                   outcome_id,
                   # sumstats_info,
                   # downloadLocation,
                   # ref_rsid,
                   instrument_region = list(chromosome = 1L,
                                            start = 1L,
                                            end = 200L),
                   pval_thresh = 5e-6,
                   rsq_thresh = 0.1,
                   bfile) {

# Validate arguments -----------------------------------------------------

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
    dplyr::filter(POS19 > instrument_region$start & POS19 < instrument_region$end) |> # Selecting the cis region only (here defined as 200kb before or after the protein-encoding region)
    dplyr::mutate(P = 10 ^ -LOG10P) |>
    dplyr::filter(P < pval_thresh)                                                                                   # Selecting "region-wide" significant cis-pQTLs (here defined as P<5e-6)

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
        dplyr::filter(chrom %in% exposure$CHROM) |> #' CHROM needs to change to chr.exposure
        dplyr::filter(pos %in% exposure$POS19)
      # Only selecting the variants that are overlapping between exposure and outcome sum stats

      if (is.null(outcome_overlap) ||
          nrow(outcome_overlap) == 0) {
        warning(paste0("Skipping ", exposure_id))
        warning("No overlap between outcome and protein exposure")
      } else {


# Reformat outcome df -----------------------------------------------------

        outcome_rsid <- outcome_overlap |>
          dplyr::select(chrom, pos, rsids)

        outcome_overlap <- outcome_overlap |>
          dplyr::mutate(phenotype = paste(outcome_id)) |>
          dplyr::mutate(id = paste(chrom, pos, alt, ref, sep = ":"))

#' This will be taken out as should happen prior to run_mr() function used
        outcome_overlap <-
          TwoSampleMR::format_data(
            outcome_overlap,
            type = "outcome",
            phenotype_col = "phenotype",
            snp_col = "rsids",
            beta_col = "beta",
            se_col = "sebeta",
            eaf_col = "af_alt",
            effect_allele_col = "alt",
            other_allele_col = "ref",
            pval_col = "pval",
            chr_col = "chrom",
            pos_col = "pos"
          )

# Reformat exposure df ----------------------------------------------------

        # Don't think we need this but have included until we make that decision
        exposure_overlap <- exposure |>
          dplyr::filter(POS19 %in% outcome_overlap$pos.outcome)

        exposure_overlap2 <- exposure_overlap |>
          dplyr::mutate(BETA = BETA* -1) |>
          dplyr::mutate(A1FREQ = 1- A1FREQ) |>
          dplyr::mutate(ALLELEX = ALLELE0,
                        ALLELE0 = ALLELE1,
                        ALLELE1 = ALLELEX)

        exposure_overlap <- dplyr::bind_rows(exposure_overlap, exposure_overlap2)
        exposure_overlap <- exposure_overlap |>
          dplyr::mutate(ID = paste(CHROM, POS19, ALLELE1, ALLELE0, sep = ":")) |>
          dplyr::mutate(phenotype = exposure_id)

        #' This will be moved outside of the function
        exposure_overlap <-
          TwoSampleMR::format_data(
            exposure_overlap,
            type = "exposure",
            phenotype_col = "phenotype",
            snp_col = "rsid",
            beta_col = "BETA",
            se_col = "SE",
            eaf_col = "A1FREQ",
            effect_allele_col = "ALLELE1",
            other_allele_col = "ALLELE0",
            pval_col = "P",
            chr_col = "CHROM",
            samplesize_col = "N",
            pos_col = "POS19",
            log_pval = T
          )

# Harmonise ---------------------------------------------------------------

        dat_u <-
          TwoSampleMR::harmonise_data(exposure_dat = exposure_overlap, outcome_dat = outcome_overlap)                                             # This is where the matching happens

    #' To move pre-function for handling X chromosomes
        # if (sumstats_info[sumstats_info$Code == exposure_id,]$chr[1] == "X") {
        #   # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
        #   dat_u <-
        #     merge(
        #       dat_u,
        #       ref_rsid[, c("V1", "V2", "V3")],
        #       by.x = "pos.exposure",
        #       by.y = "V2",
        #       all.x = T
        #     )
        #   colnames(dat_u)[colnames(dat_u) %in% c("V1", "V3")] <-
        #     c("chrom", "rsids")
        # } else {
        #   dat_u <-
        #     merge(
        #       dat_u,
        #       outcome_rsid,
        #       by.x = "pos.exposure",
        #       by.y = "pos",
        #       all.x = T
        #     )
        # }

        #' Don't need to do this unless we use IDs instead of rsIDs
        #colnames(dat_u)[colnames(dat_u) %in% c("SNP", "rsids")] <-
        #  c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column

        dat_u <- dat_u |>
          dplyr::arrange(pval.exposure)
                                                                                             # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
        dat_u <- dat_u |>
          dplyr::filter(!duplicated(SNP))

        if (nrow(dat_u) == 0) {
          warning(paste0("Skipping ", exposure_id))
          warning("No variants remaining after harmonising")
          return(NULL)
        }

# Clump -------------------------------------------------------------------
        print("Clumping")
          clump <-
            ieugwasr::ld_clump(
              dplyr::tibble(
                rsid = dat_u$SNP,
                pval = dat_u$pval.exposure,
                id = dat_u$id.exposure
              ),
              # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
              plink_bin = genetics.binaRies::get_plink_binary(),
              clump_kb = 10000,
              clump_r2 = rsq_thresh,
              bfile = bfile
            )
          dat <- dat_u |>
            dplyr::filter(dat_u$SNP %in% clump$rsid)

          # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.

# Perform MR --------------------------------------------------------------

          if (nrow(dat[dat$mr_keep,]) == 0) {
            warning(paste0("Skipping ", exposure_id))
            warning("No variants remaining after clumping")
            results <- NULL
          } else {

## Wald ratio --------------------------------------------------------------

            if (nrow(dat) == 1) {
              # If the genetic instrument includes 1 variant, you use the Wald ratio as your method
              results_mr <-
                TwoSampleMR::mr(dat, method_list = c("mr_wald_ratio"))
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
            } else if (nrow(dat) == 2) {


## 2 IVs - IVW not Egger ---------------------------------------------------

              # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
              ld <-
                ieugwasr::ld_matrix(
                  dat$SNP,
                  bfile = bfile,
                  plink_bin = genetics.binaRies::get_plink_binary()
                )
              dat2 <-
                MendelianRandomization::mr_input(
                  bx = dat$beta.exposure,
                  bxse = dat$se.exposure,
                  by = dat$beta.outcome,
                  byse = dat$se.outcome,
                  correlation = ld
                )
              output_mr_ivw_corr <-
                MendelianRandomization::mr_ivw(dat2, correl = TRUE)
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
            } else {

# >2 IVs - IVW, Egger -----------------------------------------------------

               # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
              ld <-
                ieugwasr::ld_matrix(
                  dat$SNP,
                  bfile = bfile,
                  plink_bin = genetics.binaRies::get_plink_binary()
                )
              dat2 <-
                MendelianRandomization::mr_input(
                  bx = dat$beta.exposure,
                  bxse = dat$se.exposure,
                  by = dat$beta.outcome,
                  byse = dat$se.outcome,
                  correlation = ld
                )
              output_mr_ivw_corr <-
                MendelianRandomization::mr_ivw(dat2, correl = TRUE)
              output_mr_egger_corr <-
                MendelianRandomization::mr_egger(dat2, correl = TRUE)

# Return results ----------------------------------------------------------

              results1 <-
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
              results2 <-
                data.frame(
                  exp = exposure_id,
                  outc = paste(outcome_id),
                  pvalthreshold = pval_thresh,
                  rsqthreshold = rsq_thresh,
                  nsnp = output_mr_egger_corr@SNPs,
                  method = "Egger (correlation inc)",
                  b = output_mr_egger_corr@Estimate,
                  se = output_mr_egger_corr@StdError.Est,
                  pval = output_mr_egger_corr@Pvalue.Est
                )
              results3 <-
                data.frame(
                  exp = exposure_id,
                  outc = paste(outcome_id),
                  pvalthreshold = pval_thresh,
                  rsqthreshold = rsq_thresh,
                  nsnp = output_mr_egger_corr@SNPs,
                  method = "Egger intercept (correlation inc)",
                  b = output_mr_egger_corr@Intercept,
                  se = output_mr_egger_corr@StdError.Int,
                  pval = output_mr_egger_corr@Pvalue.Int
                )
              results <- rbind(results1, results2, results3)

            }
          }

          if (is.null(results) || nrow(results) == 0) {
            warning(paste0("Skipping ", exposure_id))
            warning("No results returned from MR analysis")
          } else {
            df_sum <-
              data.frame(
                exp = NA,
                outc = NA,
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
            df_instr <- rbind(df_instr, dat[,-c(34)])

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
