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
run_mr <- function(expoure,
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
                   rsq_thresh = 0.1) {

# Validate arguments -----------------------------------------------------

validate_instrument_region_arg(instrument_region)


  result <- NULL

  print(exposure_id)
  timestamp()
  print("******")

#' TO DELETE
  #syn_code <-
  #  synGet(entity = exposure_id,
  #         downloadLocation = downloadLocation) # Downloading the summary statistics for the protein of interest

  #if (!dir.exists(fs::path_ext_remove(syn_code$path))) {
  #  untar(paste(syn_code$path),
  #        list = F,
  #        exdir = paste(syn_code$cacheDir))
  #}

#' TO DELETE CAN BE READ IN BEFORE THE FUNCTION
  # exposure <-
  #   fread(
  #     paste0(
  #       syn_code$cacheDir,
  #       "/",
  #       gsub(".tar", "", sumstats_info[sumstats_info$Code == exposure_id,]$Docname[1]),
  #       "/",
  #       "discovery_chr",
  #       sumstats_info[sumstats_info$Code == exposure_id,]$chr[1],
  #       "_",
  #       sumstats_info[sumstats_info$Code == exposure_id,]$UKBPPP_ProteinID[1],
  #       ":",
  #       sumstats_info[sumstats_info$Code == exposure_id,]$Panel[1],
  #       ".gz"
  #     )
  #   )

  #' Taken out sumstats_info here
  #' To replace with tidyverse language
  #' To move out of the run_mr function and have a specify window
  exposure <-
    exposure[exposure$GENPOS > instrument_region$start &
              # Selecting the cis region only (here defined as 200kb before or after the protein-encoding region)
              exposure$GENPOS < instrument_region$end, ]

  exposure$P <- 10 ^ -exposure$LOG10P

  exposure <-
      exposure[exposure$P < pval_thresh,]                                                                                     # Selecting "region-wide" significant cis-pQTLs (here defined as P<5e-6)

  #' This needs to be done before the function starts
  #exposure$CHROM <-
    #  ifelse(exposure$CHROM == 23, "X", exposure$CHROM)                                                                     # Renaming 23rd chromosome "X" for consistency between sum stats

    if (is.null(exposure) || nrow(exposure) == 0) {
      warning(paste0("Skipping ", exposure_id))
      warning("No significant cis pQTLs")
    } else {
      # Only selecting the chromosome of interest to speed up stuff downstream from here
      outcome_overlap <- outcome |>
        filter(`#chrom` %in% exposure$chr) |>
        filter(pos %in% exposure$GENPOS)
      # Only selecting the variants that are overlapping between exposure and outcome sum stats

      if (is.null(outcome_overlap) ||
          nrow(outcome_overlap) == 0) {
        warning(paste0("Skipping ", exposure_id))
        warning("No overlap between outcome and protein exposure")
      } else {
        outcome_rsid <-
          outcome_overlap[, c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
        # outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
        outcome_overlap$phen <- paste(outcome_id)
        outcome_overlap$id <-
          paste(
            outcome_overlap$`#chrom`,
            outcome_overlap$pos,
            outcome_overlap$alt,
            outcome_overlap$ref,
            sep = ":"
          )
        outcome_overlap <-
          format_data(
            outcome_overlap,
            type = "outcome",
            phenotype_col = "phen",
            snp_col = "id",
            beta_col = "beta",
            se_col = "sebeta",
            eaf_col = "af_alt",
            effect_allele_col = "alt",
            other_allele_col = "ref",
            pval_col = "pval",
            chr_col = "#chrom",
            pos_col = "pos"
          )

        exposure_overlap <-
          exposure[exposure$GENPOS %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
        exposure_overlap_2 <-
          exposure_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
        exposure_overlap_2$BETA <- exposure_overlap_2$BETA * -1
        exposure_overlap_2$A1FREQ  <- 1 - exposure_overlap_2$A1FREQ
        colnames(exposure_overlap_2)[colnames(exposure_overlap_2) %in% c("ALLELE0", "ALLELE1")] <-
          c("ALLELE1", "ALLELE0")
        exposure_overlap <- rbind(exposure_overlap, exposure_overlap_2)
        exposure_overlap$ID <-
          paste(
            exposure_overlap$CHROM,
            exposure_overlap$GENPOS,
            exposure_overlap$ALLELE1,
            exposure_overlap$ALLELE0,
            sep = ":"
          )
        exposure_overlap$phen <- exposure_id

        exposure_overlap <-
          format_data(
            exposure_overlap,
            type = "exposure",
            phenotype_col = "phen",
            snp_col = "ID",
            beta_col = "BETA",
            se_col = "SE",
            eaf_col = "A1FREQ",
            effect_allele_col = "ALLELE1",
            other_allele_col = "ALLELE0",
            pval_col = "LOG10P",
            chr_col = "CHROM",
            samplesize_col = "N",
            pos_col = "GENPOS",
            log_pval = T
          )

        dat_u <-
          harmonise_data(exposure_dat = exposure_overlap, outcome_dat = outcome_overlap)                                             # This is where the matching happens

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
        #     c("#chrom", "rsids")
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

        colnames(dat_u)[colnames(dat_u) %in% c("SNP", "rsids")] <-
          c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
        dat_u <-
          dat_u[order(dat_u$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
        dat_u <- dat_u[!duplicated(dat_u$SNP),]

        if (nrow(dat_u) == 0) {
          warning(paste0("Skipping ", exposure_id))
          warning("No variants remaining after harmonising")
          return(NULL)
        }

          print("Clumping")
          clump <-
            ld_clump(
              dplyr::tibble(
                rsid = dat_u$SNP,
                pval = dat_u$pval.exposure,
                id = dat_u$id.exposure
              ),
              # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
              plink_bin = genetics.binaRies::get_plink_binary(),
              clump_kb = 10000,
              clump_r2 = rsq_thresh,
              bfile = "LD_ref/g1000_eur"
            )
          dat <- dat_u[dat_u$SNP %in% clump$rsid,]

          # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.

          if (nrow(dat[dat$mr_keep,]) == 0) {
            warning(paste0("Skipping ", exposure_id))
            warning("No variants remaining after clumping")
            results <- NULL
          } else {
            if (nrow(dat) == 1) {
              # If the genetic instrument includes 1 variant, you use the Wald ratio as your method
              results_mr <-
                mr(dat, method_list = c("mr_wald_ratio"))
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
              # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
              ld <-
                ld_matrix(
                  dat$SNP,
                  bfile = "LD_ref/g1000_eur",
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
              # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
              ld <-
                ld_matrix(
                  dat$SNP,
                  bfile = "LD_ref/g1000_eur",
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

# Private fuctions --------------------------------------------------------


validate_instrument_region_arg(instrument_region) {

  if (is.null(instrument_region)) {
    invisible(TRUE)
  }

  # check is list

  # check correct names in list

  # check chromosome is type integer

  # check start and end are both integers

  invisible(TRUE)
}
