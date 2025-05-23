# Set up ------------------------------------------------------------------

# library(synapser)
library(R.utils)
library(TwoSampleMR)
library(tidyr)
library(dplyr)
library(data.table)
library(ieugwasr)
library(genetics.binaRies)
library(MendelianRandomization)

synLogin(authToken = Sys.getenv("SYNAPSE_TOKEN"))              # This is to log in to the Synapse platform, needed to gain access to the protein sum stats

sumstats_info <-
  fread("olink_protein_map_3k_v1.tsv")                   # Information on the sum stats / protein codes etc.; functions as a linker file

# Download summary statistics; altered Art's code as no Chr 23 data within the summary statistics
# pss <- fread("data/sjogrens_6098_34928.txt")
#
# pss <- pss |>
#   mutate(phen = "pss") |>
#   rename(
#     rsids = SNP,
#     "#chrom" = chr,
#     sebeta = se,
#     af_alt = eaf,
#     alt = effect_allele,
#     ref = other_allele
#   )

pss <- fread("data/h38_sjogrens_6098_34928.tsv.gz")

pss <- pss |>
  mutate(phen = "pss") |>
  rename(
    rsids = SNP,
    "#chrom" = CHR,
    sebeta = SE,
    af_alt = FRQ,
    alt = A2,
    ref = A1,
    beta = BETA,
    pval = P,
    pos = BP
  )

stopifnot(all(
  c(
    "phen",
    "beta",
    "sebeta",
    "af_alt",
    "alt",
    "ref",
    "pval",
    "#chrom",
    "pos"
  ) %in% names(pss)
))

ref_rsid <-
  fread("hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome


# set download location for synapse files
downloadLocation <- Sys.getenv("SYNAPSE_DOWNLOAD_LOCATION")

# MR function -------------------------------------------------------------

#' Performs MR
#'
#' @param synapse_id character. e.g. "syn52361761"
#' @param outcome data frame. Summary statistics for outcome, must be formatted as above
#' @param outcome_id character. Name for outcome e.g. 'sjogrens'
#' @param sumstats_info data frame.
#' @param downloadLocation path to folder where synapse files are downloaded
#' @param ref_rsid data frame
#' @param pval_thresholds number. 5e-6 by default
#'
#' @return List with 2 data frames - MR results and insruments
perform_mr <- function(synapse_id,
                       outcome,
                       outcome_id,
                       sumstats_info,
                       downloadLocation,
                       ref_rsid,
                       pval_thresholds = 5e-6) {
  result <- NULL

  print(synapse_id)
  timestamp()
  print("******")

  syn_code <-
    synGet(entity = synapse_id,
           downloadLocation = downloadLocation) # Downloading the summary statistics for the protein of interest

  if (!dir.exists(fs::path_ext_remove(syn_code$path))) {
    untar(paste(syn_code$path),
          list = F,
          exdir = paste(syn_code$cacheDir))
  }

  chrom_u <-
    fread(
      paste0(
        syn_code$cacheDir,
        "/",
        gsub(".tar", "", sumstats_info[sumstats_info$Code == synapse_id,]$Docname[1]),
        "/",
        "discovery_chr",
        sumstats_info[sumstats_info$Code == synapse_id,]$chr[1],
        "_",
        sumstats_info[sumstats_info$Code == synapse_id,]$UKBPPP_ProteinID[1],
        ":",
        sumstats_info[sumstats_info$Code == synapse_id,]$Panel[1],
        ".gz"
      )
    )

  chrom_u <-
    chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code == synapse_id,]$gene_start[1] - 200000) &
              # Selecting the cis region only (here defined as 200kb before or after the protein-encoding region)
              chrom_u$GENPOS < (sumstats_info[sumstats_info$Code ==
                                                synapse_id,]$gene_end[1] + 200000), ]
  chrom_u$P <- 10 ^ -chrom_u$LOG10P

  for (pval_thresh in pval_thresholds) {
    chrom <-
      chrom_u[chrom_u$P < pval_thresh,]                                                                                     # Selecting "region-wide" significant cis-pQTLs (here defined as P<5e-6)
    chrom$CHROM <-
      ifelse(chrom$CHROM == 23, "X", chrom$CHROM)                                                                     # Renaming 23rd chromosome "X" for consistency between sum stats

    if (is.null(chrom) || nrow(chrom) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]))
      print("No significant cis pQTLs")
    } else {
      outcome_overlap <-
        outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code == synapse_id,]$chr[1],]                              # Only selecting the chromosome of interest to speed up stuff downstream from here
      outcome_overlap <-
        outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]                                                 # Only selecting the variants that are overlapping between exposure and outcome sum stats

      if (is.null(outcome_overlap) ||
          nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]))
        print("No overlap between outcome and protein exposure")
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

        chrom_overlap <-
          chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
        chrom_overlap_2 <-
          chrom_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
        chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
        chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
        colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <-
          c("ALLELE1", "ALLELE0")
        chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
        chrom_overlap$ID <-
          paste(
            chrom_overlap$CHROM,
            chrom_overlap$GENPOS,
            chrom_overlap$ALLELE1,
            chrom_overlap$ALLELE0,
            sep = ":"
          )
        chrom_overlap$phen <-
          sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]
        chrom_overlap <-
          format_data(
            chrom_overlap,
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
        rm(chrom_overlap_2, chrom)
        dat_u <-
          harmonise_data(exposure_dat = chrom_overlap, outcome_dat = outcome_overlap)                                             # This is where the matching happens

        if (sumstats_info[sumstats_info$Code == synapse_id,]$chr[1] == "X") {
          # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
          dat_u <-
            merge(
              dat_u,
              ref_rsid[, c("V1", "V2", "V3")],
              by.x = "pos.exposure",
              by.y = "V2",
              all.x = T
            )
          colnames(dat_u)[colnames(dat_u) %in% c("V1", "V3")] <-
            c("#chrom", "rsids")
        } else {
          dat_u <-
            merge(
              dat_u,
              outcome_rsid,
              by.x = "pos.exposure",
              by.y = "pos",
              all.x = T
            )
        }

        colnames(dat_u)[colnames(dat_u) %in% c("SNP", "rsids")] <-
          c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
        dat_u <-
          dat_u[order(dat_u$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
        dat_u <- dat_u[!duplicated(dat_u$SNP),]

        if (nrow(dat_u) == 0) {
          print(paste0("Skipping ", sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]))
          print("No variants remaining after harmonising")
          break()
        }

        for (rsq_thresh in c(0.1)) {
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
          rm(chrom_overlap,
             outcome_overlap,
             outcome_rsid,
             clump)

          # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.

          if (nrow(dat[dat$mr_keep,]) == 0) {
            print(paste0("Skipping ", sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]))
            print("No variants remaining after clumping")
            results <- NULL
          } else {
            if (nrow(dat) == 1) {
              # If the genetic instrument includes 1 variant, you use the Wald ratio as your method
              results_mr <-
                mr(dat, method_list = c("mr_wald_ratio"))
              results <-
                data.frame(
                  exp = sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1],
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
                  exp = sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1],
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
                  exp = sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1],
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
                  exp = sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1],
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
                  exp = sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1],
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
              rm(results1, results2, results3)

            }
          }

          if (is.null(results) || nrow(results) == 0) {
            print(paste0("Skipping ", sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1]))
            print("No results returned from MR analysis")
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
            rm(dat,
               results,
               ld,
               dat2,
               output_mr_ivw_corr,
               output_mr_egger_corr)

            result <- list(results = df_sum,
                           instruments = df_instr)
          }
        }
      }
    }
  }

  print(paste(sumstats_info[sumstats_info$Code == synapse_id,]$Assay[1], "done"))
  # unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)                    # ATTENTON !  This step deletes the map you downloaded your sum stats in - double-check that this doesn't accidentally delete important files from your desktop/server or so

  return(result)
}

system.time(withr::with_options(list(future.globals.maxSize = 2500000000),
                                {
                                  result <- sumstats_info$Code |>
                                    purrr::set_names(\(x) sumstats_info[sumstats_info$Code == x, ]$Assay[452]) |>
                                    furrr::future_map(\(x) {
                                      synLogin(authToken = Sys.getenv("SYNAPSE_TOKEN"))
                                      mr_results <- perform_mr(
                                        synapse_id = x,
                                        outcome = pss,
                                        outcome_id = "pss",
                                        sumstats_info = sumstats_info,
                                        downloadLocation = downloadLocation,
                                        ref_rsid = ref_rsid,
                                        pval_thresholds = 5e-6
                                      )
                                      cat(
                                        paste0(x, " ", Sys.time(), "\n"),
                                        file = file.path("output", paste0(Sys.Date(), "-mr.log")),
                                        append = TRUE
                                      )
                                      mr_results
                                    },
                                    .progress = TRUE)
                                }))

result
saveRDS(result, "output/240605_pss_mr_results_raw_subset_bz.rds")

# process
no_result <- result %>%
  keep(is.null) %>%
  names()

mr_results <- result %>%
  compact() %>%
  map(\(x) x$results) %>%
  bind_rows()

mr_instruments <- result %>%
  compact() %>%
  map(\(x) x$instruments) %>%
  bind_rows()

fwrite(mr_results,
       file = "output/mr_results_pss_debug.csv")

fwrite(mr_instruments,
       file = "output/mr_instruments_pss_debug.csv")

fwrite(data.frame(protein = no_result),
       file = "output/no_mr_results_pss_debug.csv")
