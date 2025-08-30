# Library
library(data.table)
library(dplyr)
library(mrpipeline)
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(purrr)
library(furrr)

#' Date

today <- Sys.Date()

#' Set up
OUTCOME_ID <- "Seropositive RhA meta-analysis"
EXPOSURE_MAPPING_FILE_DIR <- "olink_protein_map_3k_v1.tsv"
OLINK_DIR <- "~/Shared3/rmgpbzu/ukbppp/"
OLINK_RSID_DIR <- "~/Shared3/rmgpbzu/ukbppp_metadata/"
PVAL_THRESH <- 5e-6
CIS_WINDOW_KB <- 2e+5
CLUMP_RSQ_THRESH <- 0.1
LD_REF <- "LD_ref/g1000_eur"
X_Y_CHR_FILE <- "hg38_common_chrpos_X.txt"

outcome_raw <- fread("~/Shared2/rmgpbzu/pqtl_mr/sjogren_mr/proteomic_seronegative_diseases/seroneg_rha_project/meta_analysis_gwas_files/seronegative_ra_meta_analysis_mapped_final_300725.txt")

outcome_clean <- outcome_raw %>%
  select(-sample_size, -V1) %>%
  mutate(N_cases = 17221+18019) %>%
  mutate(N_controls = 74823 + 991604) %>%
  rename(rsids = SNP) %>%
  mutate(n = N_cases)

ukbppp_protein_seropos_ra_mr <- function(
    protein_synapse_id, # Scalar: The synapse_id for the current protein (exposure)
    protein_name,   # Scalar: The assay/protein name for the current protein (exposure)
    exposure_meta_linker_file, # Data frame: The complete exposure_mapping_file, used as linker for protein metadata
    olink_dir, # Character: Path to UKBPPP base directory for GWAS files (exposure)
    olink_rsid_dir, # Character: Path to UKBPPP RSID GWAS files (exposure)
    outcome, # Data frame: Pre-formatted outcome data (e.g., AMD GWAS)
    outcome_id, # Character: Name of the outcome phenotype (e.g., OUTCOME_ID)
    exposure_pval_filter_threshold, # Numeric: P-value threshold for selecting exposure instruments
    exposure_cis_window_kb, # Numeric: Kilobase window around gene for cis-instrument selection
    clumping_r2_thresh, # Numeric: R-squared threshold for clumping instruments
    ld_ref_path_prefix, # Character: Path prefix for LD reference panel
    xy_chr_pos_filepath # Character: Path to X Y chromosome position file (for formatting pQTL data)
){
  # Filter the main exposure metadata linker file to get the row for the current protein.
  single_protein_linker_data <- exposure_meta_linker_file %>%
    filter(Code == protein_synapse_id)

  if (nrow(single_protein_linker_data) != 1) {
    stop(paste0("Error: Expected 1 row in linker file for protein_synapse_id '",
                protein_synapse_id, "', but found ", nrow(single_protein_linker_data),
                ". Check 'Code' column and uniqueness of IDs in '", EXPOSURE_MAPPING_FILE_DIR, "'."))
  }

  # Construct file paths for the UKB-PPP pQTL summary statistics for the current protein.
  exposure_file_paths <- mrpipeline::ukbppp_pqtl_file_name(
    synapse_id = protein_synapse_id, # Pass the scalar synapse ID for the current protein
    olink_linker_file = single_protein_linker_data, # Pass the filtered 1-row linker data
    olink_dir = olink_dir, # Use function argument
    olink_rsid_dir = olink_rsid_dir # Use function argument
  )

  # Read in and format the pQTL (exposure) data using the constructed file paths.
  exposure_format <- mrpipeline::format_pqtl_ukbppp(ukbppp = exposure_file_paths$ukbppp,
                                                    ukbppp_rsid = exposure_file_paths$ukbppp_rsid,
                                                    pqtl_assay = protein_name,
                                                    x_y_chr_file = xy_chr_pos_filepath) # Use function argument

  # Convert the formatted exposure data to a standard data frame.
  exposure_format <- as.data.frame(exposure_format)

  # Run the Mendelian Randomization analysis using the formatted exposure and provided outcome data.
  run_mr_result <- mrpipeline::run_mr(exposure = exposure_format, # Use the formatted exposure data
                                      exposure_id = protein_name, # Use the protein_name argument
                                      outcome = outcome,
                                      outcome_id = outcome_id,
                                      # Define the genomic region for instrument selection (cis-pQTLs).
                                      instrument_region = list(chromosome = single_protein_linker_data$chr,
                                                               start = single_protein_linker_data$gene_start,
                                                               end = single_protein_linker_data$gene_end),
                                      window = exposure_cis_window_kb, # Kilobase window around the gene for cis-pQTLs.
                                      pval_thresh = exposure_pval_filter_threshold,
                                      rsq_thresh = clumping_r2_thresh,
                                      bfile = ld_ref_path_prefix)

  return(run_mr_result) # Explicitly return the result
}

exposure_mapping_file <- fread(EXPOSURE_MAPPING_FILE_DIR) # Load the mapping file

proteomic_ukbppp_seropos_ra_mr <- purrr::map(
  .x = seq_len(nrow(exposure_mapping_file)),
  .f = ~ {
    current_protein_synapse_id <- exposure_mapping_file$Code[.x]
    current_protein_name <- exposure_mapping_file$Assay[.x]

    ukbppp_protein_seropos_ra_mr(
      protein_synapse_id = current_protein_synapse_id,
      protein_name = current_protein_name,
      exposure_meta_linker_file = exposure_mapping_file, # Pass the loaded mapping file
      olink_dir = OLINK_DIR, # Use global constant
      olink_rsid_dir = OLINK_RSID_DIR, # Use global constant
      outcome = outcome_clean, # Pass the cleaned outcome data
      outcome_id = OUTCOME_ID, # Use global constant
      exposure_pval_filter_threshold = PVAL_THRESH, # Use global constant
      exposure_cis_window_kb = CIS_WINDOW_KB, # Use global constant
      clumping_r2_thresh = CLUMP_RSQ_THRESH, # Use global constant
      ld_ref_path_prefix = LD_REF, # Use global constant
      xy_chr_pos_filepath = X_Y_CHR_FILE # Use global constant
    )
  }
)

saveRDS(proteomic_ukbppp_seropos_ra_mr, paste0("proteomic_seropositive_diseases/seroneg_rha_project/output/proteomic_rha_mr_analysis/proteomic_ukbppp_seronegative_ra_meta_analysis_gwas", today, ".rds"))

# process
no_result <- proteomic_ukbppp_seropos_ra_mr %>%
  keep(is.null) %>%
  names()

mr_results <- proteomic_ukbppp_seropos_ra_mr %>%
  compact() %>%
  map(\(x) x$results) %>%
  bind_rows()

mr_instruments <- proteomic_ukbppp_seropos_ra_mr %>%
  compact() %>%
  map(\(x) x$instruments) %>%
  bind_rows()

fwrite(mr_results,
       file = paste0("proteomic_seropositive_diseases/seroneg_rha_project/output/proteomic_rha_mr_analysis/proteomic_ukbppp_seroneg_ra_mr_results_full_pqtl_analysis_", today, ".csv"))

fwrite(mr_instruments,
       file = paste0("proteomic_seropositive_diseases/seroneg_rha_project/output/proteomic_rha_mr_analysis/proteomic_ukbppp_seroneg_ra_mr_instruments_full_pqtl_analysis__", today, ".csv"))

fwrite(data.frame(protein = no_result),
       file = paste0("proteomic_seropositive_diseases/seroneg_rha_project/output/proteomic_rha_mr_analysis/proteomic_ukbppp_seroneg_ra_no_mr_results_full_pqtl_analysis_", today, ".csv"))

#Remember to add the new number of cases and controls in at the end
