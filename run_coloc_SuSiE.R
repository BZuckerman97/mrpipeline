#' Script which will run single-cell MR analyses
#' I have a folder which contains all the outcome GWAS data I want to use
#' I want to create a function which will read in 1 outcome GWAS file
#' Ensure uniform column headings of SNP, beta, se, pval, effect_allele, other_allele, n (where available)
#' Use the format_single_cell_onek1k() to create a formatted exposure data frame
#' harmonise the data using TwoSampleMR::harmonise_data()
#' run the mr analysis using TwoSampleMR::mr(, method_list = c("inverse variance weighted",|"wald ratio"))
#' output the result into a list I can then compile later
#' I want to use purrr to map through all the onek1k files and the outcome gwas files saving each MR analysis for each outcome GWAS

# Load necessary libraries
library(TwoSampleMR)
library(purrr)
library(dplyr)
library(readr)
library(mrpipeline) # Assuming format_single_cell_onek1k is here

#' Run Single-Cell MR Analysis for one exposure and one outcome
#'
#' This function takes a single exposure GWAS file and a single outcome GWAS file,
#' formats them, harmonises them, and performs a Two-Sample MR analysis.
#'
#' @param exposure_file Path to the exposure data file (e.g., a onek1k file).
#' @param outcome_file Path to the outcome GWAS data file.
#'
#' @return A data frame with the MR results from TwoSampleMR::mr(), or NULL if an error occurs.
#'         The result is augmented with exposure and outcome file names.
run_single_mr <- function(exposure_file, outcome_file) {
  tryCatch({
    # 1. Format exposure data using the provided function
    # Assuming format_single_cell_onek1k() reads the file and formats it correctly
    # into the TwoSampleMR exposure data format.
    exposure_dat <- format_single_cell_onek1k(exposure_file)

    # 2. Read and format outcome data
    # This assumes the outcome GWAS has standard headers.
    # If not, you might need to add column mapping logic here.
    # The function `read_outcome_data` from TwoSampleMR is flexible.
    outcome_dat <- read_outcome_data(
      filename = outcome_file,
      sep = "\t", # Assuming tab-separated, adjust if needed
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      n_col = "n" # Optional, will be used if present
    )

    # 3. Harmonise exposure and outcome data
    dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat
    )

    # If no SNPs are left after harmonisation, return NULL
    if (nrow(dat) == 0) {
      warning(paste("No harmonised SNPs for exposure:", basename(exposure_file), "and outcome:", basename(outcome_file)))
      return(NULL)
    }

    # 4. Run MR analysis
    # Use "Wald ratio" for single-SNP instruments, and "Inverse variance weighted" for multiple.
    methods <- if (nrow(dat) == 1) "wald_ratio" else "mr_ivw"
    res <- mr(dat, method_list = methods)

    # 5. Augment results with file info and return
    res %>%
      mutate(
        exposure_file = basename(exposure_file),
        outcome_file = basename(outcome_file)
      )

  }, error = function(e) {
    # Log the error and return NULL for this pair
    message(sprintf("Error processing exposure '%s' and outcome '%s': %s",
                    basename(exposure_file), basename(outcome_file), e$message))
    return(NULL)
  })
}

# --- Main Script ---

# Define paths to your data folders
# Please update these paths to point to your actual data folders.
exposure_dir <- "path/to/your/onek1k/files"
outcome_dir <- "path/to/your/outcome/gwas/files"

# Get lists of all exposure and outcome files
exposure_files <- list.files(exposure_dir, full.names = TRUE, pattern = "\\.tsv\\.gz$") # Adjust pattern as needed
outcome_files <- list.files(outcome_dir, full.names = TRUE, pattern = "\\.tsv\\.gz$") # Adjust pattern as needed

# Create a grid of all combinations of exposure and outcome files
file_grid <- expand.grid(
  exposure_file = exposure_files,
  outcome_file = outcome_files,
  stringsAsFactors = FALSE
)

# Use purrr::map2 to iterate over the file pairs and run the analysis
all_mr_results <- map2(file_grid$exposure_file, file_grid$outcome_file, run_single_mr)

# Combine the list of results into a single data frame
# The `compact()` function removes any NULL elements from the list (from pairs that failed)
final_results_df <- bind_rows(compact(all_mr_results))

# You can now view or save the final results
print(head(final_results_df))

# To save the results to a CSV file:
# write_csv(final_results_df, "all_single_cell_mr_results.csv")