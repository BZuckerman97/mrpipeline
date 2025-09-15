#' Format single-cell RNA eQTL data for analysis
#'
#' This function processes single-cell eQTL data from the OneK1K project for a
#' specific cell type. It uses a mapping data frame to find the path to the
#' relevant eQTL summary statistics file, reads the data, and formats it
#' for use with the `TwoSampleMR` package.
#'
#' @param onek1k_mapping A data frame containing mapping information for OneK1K
#'   eQTL data. Must contain 'cell_type' and 'path_to_eqtl_file' columns.
#' @param onek1k_cell_type Character string. The name of the cell type to be
#'   processed (e.g., "CD4 Naive"). This must correspond to an entry in the
#'   'cell_type' column of `onek1k_mapping`.
#'
#' @return A data frame of the formatted eQTL data, ready for use as an
#'   "exposure" dataset in `TwoSampleMR`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy mapping data frame
#' mapping_df <- data.frame(
#'   cell_type = "CD4 Naive",
#'   path_to_eqtl_file = "path/to/your/eqtl_data.tsv.gz"
#' )
#' # Run the function
#' formatted_eqtls <- format_single_cell_onek1k(
#'   onek1k_mapping = mapping_df,
#'   onek1k_cell_type = "CD4 Naive"
#' )
#' }
format_single_cell_onek1k <- function(onek1k_mapping,
                                      onek1k_cell_type){

  # Input validation
  stopifnot(is.data.frame(onek1k_mapping))
  stopifnot(is.character(onek1k_cell_type) && length(onek1k_cell_type) == 1)
  stopifnot(all(c("cell_type", "path_to_eqtl_file") %in% names(onek1k_mapping)))

  # Filter mapping file for the specified cell type
  cell_type_info <- onek1k_mapping |>
    dplyr::filter(cell_type == onek1k_cell_type)

  if (nrow(cell_type_info) == 0) {
    stop("The specified `onek1k_cell_type` was not found in the mapping data frame.")
  }
  if (nrow(cell_type_info) > 1) {
    warning("Multiple entries found for the specified `onek1k_cell_type`. Using the first one.")
    cell_type_info <- cell_type_info[1, ]
  }

  eqtl_file_path <- cell_type_info$path_to_eqtl_file

  # Read eQTL data
  eqtl_data <- data.table::fread(eqtl_file_path, nThread = parallel::detectCores())

  # Rename variables, calculate standard error, and select columns
  eqtl_data_prepped <- eqtl_data |>
    dplyr::rename(
      SNP = RSID,
      beta = SPEARMANS_RHO,
      pval = P_VALUE ,
      phenotype = GENE,
      effect_allele = A1,
      other_allele = A2
    ) |>
    dplyr::mutate(
      # Calculate standard error from beta and p-value, do I need to do anything to handle p=0 or p=1 ie case_when(pval == 0 | pval == 1 ~ 0, TRUE ~
      se = dplyr::mutate(abs(beta / qnorm(pval / 2)))) |>
    dplyr::select(cell, phenotype, SNP, effect_allele, other_allele, beta, se, pval) #' Missing columns to ask why excluded are the other_allele_col

  # Format data using TwoSampleMR
  # Do I need to add in the other_allele column?
  formatted_data <- TwoSampleMR::format_data(
    dat = eqtl_data_prepped,
    type = "exposure",
    phenotype_col = "phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    pval_col = "pval"
  )

  return(formatted_data)
}

