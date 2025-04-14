
#' Formats deCODE genetics proteomic data for Mendelian Randomization analysis
#'
#' This function filters deCODE genetics GWAS by the deCODE included variants file,
#' processes the data to clean the column headings and uses TwoSampleMR::format_data() to prepare
#' the exposure dataset before run_mr()
#'
#' @param pqtl_assay String, the deCODE genetics protein assay
#' @param decode_proteomic_gwas_file_path Dataframe, the file path to the file containing deCODE GWAS data
#' @param decode_included_variants_file_path Dataframe, the file path to the file containing the deCODE GWAS data for variants to include with suitable quality control
#' @param x_y_chr_file String, file path to the file containing rsids for X and Y chromosomes
#'
#' @return A list with two elements:
#'   - `exposure`: Formatted exposure data frame (output of TwoSampleMR::format_data).
#' @export
#'
#' @examples
format_decode_pqtl <- function(decode_proteomic_gwas_file_path,
                               decode_included_variants_file_path,
                               pqtl_assay,
                               x_y_chr_file = NULL){



}
