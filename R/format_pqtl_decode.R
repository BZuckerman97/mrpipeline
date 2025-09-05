#' Formats deCODE genetics proteomic data for Mendelian Randomization analysis
#'
#' This function filters deCODE genetics GWAS by the deCODE included variants file,
#' processes the data to clean the column headings and uses TwoSampleMR::format_data() to prepare
#' the exposure dataset before run_mr()
#'
#' @param decode_proteomic_gwas_file_path Character vector of file path(s) to deCODE GWAS data,
#'   or a single pre-loaded data.frame. If multiple paths are provided, data will be combined.
#' @param decode_included_variants_file_path Character vector of file path(s) to the deCODE
#'   included variants data, or a single pre-loaded data.frame. If multiple paths are provided,
#'   data will be combined.
#' @param pqtl_assay String, name of the deCODE genetics protein assayed
#' @param x_y_chr_file String, optional file path to a file containing rsids for X and Y
#'   chromosomes. This file should be tab-separated with columns: Chromosome, Position, RSID.
#'
#' @return A list with two elements:
#'   - `exposure`: Formatted exposure data frame (output of TwoSampleMR::format_data).
#' @export
#'
#' @examples
format_pqtl_decode <- function(decode_proteomic_gwas_file_path,
                               decode_included_variants_file_path,
                               pqtl_assay,
                               x_y_chr_file = NULL){

  # Read and combine deCODE proteomic GWAS data
  if (is.character(decode_proteomic_gwas_file_path)) {
    stopifnot(all(sapply(decode_proteomic_gwas_file_path, file.exists)))
    decode_raw_data <- decode_proteomic_gwas_file_path |>
      purrr::map(\(path) data.table::fread(path, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(decode_proteomic_gwas_file_path))
    decode_raw_data <- decode_proteomic_gwas_file_path
  }

  # Read and combine deCODE included variants data
  if (is.character(decode_included_variants_file_path)) {
    stopifnot(all(sapply(decode_included_variants_file_path, file.exists)))
    included_variants_df <- decode_included_variants_file_path |>
      purrr::map(\(path) data.table::fread(path, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(decode_included_variants_file_path))
    included_variants_df <- decode_included_variants_file_path
  }

  # Join GWAS data with included variants (which contains effectAlleleFreq)
  # Assuming 'Name' is the common SNP identifier column (e.g., rsID)
  decode_filtered <- decode_raw_data |>
    dplyr::inner_join(included_variants_df, by = "Name") #' included variants data frame should only contain Name and EAF

  decode_processed <- decode_filtered |>
    dplyr::mutate(phenotype_col = pqtl_assay) |>  #' Create a phenotype_col
    #' Rename columns
    dplyr::rename(
      rsid = dplyr::all_of("rsids"),
      beta = dplyr::all_of("Beta"),
      sebeta = dplyr::all_of("SE"),
      af_alt = dplyr::all_of("effectAlleleFreq"), # This column comes from included_variants_df
      effect_allele = dplyr::all_of("effectAllele"),
      other_allele = dplyr::all_of("otherAllele"),
      pos = dplyr::all_of("Pos"),
      chr = dplyr::all_of("Chrom"),
      pval = dplyr::all_of("Pval")
    ) |>
  #' Edit chromosome variable to change it from "chr3" to 3
    dplyr::mutate(chr = gsub("chr", "", chr)) |>
  #' Rename 23rd chromosome to X for consistency
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", as.character(chr)))

  #' Handle non-Mendelian chromsosome rsIDs/SNPs
  #' Check if the chromosome is X
  if(!is.null(x_y_chr_file)){
    stopifnot(file.exists(x_y_chr_file))
    if("X" %in% unique(decode_processed$chr)){
      #' Load x_y_rsid
      x_y_info_df <- data.table::fread(x_y_chr_file)
      #' Rename columns to match
      x_y_info_df <- x_y_info_df |>
        dplyr::rename(
          chr_xy = dplyr::all_of("V1"), # Avoid name clash if 'chr' already exists
          pos = dplyr::all_of("V2"),
          rsids_xy = dplyr::all_of("V3") # Use a distinct name for rsids from this file
        )
      #' Merge with decode_corect_variants by position
      #' deCODE data is in build37 how should we handle this?
      decode_processed <- dplyr::left_join(decode_processed,
                                           x_y_info_df |> dplyr::select(dplyr::all_of("pos"), dplyr::all_of("rsids_xy")),
                                           by = "pos")

      #' Update rsid column with rsids from x_y_rsid
      decode_processed <- decode_processed |>
        dplyr::mutate(rsid = dplyr::coalesce(rsids_xy, rsid)) |>
        dplyr::select(-dplyr::all_of("rsids_xy")) # Remove temporary column
    }
  }

  #' Convert to a data frame
  decode_final_df <- as.data.frame(decode_processed)

  #' Apply format_data()
  result <- TwoSampleMR::format_data(
    decode_final_df,
    type = "exposure",
    header = TRUE,
    phenotype_col = "phenotype_col",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos",
    log_pval = FALSE
  )

# Remove the "exposure." prefix from column names to match run_mr.R expectations
  colnames(result) <- sub("^exposure\\.", "", colnames(result))

  return(list(exposure = result)) # Return as a list as per roxygen docs
  }


#' DECODE_PQTL_FILE_NAME
#'
#' @param unique_id String, unique sequence associated with each protein assayed
#' @param decode_linker_file String or Dataframe, the file path to the file containing the deCDOE metadata
#' @param decode_dir String, the file path to the directory where the deCODE summary statistics are present
#'
#' @returns a list of file paths for the deCODE proteomic summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a linker file and the deCODE data directory
#' linker_file <- "path/to/your/decode_linker.csv"
#' decode_data_dir <- "path/to/decode_summary_stats/"
#' # Get the first unique ID from the linker for demonstration
#' unique_protein_id <- data.table::fread(linker_file)$seqID[1]
#' file_paths <- decode_pqtl_file_name(unique_id = unique_protein_id,
#'                                   decode_linker_file = linker_file,
#'                                   decode_dir = decode_data_dir)
#' print(file_paths)
#' }
decode_pqtl_file_name <- function(unique_id,
                                  decode_linker_file,
                                  decode_dir){
  if (rlang::is_string(decode_linker_file)) {
    stopifnot(file.exists(decode_linker_file))
    decode_linker_file <- data.table::fread(decode_linker_file)
  } else {
    stopifnot(is.data.frame(decode_linker_file))
  }

  # get relevant metadata
  metadata <- decode_linker_file |>
    dplyr::filter(seqID == unique_id)
  stopifnot(identical(nrow(metadata), 1L))
  metadata <- as.list(metadata)

  # Construct file path using the identifier from metadata
  list(
    # Ensure components are single character strings
    decode = file.path(decode_dir, unlist(metadata$identifier)[1]),
    id = unlist(metadata$seqID)[1]
  )
}
