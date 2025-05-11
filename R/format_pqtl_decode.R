#' Formats deCODE genetics proteomic data for Mendelian Randomization analysis
#'
#' This function filters deCODE genetics GWAS by the deCODE included variants file,
#' processes the data to clean the column headings and uses TwoSampleMR::format_data() to prepare
#' the exposure dataset before run_mr()
#'
#' @param decode_proteomic_gwas_file_path Dataframe, the file path to the file containing deCODE GWAS data
#' @param decode_included_variants_file_path Dataframe, the file path to the file containing the deCODE GWAS data for variants to include with suitable quality control
#' @param pqtl_assay String, name of the deCODE genetics protein assayed
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

  #' Read from file path
  if(is.character(decode_proteomic_gwas_file_path)){
    #' Read in deCODE proteomic files using data.table::fread()
    stopifnot(file.exists(decode_proteomic_gwas_file_path))
    decode <- decode |>
      purr::map(\(x)data.table::fread(x, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  }else{
    stopifnot(is.data.frame(decode_proteomic_gwas_file_path))
  }

  #' Read from file path
  if(is.character(decode_included_variants_file_path)){
    #' Read in deCODE included variants file
    stopifnot(file.exists(decode_variants_file_path))
    decode_included_variants <- decode_included_variants |>
      purr::map(\(x)data.table::fread(x, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  }else{
    stopifnot(is.data.frame(decode_included_variants))
  }

  #' Join decode_included_variants_file so that only the correct variants are included

  decode_correct_variants <- decode |>
    dplyr::inner_join(decode_included_variants |>
                 dplyr::select(Name, effectAlleleFreq), by = "Name")

  decode_correct_variants <- decode_correct_variants |>
    dplyr::mutate(phenotype_col = pqtl_assay) |>  #' Create a phenotype_col
    #' Rename columns
    dpylr::rename(
      rsid = dplyr::all_of("rsids"),
      beta = dplyr::all_of("Beta"),
      sebeta = dplyr::all_of("SE"),
      af_alt = dplyr::all_of("effectAlleleFreq"),
      effect_allele = dplyr::all_of("effectAllele"),
      other_allele = dplyr::all_of("otherAllele"),
      pos = dplyr::all_of("Pos"),
      chr = dplyr::all_of("Chrom"),
      pval = dplyr::all_of("Pval") |>
  #' Edit chromosome variable to change it from "chr3" to 3
    dplyr::mutate(chr = gsub("chr", "", chr)) |>
  #' Rename 23rd chromosome to X for consistency
    dplyr::mutate(chr = dplyr::if_else(chr == "23", "X", as.character(chr)))

  #' Handle non-Mendelian chromsosome rsIDs/SNPs
  #' Check if the chromosome is X
  if(!is.null(x_y_chr_file)){
    stopifnot(file.exists(x_y_chr_file))
    if("X" %in% unique(decode_correct_variants$chr)){
      #' Load x_y_rsid
      x_y_rsid <- data.table::fread(x_y_chr_file)
      #' Rename columns to match
      x_y_rsid <- x_y_rsid |>
        dplyr::rename(
          pos = dplyr::all_of("V2"),
          rsids = dplyr::all_of("V3"),
          chr = dplyr::all_of("V1")
        )
      #' Merge with decode_corect_variants by position
      #' deCODE data is in build37 how should we handle this?
      decode_correct_variants <- dplyr::left_join(decode_correct_variants, x_y_rsid[,c("pos","rsids")], by = "pos")

      #' Update rsid column with rsids from x_y_rsid
      decode_correct_variants <- decode_correct_variants |>
        dplyr::mutate(rsid = dplyr::coalesce(rsids, rsid)) |>
        dplyr::select(-rsids)
    }
  }

  #' Filter for the correct chromosome - can't remember why I had this in the original script I wrote
  #' Vague memory of having instances where there were multiple chromosomes in one file, but maybe this
  #' was a fever dream

  #' Convert to a data frame

  decode_correct_variants <- as.data.frame(decode_correct_variants)
  #' Apply format_data()

  result <- TwoSampleMR::format_data(
    decode_correct_variants,
    type = "exposure",
    header = "TRUE",
    phenotype_col = "phenotpye_col",
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

  return(result)
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

  # result
  list(
    decode = file.path(decode_dir,
                       paste0(decode_linker_file$identifier[decode_linker_file$seqID == unique_id])),
    id = metadata$seqID
  )
}
