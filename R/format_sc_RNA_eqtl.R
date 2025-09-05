#' Formats single-cell RNA eQTL data frame
#'
#' Purpose of this script of functions are to create the pathway to
#' the eQTL data and then format the eQTL data ready for analysis
#' Make sure list of cell_types are created/called from the mapping file

format_sc_RNA_eqtl <- function(mapping_file,
                               cell_type,
                               eqtl_dir){
  #' Read in eQTL data for individual cell type using mapping_file pathway
  if (is.character(eqtl_dir)) {
    # Read in files using data.table::fread()
    stopifnot(file.exists(eqtl_dir))

    scRNA_data_frame <- scRNA_data_frame |>
      purrr::map(\(x) data.table::fread(x, nThread = parallel::detectCores())) |>
      dplyr::bind_rows()
  } else {
    stopifnot(is.data.frame(scRNA_data_frame))
  }

  #' Rename variables and keep only ones of interest
  scRNA_data_frame_pre_format_data <- scRNA_data_frame |>
    dplyr::mutate(
      se            = abs(`rho correlation coefficient` / qnorm(pvalue/2)),
      beta          = `rho correlation coefficient`,
      pval          = pvalue,
      cell          = `Cell type`,
      phenotype     = `Gene ID`,
      effect_allele = `SNP assessed allele`) |>
    dplyr::select(cell, phenotype, SNP, effect_allele, beta, se, pval) #' Missing columns to ask why excluded are the other_allele_col

  #' Applying format_data()
  scRNA_data_frame_formatted <- TwoSampleMR::format_data(
    dat = scRNA_data_frame_pre_format_data,
    type = "exposure",
    phenotype_col = "phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    pval_col = "pval"
  )

  return(scRNA_data_frame_formatted)
}
