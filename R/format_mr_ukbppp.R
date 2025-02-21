#' Formats data for using UKB-PPP data
#'
#' @param olink_linker_file
#' @param synapse_id
#' @param outcome
#' @param outcome_id
#' @param ref_rsid
#'
#' @returns
#' @export
#'
#' @examples
format_mr_ukbppp <- function(olink_linker_file,
                             synapse_id,
                             olink_dir,
                             outcome,
                             outcome_id,
                             ref_rsid){


}

synapse_pqtl_file_name <- function(olink_linker_file,
                              synapse_id,
                              olink_dir){
  file.path(olink_dir,
      paste0(gsub(".tar", "", olink_linker_file[olink_linker_file$Code == synapse_id,]$Docname[1]),
        "/",
        "discovery_chr",
        olink_linker_file[olink_linker_file$Code == synapse_id,]$chr[1],
        "_",
        olink_linker_file[olink_linker_file$Code == synapse_id,]$UKBPPP_ProteinID[1],
        ":",
        olink_linker_file[olink_linker_file$Code == synapse_id,]$Panel[1],
        ".gz"
      )
  )
}

synapse_rsid_file_name <- function(olink_linker_file,
                                   metadata_dir){
  file.path(metadata_dir,
            paste0(olink_linker_file$))
}
