% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_pqtl_ukbppp.R
\name{ukbppp_pqtl_file_name}
\alias{ukbppp_pqtl_file_name}
\title{UKBPPP_PQTL_FILE_NAME}
\usage{
ukbppp_pqtl_file_name(synapse_id, olink_linker_file, olink_dir, olink_rsid_dir)
}
\arguments{
\item{synapse_id}{String, synapse id to access from olink linker file}

\item{olink_linker_file}{String or Dataframe, file containing the olink linker file, or the dataframe of the linker file}

\item{olink_dir}{String, directory of the olink files}

\item{olink_rsid_dir}{String, directory of the olink rsid files}
}
\value{
a list of the 2 filepaths, one for the ukbppp_pqtl data and the other is
the corresponding rsID metadata file, as well as the name of the assay
}
\description{
UKBPPP_PQTL_FILE_NAME
}
\examples{
\dontrun{
Assuming you have a linker file and the deCODE data directory
synapse_id
olink_linker_file <- "path/to/your/olink_linker_file.csv"
olink_dir <- "path/to/olink_dir"
olink_rsid_dir <- "path/to/olink_rsid_dir"
file_paths <- ukbppp_pqtl_file_name(synapse_id = synapse_id,
                                  olink_linker_file = olink_linker_file,
                                  olink_dir = olink_dir,
                                  olink_rsid_dir = olink_rsid_dir)
print(file_paths)
}
}
