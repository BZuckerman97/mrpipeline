% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_pqtl_ukbppp.R
\name{format_pqtl_ukbppp}
\alias{format_pqtl_ukbppp}
\title{Formats pQTL data from UKB-PPP for MR analysis.}
\usage{
format_pqtl_ukbppp(ukbppp, ukbppp_rsid, pqtl_assay, x_y_chr_file = NULL)
}
\arguments{
\item{ukbppp}{Dataframe, file path to the file containing the ukbppp GWAS data}

\item{ukbppp_rsid}{Dataframe, file path to the file containing the ukbppp GWAS data rsids}

\item{pqtl_assay}{String, of the ukbppp protein assay}

\item{x_y_chr_file}{String, file path to the file containing rsids for X and Y chromosomes}
}
\value{
A list with two elements:
\itemize{
\item \code{exposure}: Formatted exposure data frame (output of TwoSampleMR::format_data).
}
}
\description{
This function processes pQTL summary statistics from UKB-PPP,
matches them with rsID information, handles non-Mendelian
chromosomes, standardizes column names, and formats the data
for use with the TwoSampleMR package.
}
\examples{
# See the test script for example usage.
}
