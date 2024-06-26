% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Haploreg_Functions.R
\name{snpToCHROMHMM_15STATE}
\alias{snpToCHROMHMM_15STATE}
\title{snpToCHROMHMM_15STATE Function}
\usage{
snpToCHROMHMM_15STATE(rsid)
}
\arguments{
\item{rsid}{A single variant with rsid}
}
\value{
Summary Statistics of the SNP,Gene and Tissue
}
\description{
Retrieve the CHROMHMM 15-state annotations for a set of SNPs from the HAPLOREG database.
}
\examples{
snpToCHROMHMM_15STATE(rsid <- "rs10")
}
