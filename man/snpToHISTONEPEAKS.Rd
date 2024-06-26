% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Haploreg_Functions.R
\name{snpToHISTONEPEAKS}
\alias{snpToHISTONEPEAKS}
\title{snpToHISTONEPEAKS Function}
\usage{
snpToHISTONEPEAKS(rsid)
}
\arguments{
\item{rsid}{genotypic variant as rsid}
}
\value{
Summary Statistics of the SNP,Gene and Tissue
}
\description{
Retrieve HISTONEPEAKS information for a set of SNPs from the HAPLOREG database.
}
\examples{
snpToHISTONEPEAKS(rsid <- "rs3");
}
