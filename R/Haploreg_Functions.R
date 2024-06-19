#' snpToAF Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' snpToAF(rsids = "rs10")
#' snpToAF(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
snpToAF <- function(rsids)
{
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("snpToAF(rsids=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
}


#' snpToCHROMHMM_15STATE Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsid genotypic variant as an rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToCHROMHMM_15STATE(rsid <- "rs10")
#' @export
snpToCHROMHMM_15STATE <- function(rsid)
{
  ex = paste0("snpToCHROMHMM_15STATE(rsid='" , rsid , "')")
  return(RSeval(expr = ex , c = connec))
}


#' snpToCHROMHMM_25STATE Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsid genotypic variant as an rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToCHROMHMM_25STATE(rsid <- "rs10");
#' @export
snpToCHROMHMM_25STATE <- function(rsid)
{
  ex = paste0("snpToCHROMHMM_25STATE(rsid='" , rsid , "')")
  return(RSeval(expr = ex , c = connec))
}



#' snpToGENCODE Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToGENCODE(rsids <- "rs10");
#' snpToGENCODE(rsids <- c("rs10","rs3","rs4","rs5"));
#' @export
snpToGENCODE <- function(rsids)
{
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("snpToGENCODE(rsids=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
}



#' snpToHISTONEPEAKS Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsid genotypic variant as an rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToHISTONEPEAKS(rsid <- "rs10");
#' @export
snpToHISTONEPEAKS <- function(rsid)
{
  ex = paste0("snpToHISTONEPEAKS(rsid='" , rsid , "')")
  return(RSeval(expr = ex , c = connec))
}



#' snpToInfo Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToInfo(rsids = "rs10")
#' snpToInfo(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
snpToInfo <- function(rsids)
{
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("snpToInfo(rsids=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
}




#' snpToMotif Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToMotif(rsid = "rs10")
#' @export
snpToMotif <- function(rsid)
{
  ex = paste0("snpToMotif(rsid='" , rsid , "')")
  return(RSeval(expr = ex , c = connec))
}





#' snpToREFSEQ Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToREFSEQ(rsids = "rs10")
#' snpToREFSEQ(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
snpToREFSEQ <- function(rsids)
{
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("snpToREFSEQ(rsids=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
}
