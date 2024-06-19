#' rsidToGRCh37 Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' rsidToGRCh37(rsids = "rs10")
#' rsidToGRCh37(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
rsidToGRCh37 <- function(rsids)
{
  base_url <- "http://172.15.1.20:8000"
  
  rsids <- toJSON(rsids)
  res <- GET(url = paste0(base_url, "/rsidToGRCh37"), query = list(rsids = rsids))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}

##################################################################################

#' rsidToGRCh38 Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' rsidToGRCh38(rsids = "rs10")
#' rsidToGRCh38(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
rsidToGRCh38 <- function(rsids)
{
  base_url <- "http://172.15.1.20:8000"
  
  rsids <- toJSON(rsids)
  res <- GET(url = paste0(base_url, "/rsidToGRCh38"), query = list(rsids = rsids))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}

##################################################################################


#' SNP_To_Allele_1KGP Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' variant_matrix_2 <- data.frame(Chromosome=c("chr1","chr1","chr1","chr1","chr1"),Position=c(10399,10420,10437,10438,10440))
#' SNP_To_Allele_1KGP(variant_matrix = variant_matrix_2)
#' @export
SNP_To_Allele_1KGP <- function(variant_matrix)
{
  base_url <- "http://172.15.1.20:8000"
  
  variant_matrix <- toJSON(variant_matrix)
  res <- GET(url = paste0(base_url, "/SNP_To_Allele_1KGP"), query = list(variant_matrix = variant_matrix))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(fromJSON(aa)))
}

##################################################################################

#' SNP_To_Allele_Gtex Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' variant_matrix <- data.frame(Chromosome=c("chr1","chr1","chr1","chr1","chr1"), Position=c(13526,13550,14451,14469,14470))
#' SNP_To_Allele_Gtex(variant_matrix)
#' @export
SNP_To_Allele_Gtex <- function(variant_matrix)
{
  base_url <- "http://172.15.1.20:8000"
  
  variant_matrix <- toJSON(variant_matrix)
  res <- GET(url = paste0(base_url, "/SNP_To_Allele_Gtex"), query = list(variant_matrix = variant_matrix))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(fromJSON(aa)))
}

##################################################################################
