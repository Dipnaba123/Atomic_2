#' GTEx Summary Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param S_G_T_Matrix  A dataframe in the format of SNP GENE TISSUE to be provides
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' S_G_T_Matrix <- data.frame(SNP = c("chr3_49771630" , "chr3_49772408" , "chr3_49773528" , "chr3_49773614" , "chr3_49773845") , Gene = c("ENSG00000182179.12" , "ENSG00000182179.12" ,"ENSG00000182179.12" , "ENSG00000182179.12" ,"ENSG00000182179.12") , Tissue = c("Whole_Blood","Lung","Liver","Pancreas","Ovary"))
#' GTEx_Summary(S_G_T_Matrix)
#' @export
GTEx_Summary <- function(S_G_T_Matrix)
{
  base_url <- "http://172.15.1.20:8000"
  
  S_G_T_Matrix <- toJSON(S_G_T_Matrix)
  res <- GET(url = paste0(base_url, "/K_gTExSumStats"), query = list(S_G_T_Matrix = S_G_T_Matrix))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(fromJSON(aa)))
}

#' Gtex_variance_finder Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param Gene Gene Of Interest in the format ENSG00000182179.12
#' @param Tissue The Tissue of Interest as per GTEx Portal
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' Gtex_variance_finder(Gene = "ENSG00000182179.12 " , Tissue = "Whole_Blood")
#' @export
Gtex_variance_finder <- function(Gene,Tissue)
{
  base_url <- "http://172.15.1.20:8000"
  
  xx <- toJSON(c(Gene,Tissue))
  res <- GET(url = paste0(base_url, "/Gtex_variance_finder"), query = list(gene_tissue = xx))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}

#' Gtex_Sample_Finder Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param Tissue The Tissue of Interest as per GTEx Portal
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' Gtex_Sample_Finder(Tissue="Whole_Blood")
#' @export
Gtex_Sample_Finder <- function(Tissue)
{
  base_url <- "http://172.15.1.20:8000"
  
  res <- GET(url = paste0(base_url, "/Gtex_Sample_Finder"), query = list(Tissue = Tissue))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}





#' Gtex_Sample_Finder Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param chr The Tissue of Interest as per GTEx Portal
#' @param start The Tissue of Interest as per GTEx Portal
#' @param end The Tissue of Interest as per GTEx Portal
#' @return set of snps
#' @examples 
#' snp_gtex(chr = 1, start = 13526, end = 14000)
#' @export
snp_gtex <- function(chr , start, end)
{
  base_url <- "http://172.15.1.20:8000"
  
  xx <- toJSON(c(chr,start,end))
  res <- GET(url = paste0(base_url, "/snp_gtex"), query = list(cumulative = xx))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}




#' snp_1kgp Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param chr The Tissue of Interest as per GTEx Portal
#' @param start The Tissue of Interest as per GTEx Portal
#' @param end The Tissue of Interest as per GTEx Portal
#' @return set of snps
#' @examples 
#' snp_1kgp(chr = 1 , start = 13500 , end = 14000)
#' @export
snp_1kgp <- function(chr , start , end)
{
  base_url <- "http://172.15.1.20:8000"
  
  xx <- toJSON(c(chr,start,end))
  res <- GET(url = paste0(base_url, "/snp_1kgp"), query = list(cumulative = xx))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}




#' snp_gtex_signif Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param gene The Tissue of Interest as per GTEx Portal
#' @param tissue The Tissue of Interest as per GTEx Portal
#' @return set of snps
#' @examples 
#' snp_gtex_signif(gene = "ENSG00000225630" , tissue= "Liver")
#' @export
snp_gtex_signif <- function(gene = "ENSG00000225630" , tissue= "Liver")
{
  base_url <- "http://172.15.1.20:8000"
  
  xx <- toJSON(c(gene,tissue))
  res <- GET(url = paste0(base_url, "/snp_gtex_signif"), query = list(cumulative = xx))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  return(fromJSON(aa))
}
