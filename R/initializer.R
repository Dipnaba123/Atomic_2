#' Initialization Function
#'
#' This Function should be called only once before using any function of the package. 
#' @param pkg_name The name of the packages
#' @return NA
#' @examples 
#' g_initialize("Gtex_Summary");
#' g_initialize( c("Gtex_Summary", "1000GP_LD_Summary", "Gtex_LD_Summary" , "UKBB_BP_Enrichment" , "Haploreg" , "Genomic_Regions" ));
#' @export
g_initialize <- function(func_name){
  
  if (!require("RSclient", quietly = TRUE))
    install.packages("RSclient")
  library("RSclient")
  
  connec <<- RSconnect(host = "172.15.1.20" , port = 8080)
  
  if(length(grep(func_name,pattern = "Gtex_Summary")) != 0 )
  {
    print("----- Dependencies Of GTEx Summary Functions uploaded ------")
    RSeval(expr = "source('/home/diptanil/Rcode/gtex/k_gtex_fast_v4.R')" , connec)
  }
  
  if(length(grep(func_name,pattern = "1000GP_LD_Summary")) != 0 )
  {
    print("----- Dependencies Of 1000GP_LD_Summary Functions uploaded ------")
    RSeval(expr = "source('/home/diptanil/Rcode/LD/codes/calculate_ld.R')" , connec)
  }
  
  if(length(grep(func_name,pattern = "Gtex_LD_Summary")) != 0 )
  {
    print("----- Comming Soon ------")
  }
  
  if(length(grep(func_name,pattern = "UKBB_BP_Enrichment")) != 0 )
  {
    print("----- Comming Soon ------")
  }
  
  if(length(grep(func_name,pattern = "Haploreg")) != 0 )
  {
    print("----- Dependencies Of Haploreg Functions uploaded ------")
    RSeval(expr = "source('~/Rcode/haploreg_parser/final_codes/haploreg_initializer.R')" , connec)
  }
  
  if(length(grep(func_name,pattern = "Genomic_Regions")) != 0 )
  {
    print("----- Dependencies Of Genomic_Regions Functions uploaded ------")
    RSeval(expr = "source('~/Rcode/hash_functions/codes/rsidToGRCh37.R')" , connec)
    RSeval(expr = "source('~/Rcode/hash_functions/codes/rsidToGRCh38.R')" , connec)
    RSeval(expr = "source('~/Rcode/Genomic_Regions/SNP_To_Allele.R')" , connec)
    RSeval(expr = "source('/home/diptanil/Rcode/Genomic_Regions/Lifthg19Tohg38.R')" , connec)
  }
  return(0)
}

#' Closing Function
#'
#' This Function should be called only once at the end of the program
#' @param con The connection to the server.
#' @return NA
#' @examples 
#' g_free(connec)
#' @export
g_free <- function(con)
{
  RSclose(con)
}