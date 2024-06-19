#' LD Generation Function
#'
#' A function in R that retrieves a summary of SNP-gene-tissue associations from the Genotype-Tissue Expression GTEx dataset.
#' @param snp_set_1 set of variant 1
#' @param snp_set_2 set of variant 2
#' @param choice "rsids" or "variant_matrix"
#' @param population Populations like "EUR","SAS","AFR","AMR","EAS","ALL"
#' @return Linkage
#' @examples
#' snp_set_1 = c("rs112155239","rs62636508","rs58108140","rs574746232","rs571093408")
#' snp_set_2 = c("rs112155239","rs62636508","rs58108140","rs574746232","rs571093408")
#' population = "EUR"
#' choice = "rsids"
#' calculate_LD(snp_set_1 = snp_set_1 , snp_set_2 = snp_set_2, choice="rsids" ,population = "EUR")
#' 
#' 
#' snp_set_1 = c("1_13011","1_10440","1_10583","1_13011")
#' snp_set_2 = c("1_13011","1_10440","1_10583","1_13011")
#' population = "EUR"
#' choice = "variant_matrix"
#' calculate_LD(snp_set_1 = snp_set_1 , snp_set_2 = snp_set_2, choice="variant_matrix" ,population = "EUR")
#' @export
calculate_LD <- function(snp_set_1 , snp_set_2 , choice , population)
{
  snp_set_1 <- paste0(snp_set_1 , collapse = "','")
  snp_set_2 <- paste0(snp_set_2 , collapse = "','")
  ex = paste0("calculate_ld(snp_set_1 = c('" , snp_set_1 ,"') , snp_set_2 = c('", snp_set_2, "'),choice='" , choice , "' ,population = '",population,"')")
  return(RSeval(expr = ex , c = connec))
}

