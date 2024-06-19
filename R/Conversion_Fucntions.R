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
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("rsidToGRCh37(rsid=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
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
  rsids <- paste0(rsids , collapse = "','")
  ex = paste0("rsidToGRCh38(rsid=c('" , rsids , "'))")
  return(RSeval(expr = ex , c = connec))
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
  Chr <- paste0(variant_matrix[,1] , collapse = "','")
  Pos <- paste0(variant_matrix[,2] , collapse = ",")
  ex = paste0("SNP_To_Allele_1KGP(variant_matrix <- data.frame(Chromosome=c('" , Chr ,"') , Position=c(" ,Pos, ")))")
  return(RSeval(expr = ex , c = connec))
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
  Chr <- paste0(variant_matrix[,1] , collapse = "','")
  Pos <- paste0(variant_matrix[,2] , collapse = ",")
  ex = paste0("SNP_To_Allele_Gtex(variant_matrix <- data.frame(Chromosome=c('" , Chr ,"') , Position=c(" ,Pos, ")))")
  return(RSeval(expr = ex , c = connec))
}

##################################################################################


#' Lifthg19Tohg38 Function
#'
#' Lift-over genomic coordinates from hg19 to hg38 position
#' @param rsids A vector of Chromosome and a vector of Position beloning to hg19 reference Frame.
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' Chromosome <- c(8,8,8,8,8,8,8,8,8,8)
#' Position <- c(101592213,106973048,108690829,108681675,103044620,104152280,100479917,103144592,108166508,105201080)
#' Lifthg19Tohg38(Chromosome , Position)
#' @export
Lifthg19Tohg38 <- function(Chromosome , Position)
{
  Chromosome <- paste0(Chromosome , collapse = ",")
  Position <- paste0(Position , collapse = ",")
  ex = paste0("hg19Tohg38(seqnames <- c(" , Chromosome ,") , ranges <- c(" ,Position, "))")
  return(RSeval(expr = ex , c = connec))
}


