#' snpToAF Function
#'
#' Retrieve allele frequencies for a set of SNPs from the HAPLOREG database.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples
#' snpToAF(rsids = "rs10")
#' snpToAF(rsids = c("rs10","rs3","rs4","rs10047249"))
#' @export
snpToAF <- function(rsids)
{
  
  # Base URL of the API
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToAF"), query = list(rsids = paste0(rsids , collapse = "_")))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  
  aan <- str_extract_all(aa, "\\d+\\.\\d+") %>%  # Extract all sequences of digits with decimal points
    unlist() %>%                                   # Flatten the list to a vector
    as.numeric() 
  
  aan[aan == 1.57] <- NA
  aand <- matrix(aan, ncol = 4, byrow = TRUE) %>% 
    as.data.frame()%>%
    setNames(c("AFR","AMR","ASN","EUR"))
  return(aand)
}


#' snpToCHROMHMM_15STATE Function
#'
#' Retrieve the CHROMHMM 15-state annotations for a set of SNPs from the HAPLOREG database.
#' @param rsid A single variant with rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToCHROMHMM_15STATE(rsid <- "rs10")
#' @export
snpToCHROMHMM_15STATE <- function(rsid)
{
  # Base URL of the API
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToCHROMHMM_15STATE"), query = list(rsid = rsid))
  aa <- content(res, as ="text" , encoding = "UTF-8")

  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    filter(raw_data != "") %>%  # Remove empty strings
    separate(raw_data, into = c("CELL","State_Number"), sep = "\\|")  # Split the string into two columns
  
  
  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    filter(raw_data != "") %>%  # Remove empty strings
    separate(raw_data, into = c("CELL", "State_Number"), sep = "\\|") %>%  # Split the string into two columns
    separate(State_Number, into = c("State_Number", "State_Description"), sep = "_")  # Split the Category column into two columns
  
  return(df)
}


#' snpToCHROMHMM_25STATE Function
#'
#' Retrieve the CHROMHMM 25-state annotations for a set of SNPs from the HAPLOREG database.
#' @param rsid A single variant with rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToCHROMHMM_25STATE(rsid <- "rs10");
#' @export
snpToCHROMHMM_25STATE <- function(rsid)
{
  # Base URL of the API
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToCHROMHMM_25STATE"), query = list(rsid = rsid))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  
  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    filter(raw_data != "") %>%  # Remove empty strings
    separate(raw_data, into = c("CELL","State_Number"), sep = "\\|")  # Split the string into two columns
  
  
  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    filter(raw_data != "") %>%  # Remove empty strings
    separate(raw_data, into = c("CELL", "State_Number"), sep = "\\|") %>%  # Split the string into two columns
    separate(State_Number, into = c("State_Number", "State_Description"), sep = "_")  # Split the Category column into two columns
  
  return(df)
}



#' snpToGENCODE Function
#'
#' Retrieve GENCODE information for a set of SNPs from the HAPLOREG database.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToGENCODE(rsids <- "rs10");
#' snpToGENCODE(rsids <- c("rs10","rs3","rs4","rs5"));
#' @export
snpToGENCODE <- function(rsids)
{
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToGENCODE"), query = list(rsids = paste0(rsids , collapse = "_")))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  df <- matrix(fromJSON(aa), ncol = 4, byrow = TRUE) %>%
    as.data.frame() %>%
    setNames(c("Distance","Direction","ID","Name"))
  return(df)
}



#' snpToHISTONEPEAKS Function
#'
#' Retrieve HISTONEPEAKS information for a set of SNPs from the HAPLOREG database.
#' @param rsid genotypic variant as rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToHISTONEPEAKS(rsid <- "rs3");
#' @export
snpToHISTONEPEAKS <- function(rsid)
{
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToHISTONEPEAKS"), query = list(rsid = rsid ))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    filter(raw_data != "") %>%  # Remove empty strings
    separate(raw_data, into = c("CELL", "MARK_NAME"), sep = "\\|") %>%  # Split the string into two columns
    separate(MARK_NAME, into = c("MARK_NAME", "MARK_DESCRIPTION"), sep = "_")  # Split the Category column into two columns
  
  return(df)
}




#' snpToMotif Function
#'
#' Retrieve MOTIF information for a set of SNPs from the HAPLOREG database.
#' @param rsid genotypic variant as rsid
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToMotif(rsid = "rs10")
#' @export
snpToMotif <- function(rsid)
{
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToMotif"), query = list(rsid = rsid ))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  # Convert the vector to a data frame
  df <- data.frame(raw_data = fromJSON(aa)) %>%
    separate(raw_data, into = c("PWM","dd","DELTA","POSITION","STRAND","REFSCORE","ALTSCORE"), sep = "\\|") %>%
    select(-dd)  # Remove the middle column
  
  return(df)
}


#' snpToREFSEQ Function
#'
#' Retrieve GENCODE information for a set of SNPs from the HAPLOREG database.
#' @param rsids A vector of rsids
#' @return Summary Statistics of the SNP,Gene and Tissue
#' @examples 
#' snpToREFSEQ(rsids = "rs10")
#' snpToREFSEQ(rsids = c("rs10","rs3","rs4","rs1"))
#' @export
snpToREFSEQ <- function(rsids)
{
  base_url <- "http://172.15.1.20:8000"
  res <- GET(url = paste0(base_url, "/snpToREFSEQ"), query = list(rsids = paste0(rsids , collapse = "_")))
  aa <- content(res, as ="text" , encoding = "UTF-8")
  df <- matrix(fromJSON(aa), ncol = 4, byrow = TRUE) %>%
    as.data.frame() %>%
    setNames(c("Distance","Direction","ID","Name"))
  return(df)
}

