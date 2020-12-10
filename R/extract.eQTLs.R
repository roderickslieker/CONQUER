#' Extract eQTLs
#' @directory The directory of the SNP files
#' @SNPs A list of SNPs if `NULL` eQTLs for all SNPs are extracted
#' @return A `dataframe` of nominal significant eQTLs
#' @export
extract.eQTLs <- function(directory, SNPs)
{
  extractedData <- lapply(SNPs,function(x){
    load(paste0(directory,"/",x,".RData"))
    get(x)
  })
  #Extract QTLs from extractedData
  message("Loading eQTLs...")
  CiseQTL <- lapply(1:length(extractedData),function(i){
    if(nrow(extractedData[[i]]$eQTLs) != 0) return(extractedData[[i]]$eQTLs)
  }) %>% do.call(what=rbind)
  CiseQTL$Type <- "Cis"

  TranseQTL <- lapply(1:length(extractedData),function(i){
    if(nrow(extractedData[[i]]$eQTLsTrans) != 0) return(extractedData[[i]]$eQTLsTrans)
  }) %>% do.call(what=rbind)
  TranseQTL$Type <- "Trans"

  eQTLs <- dplyr::bind_rows(CiseQTL,TranseQTL); rm(CiseQTL,TranseQTL)
  eQTLs <- eQTLs[eQTLs$pValue <= 0.05,]
  return(eQTLs)
}


extract.eQTLs.noLoad <- function(extractedData)
{
  #Extract QTLs from extractedData
  message("Loading eQTLs...")
  CiseQTL <- lapply(1:length(extractedData),function(i){
    if(nrow(extractedData[[i]]$eQTLs) != 0) return(extractedData[[i]]$eQTLs)
  }) %>% do.call(what=rbind)
  CiseQTL$Type <- "Cis"

  TranseQTL <- lapply(1:length(extractedData),function(i){
    if(nrow(extractedData[[i]]$eQTLsTrans) != 0) return(extractedData[[i]]$eQTLsTrans)
  }) %>% do.call(what=rbind)
  TranseQTL$Type <- "Trans"

  eQTLs <- dplyr::bind_rows(CiseQTL,TranseQTL); rm(CiseQTL,TranseQTL)
  eQTLs <- eQTLs[eQTLs$pValue <= 0.05,]
  return(eQTLs)
}
