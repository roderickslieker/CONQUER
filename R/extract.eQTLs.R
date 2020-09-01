#' Extract eQTLs
#' @directory The directory of the SNP files
#' @SNPs A list of SNPs if `NULL` eQTLs for all SNPs are extracted
#' @return A `dataframe` of nominal significant eQTLs
#' @export
extract.eQTLs <- function(directory, SNPs)
{
  abstractData <- lapply(SNPs,function(x){
    load(paste0(directory,"/",x,".RData"))
    get(x)
  })
  #Extract QTLs from abstractData
  message("Loading eQTLs...")
  CiseQTL <- lapply(1:length(abstractData),function(i){
    if(nrow(abstractData[[i]]$eQTLs) != 0) return(abstractData[[i]]$eQTLs)
  }) %>% do.call(what=rbind)

  TranseQTL <- lapply(1:length(abstractData),function(i){
    if(nrow(abstractData[[i]]$eQTLsTrans) != 0) return(abstractData[[i]]$eQTLsTrans)
  }) %>% do.call(what=rbind)

  eQTLs <- dplyr::bind_rows(CiseQTL,TranseQTL); rm(CiseQTL,TranseQTL)
  eQTLs <- eQTLs[eQTLs$pValue <= 0.05,]
  return(eQTLs)
}
