#' Extract LD
#' @directory The directory of the SNP files
#' @SNPs A list of SNPs if `NULL` eQTLs for all SNPs are extracted
#' @return A `dataframe` with SNPs in LD
#' @export
extract.LD <- function(directory, SNPs)
{
  abstractData <- lapply(SNPs,function(x){
    load(paste0(directory,"/",x,".RData"))
    get(x)
  })
  #Extract QTLs from abstractData
  message("Loading LD...")
  LD <- lapply(1:length(abstractData),function(i){
    if(nrow(abstractData[[i]]$LD) != 0){
    	out <- abstractData[[i]]$LD[c("variation","chr","start","r2","consequence_type")]
    	return(out)
    }
  }) %>% do.call(what=rbind)
  LD <- LD[LD$r2 >= 0.8,]
  return(LD)
}
