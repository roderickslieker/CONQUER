#' @importFrom reshape2 colsplit
geteQTLdata <- function(lead, lead.pos, Genes, parallel, allTissues, precalculated) {

  if(length(Genes) == 0) {
    return(data.frame())
  }

  if(!is.null(allTissues)){
  }else{
    allTissues <- conquer.db::gtexTissuesV8
  }

  if(precalculated)
  {
    eQTLs <- get_eQTL_precalculated(lead, tissues = allTissues)
  }else{
    eQTLs <- lapply(X = Genes$gene_id, FUN = get_eQTL_bulk, lead=lead, tissues = allTissues) %>% do.call(what=rbind)
  }





  if(is.null(eQTLs)){
    return(data.frame())
  }else if(nrow(eQTLs) == 0) {
    return(data.frame())
  }else{
    colnames(eQTLs)[colnames(eQTLs) == "snpId"] <- "SNP"
    colnames(eQTLs)[colnames(eQTLs) == "tissueSiteDetailId"] <- "tissue"
    colnames(eQTLs)[colnames(eQTLs) == "geneSymbol"] <- "gene"
    eQTLs$Pval.ratio <- eQTLs$pValue / eQTLs$pValueThreshold
    return(eQTLs)
  }
}






