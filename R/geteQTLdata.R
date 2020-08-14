geteQTLdata <- function(lead, Genes, parallel, allTissues) {

  if(length(Genes) == 0) {
    return(data.frame())
  }

  if(!is.null(allTissues)){
  }else{
    allTissues <- conquer.db::gtexTissuesV8
  }

  #if(parallel){
  #  res <- BiocParallel::bplapply(X = Genes$gene_id,
  #                    FUN = function(x, get_eQTL_bulk, lead, allTissues){
  #                      get_eQTL_bulk(x, lead, allTissues)
  #                    },
  #                    get_eQTL_bulk=get_eQTL_bulk,
  #                    lead=lead,
  #                    allTissues=allTissues)

  #}else{
  res <- lapply(X = Genes$gene_id, FUN = get_eQTL_bulk, lead=lead, tissues = allTissues)
  #}
  eQTLs <- do.call(rbind, res)
  colnames(eQTLs)[colnames(eQTLs) == "snpId"] <- "SNP"
  colnames(eQTLs)[colnames(eQTLs) == "tissueSiteDetailId"] <- "tissue"
  colnames(eQTLs)[colnames(eQTLs) == "geneSymbol"] <- "gene"

  if (nrow(eQTLs) == 0) {
    return(data.frame())
  }
  #add pVal.ratio for significant cis genes
  eQTLs$Pval.ratio <- eQTLs$pValue / eQTLs$pValueThreshold
  return(eQTLs)
}






