#' Extract enriched pathways
#' @SNPSummary A SNPSummary object
#' @Database Optional string. One of the following: GO_Biological_Process_2018,
#' GO_Cellular_Component_2018, GO_Molecular_Function_2018, GWAS_Catalog_2019,
#' KEGG_2019_Human, Panther_2016, Reactome_2016, OMIM_Disease.
#' @return A `dataframe` of enriched pathways
#' @export
extract.pathways <- function(SNPSummary, subset=NULL)
{
  tissues <- colnames(SNPSummary)
  all.pathways <- lapply(tissues, function(t.sel)
  {
    snp.sub <- SNPSummary[["allOR",t.sel]]

    keep <- lapply(seq_along(snp.sub), function(i){!is.null(snp.sub[[i]])}) %>% do.call(what=c)
    snp.sub <- snp.sub[keep]

    for(i in 1:length(snp.sub)){
      snp.sub[[i]]$Module <- names(snp.sub)[i]
    }

    all.pathway <- do.call(rbind, snp.sub)
    all.pathway <- all.pathway[all.pathway$Adjusted.P.value <= 0.05,]
    all.pathway$Tissue <- t.sel
    all.pathway
  }) %>% do.call(what=rbind)


  if(!is.null(subset))
  {
    all.pathways <- all.pathways[all.pathways$DB %in% subset,]
  }
  all.pathways <- all.pathways[order(all.pathways$Adjusted.P.value, decreasing = F),]

  return(all.pathways)
}
