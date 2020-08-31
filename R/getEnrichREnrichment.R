#' importFrom httr POST
#' import enrichR
#' import ggplot2
getEnrichREnrichment <- function(tissue, SNPSummary)
{

  genSets <- SNPSummary[["Module_Genes",tissue]]

  if(class(genSets) == "character")
  {
    genSets <- list(`1` = genSets)
  }
  ModuleOdds <- lapply(genSets, function(genes.ens){
    if(length(genes.ens) > 2)
    {
      #Convert genes
      r <- httr::POST("https://biotools.fr/human/ensembl_symbol_converter/",
                      body = list(api=1, ids=toJSON(genes.ens)))
      genes.symbol <- rawToChar(r$content) %>% jsonlite::fromJSON() %>% as.character()

      pws <- enrichR::enrichr(genes = genes.symbol, databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018","GO_Molecular_Function_2018","KEGG_2019_Human","Panther_2016","GWAS_Catalog_2019","Reactome_2016","OMIM_Disease"))

      pws <- lapply(1:length(pws), function(x){
        #cat(x)
        temp <- pws[[x]]
        if(nrow(temp) >= 1)
          temp$DB <- names(pws)[x]
        temp
      }) %>% do.call(what=rbind)
      pws$DB <- factor(pws$DB)
      return(pws)
    }

  })
  names(ModuleOdds) <- 1:length(ModuleOdds)
  return(ModuleOdds)
}


