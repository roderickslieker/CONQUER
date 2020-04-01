#' @importFrom BiocGenerics toTable
entrezToENSEMBL <- function(entrez){
  tbl <- BiocGenerics::toTable(org.Hs.egENSEMBL)
  ens.pathway <- tbl[match(entrez, tbl$gene_id),2]
  if(ens.pathway %>% is.na() %>% sum() == length(entrez)){
    message("Could not map any entrez IDs to ENSEMBL ID")
    return(NULL)
  }
  return(ens.pathway)
}
