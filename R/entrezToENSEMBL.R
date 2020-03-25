#' entrezToENSEMBL
#'
#' @param entrez
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[data.frame]]
entrezToENSEMBL <- function(entrez){
  tbl <- toTable(org.Hs.egENSEMBL)
  ens.pathway <- tbl[match(entrez, tbl$gene_id),2]
  if(ens.pathway %>% is.na() %>% sum() == length(entrez)){
    message("Could not map any entrez IDs to ENSEMBL ID")
    return(NULL)
  }
  return(ens.pathway)
}
