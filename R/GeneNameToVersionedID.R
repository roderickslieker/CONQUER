#' GeneNameToVersionedID
#'
#' @param geneNames
#'
#' @keywords internal
#' @usage NULL
#' @import rio
#' @return [[data.frame]]
GeneNameToVersionedID <- function(geneNames){
  query <- sprintf("https://gtexportal.org/rest/v1/reference/gene?geneId=%s&gencodeVersion=v26&genomeBuild=GRCh38%%2Fhg38&pageSize=250&format=tsv",
                   paste(geneNames, collapse = "%2C"))
  tryCatch({
    geneCodesDF <- rio::import(query, format = '\t')
    IDs <- geneCodesDF$gencodeId
  },
  error=function(e) {
    errorNum <- gsub(".*?([0-9]+).*", "\\1", e)
    httr::warn_for_status(as.numeric(errorNum),
                          task = "get median expression data")
    return(IDs <- NULL) #if API error, return NULL
  })


  return(IDs)
}
