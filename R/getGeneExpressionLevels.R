#' geteQTLdata
#'
#' @param allIDs
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[data.frame]]
getGeneExpressionLevels <- function(allIDs) {
  queryGenes <- paste(allIDs, collapse = "%2C")
  #get median expression
  query <- paste0("https://gtexportal.org/rest/v1/expression/",
                  "medianGeneExpression?datasetId=gtex_v8&gencodeId=",
                  queryGenes,
                  "&hcluster=false&format=json")
  tryCatch({
    geneExpr <- as.data.frame(jsonlite::fromJSON(query, flatten = T)$medianGeneExpression)
  },
  error=function(e) {
    errorNum <- gsub(".*?([0-9]+).*", "\\1", e)
    httr::warn_for_status(as.numeric(errorNum),
                          task = "get median expression data")
    return(geneExpr <- data.frame()) #if API error, return NULL
  })
  return(geneExpr)
}
