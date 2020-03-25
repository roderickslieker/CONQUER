#' getCoexpression
#'
#' @param geneSet Vector of genes
#' @param expressionData data.frame with expression data with samples as rows and genes as columns
#' @param method Which correlation coefficient to calculate "pearson", "spearman" or "kendall"
#' @param threshold Threshold of minimum correlation coefficient on which to base the co-expression
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[vector]]
getCoexpression <- function(geneSet,expressionData, method = "spearman", threshold = 0.90){
  #Subset geneSet from expression data
  originalSet <- expressionData[,geneSet]
  #Perform correlation analysis
  correlationMatrix <- cor(originalSet, expressionData, method = method)
  #Extract columns in which there is a coefficient higher then the threshold
  coExpressed <- correlationMatrix[,apply(correlationMatrix,2,max) >= threshold] %>% colnames()
  #Combine coexpressed genes with original geneSet
  coExpressed <- c(coExpressed,geneSet) %>% unique()
  return(coExpressed)
}
