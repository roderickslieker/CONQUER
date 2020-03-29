#' agnesPearson
#'
#' @param x [[matrix]] A matrix containg the gene expression.
#' @param k [[numeric]] A integer number indicating the number of clusters.
#'
#' @keywords internal
#' @usage NULL
#' @import cluster
#' @importFrom stats cutree
#' @return [[list]]
agnesPearson <- function(x,k){
  distance <- 1 -  cor(x = t(x) ,method = "spearman")
  agnesObj <- cluster::agnes(x = distance, diss = T, method="gaverage")
  assignments <- stats::cutree(agnesObj,k=k)
  names(assignments) <- rownames(distance)
  return(list("cluster" = assignments))
}
