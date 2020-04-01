#' ClusterPam
#'
#' @param expressionData
#'
#' @keywords internal
#' @usage NULL
#' @import cluster
#' @return [[data.frame]]
ClusterPam <- function(expressionData){
  set.seed(123)
  message("Calculating Gap statistic to determine optimal K...")
  Data <-  1 - stats::cor(expressionData,method="spearman")
  gap_stat <- cluster::clusGap(Data, FUN = cluster::pam,
                               K.max = 40, B = 20, diss=T)

  k <- cluster::maxSE(gap_stat$Tab[, "gap"],
                      gap_stat$Tab[, "SE.sim"],
                      method="globalSEmax")
  message(paste0("Optimal K = ",k))
  pam.res <- stats::kmeans(Data, centers = k)
  return(pam.res)
}
