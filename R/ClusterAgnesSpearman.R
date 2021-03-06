#' Cluster Agnes /Spearman
#' @import cluster
#' @return [[data.frame]]
ClusterAgnesSpearman <- function(expressionData){
  message("Performing hierarchical clustering with flexible UPGMA...")
  message("Calculating Gap statistic to determine optimal K...")
  Data <- log10(expressionData+1) %>% t() %>% scale()
  K.max <- ifelse(nrow(Data) < 30, nrow(Data)-1, 30)
  gap_stat <- cluster::clusGap(Data, FUN = agnesPearson,
                               K.max = K.max, B = 5)

  k <- cluster::maxSE(gap_stat$Tab[,3],
                      gap_stat$Tab[,4],
                      method="globalSEmax")
  message(paste0("Optimal K = ",k))
  ag.res <- agnesPearson(x = Data, k = k)
  return(ag.res)

}
