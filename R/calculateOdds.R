#' calculateOdds
#'
#' @param pathway [[character]] A character vector containing the genes that are included in a pathway.
#' @param geneSet [[character]] A character vector containing the found genes.
#' @param background [[character]] A character vector containing a set of background genes.
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[vector]]
calculateOdds <- function(pathway,geneSet,background){
  tp <- length(intersect(geneSet,pathway))
  fp <- length(pathway) - tp
  fn <- length(geneSet) - tp
  tn <- background - length(pathway) - length(geneSet) + tp
  table <- matrix(c(tp,fp,fn,tn),nrow=2)
  res <- stats::fisher.test(table)
  odds <- res$estimate
  ci_1 <- res$conf.int[1]
  ci_2 <- res$conf.int[2]
  pval <- res$p.value
  return(c(odds,ci_1,ci_2,background,pval,tp))
}
