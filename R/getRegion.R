#' getRegion
#'
#' @param ChromatinInteractionData
#' @param mainSNP
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[vector]]
getRegion <- function(ChromatinInteractionData, mainSNP){
  startPosition <- as.integer(mainSNP$start)
  if (is.null(ChromatinInteractionData)) {
    geneStart <- ifelse(startPosition - 500000 < 0, 0, startPosition - 500000)
    geneEnd <- startPosition + 500000
  } else {
    #use chromosome interaction data to get region for genes (minimum 1Mb region)
    geneStart <- ChromatinInteractionData[[3]]$start
    geneEnd <- ChromatinInteractionData[[3]]$end
    #if start 1Mb window below zero, make it zero
    oneMbLimit <- ifelse(startPosition - 500000 < 0, 0, startPosition-500000)
    geneStart <- ifelse(geneStart < oneMbLimit, geneStart, oneMbLimit)
    #if chromosome interaction region ends before 1Mb window, use 1mb window end
    geneEnd <- ifelse(geneEnd > startPosition + 500000, geneEnd, startPosition + 500000)
  }
  return(c("startPos" = geneStart, "endPos" = geneEnd))
}

