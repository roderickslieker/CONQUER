#' getChromatinData
#'
#' @param chr
#' @param startPos
#' @param endPos
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[IRanges]]
getChromatinData <- function(chr, startPos, endPos,Chromatin){
  query <- GenomicRanges::GRanges(seqnames=paste0("chr", chr),
                   ranges=IRanges::IRanges(start = startPos, end = endPos))
  ChromatinData <- IRanges::subsetByOverlaps(Chromatin,query)
  return(ChromatinData)
}

