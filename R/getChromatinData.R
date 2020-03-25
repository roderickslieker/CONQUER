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
getChromatinData <- function(chr, startPos, endPos){
  query <- GenomicRanges::GRanges(seqnames=paste0("chr", chr),
                   ranges=IRanges::IRanges(start = startPos, end = endPos))
  ChromatinData <- IRanges::subsetByOverlaps(conquer.db::Chromatin,query)
  return(ChromatinData)
}

