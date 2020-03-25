#' getChromatinInteractionData
#'
#' @param chr
#' @param topHits
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[IRanges]]
getChromatinInteractionData <- function(chr, topHits){
  RangeOfInterest <- GenomicRanges::GRanges(seqnames = paste0("chr", chr),
                                            ranges = IRanges::IRanges(min(topHits$start),
                                                                      end=max(topHits$start)))

  from.overlap <- GenomicRanges::findOverlaps(conquer.db::ChromatinInteractions, RangeOfInterest)@from

  to.overlap <- GenomicRanges::findOverlaps(conquer.db::ChromatinInteractions$to, RangeOfInterest)@from

  overlap <- conquer.db::ChromatinInteractions[unique(c(from.overlap, to.overlap)),]

  if (nrow(GenomicRanges::as.data.frame(overlap)) > 0) {
    overlap.bed <- chromInteractionsToBED(overlap, RangeOfInterest)
    return(overlap.bed)
  } else {
    return(NULL)
  }
}


