#' chromInteractionsToBED
#'
#' @param interactions [[GRanges]] A Genomic ranges object containing the chromatin interactions.
#' @param RangeOfInterest [[]]
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[list]]
chromInteractionsToBED <- function(interactions, RangeOfInterest) {
  #split interactions GRanges into 2 bed format dataframes
  #GRanges to dataframe
  df <- GenomicRanges::as.data.frame(interactions)
  #make dataframe of region to plot
  region <- data.frame(chr=RangeOfInterest@seqnames[1],
                       start=min(c(df$start, df$to.start,
                                   GenomicRanges::start(RangeOfInterest))),
                       end=max(c(df$end, df$to.end,
                                 GenomicRanges::end(RangeOfInterest))))
  #make bed format dataframe of from ranges
  from <- data.frame(chr=df$seqnames, start=df$start, end=df$end,
                     gene=df$from.gene, tissue=df$cell.tissue, id=1:nrow(df))
  #make bed format dataframe of to ranges
  to <- data.frame(chr=df$to.seqnames, start=df$to.start, end=df$to.end,
                   gene=df$to.gene, id=1:nrow(df))
  #format chr columns
  from$chr <- as.character(from$chr)
  to$chr <- as.character(to$chr)
  #return dataframes in a list
  return(list(from, to, region))
}
