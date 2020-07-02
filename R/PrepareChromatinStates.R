#' PrepareChromatinStates
#'
#' @param SNPData
#' @param sampleSel
#'
#' @keywords internal
#' @importFrom reshape2 melt
#'
#' @return [[list]]
PrepareChromatinStates <- function(SNPData, sampleSel = NULL){
  chr <- SNPData$SNP$chr
  LD80 <- SNPData$topHits
  States <- SNPData$chromatin
  if(nrow(LD80) > 5) {
    #make GRanges object of LD range
    SNPsrange <- GenomicRanges::GRanges(
      seqnames = paste0("chr", chr),
      ranges = IRanges::IRanges(start = min(LD80$start),
                                end = max(LD80$start)))
  }else {
    SNPsrange <- GenomicRanges::GRanges(
      seqnames = paste0("chr", chr),
      ranges = IRanges::IRanges(start = as.numeric(SNPData$SNP$start) - 10000,
                                end = as.numeric(SNPData$SNP$start) + 10000))
  }
  ###Find overlap between SNP range and chromatin data###
  StatesinRange <- GenomicRanges::findOverlaps(States, SNPsrange)
  StatesinRange <- States[StatesinRange@from,]
  StatesPlotData <- GenomicRanges::as.data.frame(StatesinRange)


  ###Define Y value for sample###
  sampley <- data.frame(sample = sort(as.character(unique(StatesPlotData$sample)),decreasing = T))
  sampley$Y <- 1:nrow(sampley)
  StatesPlotData$Y <- sampley[match(StatesPlotData$sample,sampley$sample),"Y"]
  StatesPlotData$sample <- factor(as.character(StatesPlotData$sample), levels=sampley$sample)


  if(!is.null(sampleSel)){
    StatesPlotData <- StatesPlotData[StatesPlotData$sample %in% sampleSel,]
  }
  ### SNP positions in plot###
  jit <- lapply(unique(StatesPlotData$Y),function(num){LD80$start})
  jit <- as.data.frame(do.call(rbind,jit))
  jit$index <- unique(StatesPlotData$sample)
  jit <- reshape2::melt(jit,id.vars="index")
  jit$group <- StatesPlotData[match(jit$index,StatesPlotData$sample),"group"]
  ###Color scale###
  fillScale <- as.character(StatesPlotData[match(unique(StatesPlotData$target),StatesPlotData$target),"colr"])
  names(fillScale) <- unique(StatesPlotData$target)
  output <- list(StatesPlotData, fillScale, jit)
  return(output)
}
