#' geteQTLdata
#'
#' @param chr
#' @param startPos
#' @param endPos
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[GRanges]]
getGenesInRegion <- function(chr, startPos, endPos) {
  #format region
  region <- paste(chr,
                  paste(startPos,
                        endPos,
                        sep = '..'),
                  sep = ':')
  #GET genes in region
  tryCatch({
    GENES <- jsonlite::fromJSON(
      paste0("http://rest.ensembl.org/overlap/region/human/",
             region, "?feature=gene"),
      flatten = T)},
    error=function(e){
      errorNum <- gsub(".*?([0-9]+).*", "\\1", e)
      stop_for_status(as.numeric(errorNum),
                      paste("get genes in", region))
    })
  if (is.null(nrow(GENES))) {return(GenomicRanges::GRanges())}
  GENES <- GENES[is.element(GENES$biotype,
                            c("protein_coding", "snRNA", "lincRNA", "snoRNA",
                              "ncRNA", "sRNA", "macro_lncRNA", "tRNA",
                              "scaRNA", "miRNA", "lncRNA", "rRNA")),]
  if (nrow(GENES)==0) {return(GenomicRanges::GRanges())}
  #to GRanges object
  gr.Genes <- GenomicRanges::GRanges(
    seqnames = paste0("chr",GENES$seq_region_name),
    ranges = IRanges::IRanges(start = GENES$start, end = GENES$end),
    strand = GENES$strand,
    name = GENES$external_name,
    gene_id = GENES$gene_id,
    biotype = GENES$biotype
  )
  return(gr.Genes)
}
