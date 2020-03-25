#' getExons
#'
#' @param Genes
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[data.frame]]
getExons <- function(Genes) {
  #make columns to fill for data.frame
  parent <- character()
  exend <- integer()
  exstart <- integer()

  #GET exon info for every gene
  for(gene in Genes$gene_id) {
    #get exon info
    EXON <- jsonlite::fromJSON(paste0("http://rest.ensembl.org/lookup/id/",
                                      gene, "?expand=1"),
                               flatten = T)
    #take the transcript that is canonical and add to column variables
    canon <- match(1,EXON$Transcript$is_canonical)
    parent <- as.character(c(parent, rep(
      EXON$Transcript$Parent[1],
      (length(EXON$Transcript$Exon[[canon]]$end)+2))))
    exend <- c(exend, EXON$Transcript$Exon[[canon]]$end)
    exstart <- c(exstart, EXON$Transcript$Exon[[canon]]$start)
    #add UTR regions to exon lists
    if (!is.null(EXON$Transcript$Translation.end[canon])) {
      exend <- c(exend,
                 Genes[Genes$gene_id==gene,]@ranges@start,#5'UTR end
                 EXON$Transcript$Translation.end[canon])#3'UTR end
      exstart <- c(exstart,
                   EXON$Transcript$Translation.start[canon],#5'UTR start
                   GenomicRanges::end(Genes[Genes$gene_id==gene,]))#3'UTR start
    } else {parent <- head(parent, -2)}
  }
  #add column variables to data.frame
  expos <- data.frame(parent, exstart=exstart, exend=exend)
  return(expos)
}
