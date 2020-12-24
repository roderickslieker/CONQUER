#' @importFrom BiocGenerics start end
getDataforSingleSNP <- function(variant, precalculated, directory=NULL, token=NULL, population="CEU",Chromatin, allTissues){
  message(sprintf("Retrieving data for: %s",variant))

  tryCatch({
    #Retrieve Ensembl information
    message("getVariantInfoEnsembl")
    mainSNP <- getVariantInfoEnsembl(variant,
                                     population = population)

    #Retrieve Linkage Disequilibrium Block from NIH LDlink
    message("getLDLink")
    LinkDis <- getLDLink(SNP = mainSNP,
                         token = token,
                         population = population)

    #Select top hits from LinkDis table
    topHits <- LinkDis[LinkDis$r2>=0.8,]

    #Retrieve Chromatin interaction Data
    message("getChromatinInteractionData")
    ChromatinInteractionData <- getChromatinInteractionData(chr = mainSNP$chr,
                                                            topHits = topHits)

    #Get Start and End position of region of interest
    region <- getRegion(ChromatinInteractionData = ChromatinInteractionData,
                        mainSNP =  mainSNP)

    #Get genes within region of interest
    Genes <- getGenesInRegion(chr = mainSNP$chr,
                              startPos = region["startPos"],
                              endPos = region["endPos"])

    #Get ChromatinData
    ChromatinData <- getChromatinData(chr = mainSNP$chr,
                                      startPos = min(BiocGenerics::start(Genes)),
                                      endPos = max(BiocGenerics::end(Genes)),
                                      Chromatin=Chromatin)

    #Get exons for locusZoomPlot
    query <- GenomicRanges::GRanges(seqnames=paste0("chr", mainSNP$chr),
                                    ranges=IRanges::IRanges(start = min(LinkDis$start), end = max(LinkDis$start)))
    locusGenes <- IRanges::subsetByOverlaps(Genes,query)
    locusExons <- getExons(locusGenes)

    #Transcription factors
    TFquery <- GenomicRanges::GRanges(seqnames=paste0("chr", mainSNP$chr),
                                      ranges=IRanges::IRanges(start = min(topHits$start) - 10000,
                                                              end = max(topHits$start) + 10000))
    TFdata <- IRanges::subsetByOverlaps(conquer.db::TranscriptionFactorsAll,TFquery)

    #cis eQTLs
    cisQuery <-  GenomicRanges::GRanges(seqnames=paste0("chr", mainSNP$chr),
                                        ranges=IRanges::IRanges(start = ifelse(as.integer(mainSNP$start) - 500000 < 0, 0,
                                                                               as.integer(mainSNP$start) - 500000),
                                                                end = as.integer(mainSNP$start) + 500000))

    cisGenes <- IRanges::subsetByOverlaps(Genes,cisQuery)
    message("geteQTLdata")
    lead.pos <- paste0(mainSNP$chr, "_", mainSNP$start, "_")
    eQTLdata <- geteQTLdata(lead = mainSNP$variation,
                            lead.pos = lead.pos,
                            Genes = cisGenes,
                            allTissues=allTissues,
                            precalculated=precalculated)

    #trans eQTLs
    transGenes <- Genes[!Genes$name %in% cisGenes$name,]
    if(length(transGenes) != 0) {
      lead.pos <- paste0(mainSNP$chr,"_",mainSNP$start,"_")
      transeQTLdata <- geteQTLdata(lead = mainSNP$variation,
                                   lead.pos = lead.pos,
                                   Genes = transGenes,
                                   parallel = parallel,
                                   allTissues=allTissues,
                                   precalculated = precalculated)
    }else {
      transeQTLdata <- data.frame()
    }

    #Gene expression
    message("Gene expression")
    VersionedIDs <- GeneNameToVersionedID(Genes$gene_id)
    blocks <- split(VersionedIDs, ceiling(seq_along(VersionedIDs) / 90))
    blockExpressions <- lapply(X = blocks, FUN = getGeneExpressionLevels)
    geneExpression <-  do.call(rbind, blockExpressions)

    #Combine all data
    assign(variant,list("SNP" = mainSNP,
                        "topHits" = topHits,
                        "LD" = LinkDis,
                        "Exons" = locusExons,
                        "chromInt" = ChromatinInteractionData,
                        "chromatin" = ChromatinData,
                        "TFs" = TFdata,
                        "eQTLs" = eQTLdata,
                        "eQTLsTrans" = transeQTLdata,
                        "geneExpr" = geneExpression,
                        "genes" = Genes
                        )
    )

    save(list = variant, file = sprintf("%s/%s.RData",directory,variant))
    message(sprintf("%s successful!",variant))
    return(c("SNP" = variant,"status" = TRUE))
  },
  error=function(e) {
    message(sprintf("%s failed!",variant))
    return(c("SNP" = variant,"status" = FALSE))
  })

}
