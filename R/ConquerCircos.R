#' @import viridis
ConquerCircos <-  function(SNPData, tissue=NULL){
  chromatinInteractions <- conquer.db::ChromatinGroups$Interactions
  chromatinStates <- conquer.db::ChromatinGroups$States
  if(is.null(tissue)){
    tissue <- chromatinInteractions$name
    subtissue <- NULL
  }else{
    subtissue <- chromatinStates[chromatinStates$group == tissue, "name"]
    tissue <- chromatinInteractions[chromatinInteractions$group == tissue, "name"]
  }
  mainSNP <- SNPData$SNP$variation
  #Generate link tracks

  chr.start <- SNPData$chromInt[[1]]
  chr.end <- SNPData$chromInt[[2]]
  chr.range <- SNPData$chromInt[[3]]



  chr.range$start <- min(start(SNPData$genes))
  chr.range$end <- max(end(SNPData$genes))
  if(!is.null(tissue)){
    chr.start <- chr.start[chr.start$tissue %in% tissue, ]
    chr.end <- chr.end[match(chr.start$id, chr.end$id), ]
  }
  #Generate LD
  LD <- SNPData$topHits
  LD.min <- min(LD$start)
  LD.max <- max(LD$start)
  LD$color <- ifelse(LD$variation == mainSNP,"red","black")
  # Background track
  link.list.sum <- BioCircos::BioCircosArcTrack(trackname = "background",
                                                chromosomes = as.character(chr.range$chr),
                                                starts = LD.min,
                                                ends = LD.max,
                                                colors = "#A4A4A4",opacities = 0.5,
                                                maxRadius = 1.14,
                                                minRadius = 1.00)
  link.list.sum <- link.list.sum + BioCircos::BioCircosSNPTrack(trackname = "snps",
                                                                chromosomes = as.character(chr.range$chr),
                                                                positions = LD$start,
                                                                values = LD$r2,
                                                                colors = LD$color,
                                                                fill = "grey",
                                                                size = 3,
                                                                labels = sprintf("SNP:%s<br>r2:%s<br>Consequence:%s",
                                                                                 LD$variation, LD$r2, LD$consequence_type),
                                                                maxRadius = 1.13,
                                                                minRadius = 1.01)

  #Generate Gene
  genes <- SNPData$genes
  link.list.sum <- link.list.sum + BioCircos::BioCircosArcTrack(trackname = "genes",
                                                                chromosomes = as.character(chr.range$chr),
                                                                starts = start(genes),
                                                                ends = end(genes),
                                                                labels = sprintf("Name:%s<br>ID:%s",
                                                                                 genes$name, genes$gene_id),
                                                                colors = "#A4A4A4",opacities = 0.5,
                                                                maxRadius = 0.99,
                                                                minRadius = 0.94)


  #Generate histone mods
  if(!is.null(subtissue)){
    states <- SNPData$chromatin
    range.grr <- GenomicRanges::GRanges(chr.range$chr, IRanges(chr.range$start, chr.range$end))
    ol.grr <- GenomicRanges::findOverlaps(query = range.grr, subject = states) %>% as.matrix()
    states <- states[ol.grr[,2],]
    states <- states[states$sample %in% subtissue,]
    targets <- as.character(states$target) %>% unique()
    celllines <- as.character(states$sample) %>% unique()
    nocel <- length(celllines)
    state.rings <- lapply(1:nocel, function(i,nocel)
    {
      st <- 0.93 - (i-1)*0.05 + (-0.01*(i-1))
      en <- 0.88 - (i-1)*0.05 + (-0.01*(i-1))
      out <- data.frame(cellline = celllines[i],st, en)
      return(out)
    }, nocel=nocel) %>% do.call(what = rbind)
    res.chromstates <- lapply(1:nrow(state.rings), function(i, state.rings, states){

      state.rings.sel <- state.rings[i,]
      states.sub <- states[states$sample %in% state.rings.sel$cellline,]

      states.unique <- unique(states.sub$target) %>% as.character()
      states.sub <- states.sub[start(states.sub) >= chr.range$start,]
      states.sub <- states.sub[end(states.sub) <= chr.range$end,]

      link.list.sum <- BioCircos::BioCircosArcTrack(trackname = "ChromatinStates", maxRadius = state.rings.sel$st, minRadius = state.rings.sel$en, label=as.character(states.sub$target),
                                         chromosomes = as.character(chr.range$chr), starts = start(states.sub), ends = end(states.sub), colors = as.character(states.sub$colr))
      return(link.list.sum)
    }, state.rings = state.rings, states=states)

    if(length(res.chromstates) == 1){
      link.list.sum <- link.list.sum + res.chromstates[[1]]
    }else if(length(res.chromstates) == 2){
      link.list.sum <- link.list.sum + res.chromstates[[1]]+res.chromstates[[2]]
    }else{
      for(i in 1:length(res.chromstates)){
        link.list.sum <- link.list.sum + res.chromstates[[i]]
      }
    }
  }

  if(length(chr.start$tissue) == 0){

  }else{
    cols.links <- viridis::viridis_pal()(length(unique(chr.start$tissue)))
    names(cols.links) <- unique(chr.start$tissue)
    chr.start$cols <- cols.links[match(chr.start$tissue, names(cols.links))]
    link.list <-lapply(unique(chr.start$tissue), function(tissue, chr.start, chr.end){
      chr.start.sub <- chr.start[chr.start$tissue %in% tissue,]
      chr.end.sub <- chr.end[chr.start$tissue %in% tissue,]
      tissues <- as.character(chr.start.sub$tissue) %>% as.character()
      labels <- sprintf("<br>Tissue:%s<br>From:%s<br>To:%s",tissues,chr.start.sub$gene, chr.end.sub$gene)
      cols <- substr(x = chr.start.sub$cols, start = 1,7) %>% as.character() %>% unique()
      if(exists("state.rings")){
        maxRadius.links <- min(state.rings$en) - 0.03
      }else{
        maxRadius.links <- 0.8
      }
      links <- BioCircos::BioCircosLinkTrack(trackname = "Links",
                                             gene1Chromosomes = chr.start.sub$chr,
                                             gene1Starts = chr.start.sub$start,
                                             gene1Ends = chr.start.sub$end,
                                             gene2Chromosomes = chr.end.sub$chr,
                                             gene2Starts = chr.end.sub$start,
                                             gene2Ends = chr.end.sub$end,
                                             axisPadding = 6,color = cols, width = "0.1em",
                                             labels = labels, displayLabel = F,
                                             maxRadius = maxRadius.links)
      links

    }, chr.start=chr.start, chr.end=chr.end)
    if(exists(x = "link.list.sum")){
      for(i in 1:length(link.list)){
        link.list.sum <- link.list.sum + link.list[[i]]
      }
    }else{
      if(length(link.list) == 1){
        link.list.sum <- link.list[[1]]
      }
      else if(length(link.list) == 2){
        link.list.sum <- link.list[[1]]+link.list[[2]]
      }else{
        link.list.sum <- link.list[[1]]+link.list[[2]]
        for(i in 3:length(link.list)){
          link.list.sum <- link.list.sum + link.list[[i]]
        }
      }


  }

}
  #Length
  if(!is.null(subtissue))
  {
    min.value <- min(c(start(states), start(genes)))
    max.value <- max(c(end(states), end(genes)))
  }else {
    min.value <- min(start(genes))
    max.value <- max(end(genes))
  }
  len <- max.value - min.value
  genome <- list(len)
  names(genome) <- chr.range[1,1]
  #Plot tracks
  BioCircos::BioCircos(width = 1000,
                       height=1000,
                       genome = genome,
                       genomeTicksDisplay = TRUE,
                       genomeTicksScale = 1e+5,
                       tracklist = link.list.sum,
                       genomeFillColor = c("white","white"),
                       genomeBorderColor = "black",
                       genomeLabelTextSize = "14pt",
                       genomeLabelOrientation = -180,
                       genomeLabelDy = 40,
                       genomeBorderSize = 1,
                       zoom = FALSE,
                       genomeTicksTextSize = 14,
                       chrPad = 0.002,
                       LINKMouseOverStrokeWidth = 2,
                       LINKMouseOverOpacity = 0.9,
                       LINKMouseOverStrokeColor = c("red","red"),
                       LINKMouseOverDisplay = TRUE)


}
