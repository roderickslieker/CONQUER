#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom viridis viridis_pal
plotColoc <- function(rsID, all.coloc=ColocSummary, loadedSNPs=loadedSNPs, filter=TRUE)
{
  snp.data <- loadedSNPs[[rsID]]
  all.coloc.sel <- all.coloc[[rsID]]
  if(!is.null(all.coloc.sel))
  {
    all.coloc.sel$Id <- paste0(all.coloc.sel$Tissue,"_", all.coloc.sel$genecodeId)
    maxx <- by(all.coloc.sel$SNP.PP, all.coloc.sel$Id, max) %>% as.list() %>% do.call(what=c)
    if(filter)
    {
      all.coloc.sel <- all.coloc.sel[all.coloc.sel$Id %in% names(maxx[maxx >= .10]),]
    }

    if(nrow(all.coloc.sel) >= 1)
    {
      all.coloc.sel$Lead.name <- rsID
      #all.coloc.sel <- all.coloc.sel[all.coloc.sel$SNP.PP >= 0.1,]
      # Range
      start.coloc <- min(all.coloc.sel$Position)
      end.coloc <- max(all.coloc.sel$Position)
      coloc.range <- GRanges(paste0("chr", snp.data$SNP$chr), IRanges(start.coloc, end.coloc))

      # Extract genes
      genes <- snp.data$genes
      genes <- genes[as.matrix(findOverlaps(coloc.range, genes))[,2],]
      genes <- data.frame(start = start(genes), end = end(genes),
                          strand = ifelse(strand(genes) == "+", "forward", "reverse"),
                          direction = ifelse(strand(genes) == "+", 1, -1),
                          molecule = "Gene",
                          gene = elementMetadata(genes)$name,
                          ENSEMBL = elementMetadata(genes)$gene_id)
      genes$gene <- factor(genes$gene)

      cols <- sample(c(viridis::viridis_pal(option = "A")(50), viridis::viridis_pal(option = "D")(50)), 50)

      p1 <-  ggplot(all.coloc.sel, aes(x=Position, y=SNP.PP, col=Tissue,label=Lead.name))+
        geom_point()+
        geom_line(lwd=.5, alpha=.5)+
        theme(legend.position = "none")+
        facet_grid(Symbol~.)+
        ylab("Posterior probability of the SNP")+
        xlab("Position")+
        ylim(0,1)+
        xlim(start.coloc, end.coloc)+
        scale_colour_manual(values = cols)+
        theme(legend.position = "none")+
        geom_point(aes(x=Position, y=Lead), pch=8, col="black")


      p2 <- ggplot(genes, aes(label=gene)) +
        geom_linerange(aes(ymin = start, ymax = end, x = gene, col=gene), position = position_dodge(.5), lwd=2) +
        theme(legend.position = "none")+
        coord_flip(ylim=c(start.coloc, end.coloc))+
        scale_colour_manual(values = viridis::viridis_pal()(length(unique(genes$gene))))+
        ylab("Position")+
        xlab("")

      if(sum(genes$strand == "forward")>=1){p2 <- p2+geom_text(data=genes[genes$strand == "forward",],
                                                               aes(y=end, x=gene,size = 3, label = '>', family = "Arial", col=gene),position = position_dodge(.5))}
      if(sum(genes$strand == "reverse")>=1){p2 <- p2 + geom_text(data=genes[genes$strand == "reverse",],
                                                                 aes(y=start, x=gene,size = 3, label = '<', family = "Arial", col=gene),position = position_dodge(.5))}

      out <- plotly::subplot(p1,p2, nrows = 2, shareX = T, heights = c(0.8,0.2))
    }else{

      out <- plotly::ggplotly(ggplot2::ggplot(data.frame(x=1,y=1, label="Cannot test colocalization for this gene/SNP/tissue set"),
                                              aes(x=x,y=y, label=label))+
                                ggplot2::geom_text()+
                                ggplot2::theme_void())
    }
  }else{
    out <- plotly::ggplotly(ggplot2::ggplot(data.frame(x=1,y=1, label="Cannot test colocalization for this gene/SNP/tissue set"),
                                            aes(x=x,y=y, label=label))+
                              ggplot2::geom_text()+
                              ggplot2::theme_void())
  }
  return(out)
}
