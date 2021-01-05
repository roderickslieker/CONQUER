#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom viridis viridis_pal
plotColoc <- function(rsID, all.coloc=ColocSummary, loadedSNPs=loadedSNPs, filter=FALSE, interactive=TRUE,tissues=NULL)
{
  snp.data <- loadedSNPs[[rsID]]
  all.coloc.sel <- all.coloc[[rsID]]



  if(!is.null(tissues))
  {
    all.coloc.sel <- all.coloc.sel[all.coloc.sel$Tissue %in% tissues,]
  }

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
      # Range
      start.coloc <- min(all.coloc.sel$Position)
      end.coloc <- max(all.coloc.sel$Position)
      coloc.range <- GRanges(paste0("chr", snp.data$SNP$chr), IRanges(start.coloc, end.coloc))

      # Extract genes
      genes <- snp.data$genes
      genes <- genes[as.matrix(findOverlaps(coloc.range, genes))[,2],]

      check.genes <- length(genes) == 0

      if(check.genes)
      {
        genes <- data.frame(start = NA,
                            end = NA,
                            strand = NA,
                            direction = NA,
                            molecule = "Gene",
                            gene = "No genes",
                            ENSEMBL = NA)
      }else{
        genes <- data.frame(start = start(genes), end = end(genes),
                            strand = ifelse(strand(genes) == "+", "forward", "reverse"),
                            direction = ifelse(strand(genes) == "+", 1, -1),
                            molecule = "Gene",
                            gene = elementMetadata(genes)$name,
                            ENSEMBL = elementMetadata(genes)$gene_id)
        genes$gene <- factor(genes$gene)

      }


      cols <- sample(c(viridis::viridis_pal(option = "A")(50), viridis::viridis_pal(option = "D")(50)), 50)

      #Label
      all.coloc.sel$Label <- sprintf("SNP:%s<br>Lead SNP:%s", all.coloc.sel$snp, all.coloc.sel$Lead.name)

      LD <- snp.data$LD
      LD <- LD[LD$start >= start.coloc & LD$start <= end.coloc,]
      Lead <- do.call(cbind, snp.data$SNP) %>% as.data.frame()

      p0 <-  ggplot(LD, aes(x=start, y=r2, col=r2,SNP=variation, consequence_type=consequence_type))+
        geom_point()+
        theme(legend.position = "none")+
        ylab("LD (r2)")+
        xlab("Position")+
        ylim(0,1)+
        xlim(start.coloc, end.coloc)+
        scale_colour_gradientn(colours = viridis::viridis_pal(option = "D")(10))+
        theme(legend.position = "none")+
        geom_point(data=Lead, aes(x=as.numeric(start), y=1), pch=8, col="black")


      p1 <-  ggplot(all.coloc.sel, aes(x=Position, y=SNP.PP, col=Tissue,label=Label))+
        geom_point()+
        geom_line(lwd=.5, alpha=.5)+
        theme(legend.position = "none")+
        facet_grid(Symbol~.)+
        ylab("Posterior probability of the SNP")+
        xlab("Position")+
        ylim(0,1)+
        xlim(start.coloc, end.coloc)+
        scale_colour_manual(values = cols)+
        theme(legend.position = "none")


      if(sum(!is.na(all.coloc.sel$Lead)) != 0)
      {
        p1 <- p1 + geom_point(aes(x=Position, y=Lead), pch=8, col="black")
      }

      p2 <- ggplot(genes, aes(label=gene)) +
        geom_linerange(aes(ymin = start, ymax = end, x = gene, col=gene), position = position_dodge(.5), lwd=2) +
        theme(legend.position = "none")+
        coord_flip(ylim=c(start.coloc, end.coloc))+
        scale_colour_manual(values = viridis::viridis_pal()(length(unique(genes$gene))))+
        ylab("Position")+
        xlab("")

      if(!check.genes)
      {
        if(sum(genes$strand == "forward")>=1){p2 <- p2+geom_text(data=genes[genes$strand == "forward",],
                                                                 aes(y=end, x=gene,size = 3, label = '>', col=gene),position = position_dodge(.5))}
        if(sum(genes$strand == "reverse")>=1){p2 <- p2 + geom_text(data=genes[genes$strand == "reverse",],
                                                                   aes(y=start, x=gene,size = 3, label = '<', col=gene),position = position_dodge(.5))}
      }

      out <- list(p1,p2,p0)
    }else{

      out <- ggplot2::ggplot(data.frame(x=1,y=1, label="Cannot test colocalization for this gene/SNP/tissue set"),
                                              aes(x=x,y=y, label=label))+
                                ggplot2::geom_text()+
                                ggplot2::theme_void()
    }
  }else{
    out <- ggplot2::ggplot(data.frame(x=1,y=1, label="Cannot test colocalization for this gene/SNP/tissue set"),
                                            aes(x=x,y=y, label=label))+
                              ggplot2::geom_text()+
                              ggplot2::theme_void()
  }

  if(interactive)
  {
    if(class(out) == "list")
    {
      out <- plotly::subplot(p1, p0,p2, nrows = 3, shareX = T, heights = c(0.6,0.25,0.15))
    }else{
      out <- plotly::ggplotly(out)
    }
  }else{
    if(class(out) == "list")
    {
      p1 <- out[[1]]
      p2 <- out[[2]]
      p0 <- out[[3]]
      p1 <- p1 + theme(legend.position = "bottom")
      out <- patchwork::wrap_plots(p1,p0, p2, nrow=3, heights =c(0.6,0.25,0.15))
    }else{
    }
  }
  return(out)
}
