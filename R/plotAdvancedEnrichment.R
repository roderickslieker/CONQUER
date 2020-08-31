plotAdvancedEnrichment <- function(tissue, interactive=TRUE, SNPSummary)
{
  pws <- SNPSummary[["allOR",tissue]]
  for(k in names(pws)){
    temp <- pws[[k]]
    if(is.null(temp)) next
    if(nrow(temp) >= 1){
      temp$Module <- paste0("Module",k)
    }
    pws[[k]] <- temp
  }
  pws <- do.call(rbind, pws)
  pws$Pvalue <- -log10(pws$Adjusted.P.value)
  px <- ggplot(pws, aes(x=Odds.Ratio, y=Pvalue, col=DB, label=Term, Overlap=Overlap,
                        Pvalue=Adjusted.P.value, Module = Module))+
    geom_point()+
    scale_x_continuous(trans="log2")+
    theme(legend.position="none")+
    scale_colour_manual(values = viridis::viridis_pal()(length(levels(pws$DB))))+
    ylab("Adjusted P-value")+
    xlab("Odds ratio")

  if(interactive)
  {
    plotly::ggplotly(px)
  }else{
    px
  }

}
