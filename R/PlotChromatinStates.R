#' @importFrom plotly ggplotly
#' @import ggplot2
PlotChromatinStates <- function(StatesPlotData, fillScale, jit, interactive=TRUE){

  ###Plot###
  chromatinplot <- ggplot2::ggplot(data = StatesPlotData, ggplot2::aes(x=start,y=Y)) +
    ggplot2::geom_segment(ggplot2::aes(x=-Inf,xend=Inf, y=sample, yend=sample),linetype=3, color="black") +
    ggplot2::geom_rect(ggplot2::aes(fill=target, xmin=start, xmax=end),
                       ymin=StatesPlotData$Y-0.4, ymax=StatesPlotData$Y+0.4, alpha=0.7) +
    ggplot2::geom_point(data = jit, ggplot2::aes(x = value, y=index), size=1, shape=3, alpha=0.5) +
    ggplot2::scale_fill_manual(values = fillScale) + ggplot2::theme_minimal()

  if(interactive)
  {
    plotly::ggplotly(chromatinplot)
  }else{
    chromatinplot
  }

}

