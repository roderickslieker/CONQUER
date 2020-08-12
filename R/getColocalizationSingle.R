getColocalizationSingle <- function(gencodeId, leadSNP, tissue, loadedSNPs)
{
  SNP.sel <- loadedSNPs[[leadSNP]]
  LD <- SNP.sel$LD
  SNPs <- LD[LD$r2 >= 0.8,"variation"]
  qtls.data <- CONQUER:::get_eQTL_bulk(gencodeId, SNPs, tissue)

  if(nrow(qtls.data) >= 1)
  {
    abf.data <- coloc::finemap.abf(dataset=list(beta = as.numeric(qtls.data$nes), sdY=1,
                                                varbeta = as.numeric(qtls.data$error)^2, snp=qtls.data$snpId, type="quant"))


    abf.data <- na.omit(abf.data)
    LD.sel <- LD[match(abf.data$snp, LD$variation),]
    abf.data$Tissue <- tissue
    abf.data$Position <- LD.sel$start
    abf.data$Symbol <- gencodeId
    abf.data$Lead <- ifelse(abf.data$snp == leadSNP,1,NA)
    abf.data <- list(abf.data)
    names(abf.data) <- leadSNP
    out <- plotColoc(rsID = leadSNP, all.coloc=abf.data,
                     loadedSNPs = loadedSNPs, filter=FALSE)
  }else{
    out <- plotly::ggplotly(ggplot2::ggplot(data.frame(x=1,y=1, label="Cannot test colocalization for this gene/SNP/tissue set"),
                                   aes(x=x,y=y, label=label))+
                              ggplot2::geom_text()+
                              ggplot2::theme_void())
  }
  return(out)
}


