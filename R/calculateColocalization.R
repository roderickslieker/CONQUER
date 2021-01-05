#' @importFrom reshape2 colsplit
#' @importFrom coloc finemap.abf
calculateColocalization <- function(genecode.id, data.qtls=data.qtls, LD=LD)
{
  #cat(genecode.id)
  data.in <- data.qtls[[genecode.id]]
  if(nrow(data.in) != 0){
    vid <- data.in[,"variantId"]
    vid <- reshape2::colsplit(vid, "_", names=c("Chromsome","Position","REF","ALT","Build"))
    vid <- data.frame(vid, rsID = data.in[,"SNP"])
    max.in <- grep("datasetId", colnames(data.in)) - 1
    tissue.eqtls <- data.in[,1:max.in]

    newData <- lapply(colnames(tissue.eqtls), function(tissue)
    {
      #cat(tissue)
      td.sub <- tissue.eqtls[,tissue]
      td.sub$Tissue <- tissue
      td.sub <- data.frame(td.sub, vid)
      #td.sub <- td[td$Position %in% LD$start,]
      #td.sub$rsID <-LD[match(td.sub$Position, LD$start),"variation"]
      td.sub <- na.omit(td.sub)

      if(nrow(td.sub) !=0)
      {
        abf.data <- coloc::finemap.abf(dataset=list(pvalues = td.sub$pValue, beta = td.sub$nes, sdY=1,
                                                    varbeta = td.sub$se^2, snp=td.sub$rsID, type="quant"))

        abf.data <- na.omit(abf.data)
        td.subx <- td.sub[match(abf.data$snp, td.sub$rsID),]
        abf.data$NES <- td.subx$nes
        abf.data$pValue <- td.subx$pValue

        abf.data$SE <- td.subx$se
        abf.data$Position <- td.subx$Position


        abf.data$Tissue <- tissue
        abf.data$genecodeId <- genecode.id
        return(abf.data)

      }else{}
    }) %>% do.call(what=rbind)
    return(newData)
  }else{
    NULL
  }
}
