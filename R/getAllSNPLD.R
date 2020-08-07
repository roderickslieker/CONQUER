#' LDsnps
getAllSNPLD <- function(SNPs){
      LDSNPs <- lapply(names(SNPs), function(SNP){
        LDtable <- SNPs[[SNP]]$LD
        LDSNPs <- LDtable[LDtable$r2 >= 0.8,c("variation","chr","start")]
        colnames(LDSNPs)[1] <- "LDSNP"
        LDSNPs$leadingSNP <- SNP
        return(LDSNPs)
      }) %>% do.call(what=rbind)
      LDSNPs$id <- paste0("chr", LDSNPs$chr,"_", LDSNPs$start)
      return(LDSNPs)
}
