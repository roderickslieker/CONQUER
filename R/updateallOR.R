#' Add eQTLs to allOR
updateallOR <- function(allOR.in, eQTLs.in)
{
  eQTLs.in <- eQTLs.in[eQTLs.in$pValue <= 0.05,]
  for(t in 1:length(allOR.in))
  {
    tissue <- names(allOR.in)[t]
    eQTLs.sub <- eQTLs.in[eQTLs.in$tissue %in% tissue,]
    allOR.sub <- allOR.in[[tissue]]
    for(i in 1:length(allOR.sub))
    {
      temp <- allOR.sub[[i]]
      if(is.null(temp)) next
      temp$eGenes <- NA
      temp$eQTLs <- NA

      for(k in 1:nrow(temp))
      {
        gns <- strsplit(temp[k,"Genes"],";")[[1]]
        gns.mod <- gns[gns %in% eQTLs.sub$gene]
        if(length(gns.mod) >= 1)
        {
          gns.mod <- paste0(gns.mod, collapse = ";")
          eQTLs.sub.mod <- eQTLs.sub[eQTLs.sub$gene %in% gns,]
          snp.mod <- paste0(unique(eQTLs.sub.mod$SNP), collapse = ";")
          temp[k,"eGenes"] <- gns.mod
          temp[k,"eQTLs"] <- snp.mod
        }else{
          next
        }

      }
      allOR.sub[[i]] <- temp
    }
    allOR.in[[tissue]] <- allOR.sub
  }
  return(allOR.in)
}

