#' @import stringr
#' @import enrichR
#' @importFrom dplyr bind_rows first
abstractAnalyze <- function(variants, directory, tissues, clustering = "PAM", pcutoff){
  abstractData <- lapply(variants,function(x){
    load(paste0(directory,"/",x,".RData"))
    get(x)
  })

  options(stringsAsFactors = F)
  #Load extended kegg pathways
  message("Loading Pathways...")
  #Load canonical kegg pathways
  KEGG_DATA <- prepare_KEGG(species = "hsa",
                                              KEGG_Type = "KEGG",
                                              keyType = "kegg")
  #Map entrez to ENSG IDs
  KEGG_DATA$ENSG <- lapply(X = KEGG_DATA$PATHID2EXTID,
                           FUN = entrezToENSEMBL)

  #Extract QTLs from abstractData
  message("Loading eQTLs...")
  CiseQTL <- lapply(1:length(abstractData),function(i){
   if(nrow(abstractData[[i]]$eQTLs) != 0) return(abstractData[[i]]$eQTLs)
    }) %>% do.call(what=rbind)

  TranseQTL <- lapply(1:length(abstractData),function(i){
    if(nrow(abstractData[[i]]$eQTLsTrans) != 0) return(abstractData[[i]]$eQTLsTrans)
  }) %>% do.call(what=rbind)

  eQTLs <- dplyr::bind_rows(CiseQTL,TranseQTL); rm(CiseQTL,TranseQTL)

  eQTLs$gencodeId <- sapply(X = eQTLs$gencodeId,
                            FUN = function(ID){
                              stringr::str_split(ID, "[.]", simplify = TRUE) %>% dplyr::first()
                            }
  )

  #Rename tissues for eQTL data
  eQTLTissue <- gsub("[.]+", "_", tissues, perl=T)
  eQTLTissue <- ifelse(stringr::str_sub(eQTLTissue,-1) == "_", substr(eQTLTissue,start = 1, stop = str_length(eQTLTissue)-1),eQTLTissue)

  eQTLs <- eQTLs[eQTLs$tissue %in% eQTLTissue,]

  if(pcutoff == 'stringent')
  {
    eQTLs <- eQTLs[eQTLs$pValue <= eQTLs$pValueThreshold,]
  }else if(pcutoff == 'liberal'){
    eQTLs <- eQTLs[eQTLs$pValue <= 0.001 & eQTLs$Pval.ratio <= 2,]
  }else if(pcutoff == 'veryliberal')
  {
    eQTLs <- eQTLs[eQTLs$pValue <= 0.05,]
  }else{
    stop("pcutoff should be one of the following: stringent, liberal, veryliberal")
  }

  eQTLs_list <- split(x = eQTLs, f = eQTLs$tissue)
  #tissues <- sort(tissues)
  eQTLs_list <- eQTLs_list[eQTLTissue]
  eQTLs_list <- eQTLs_list[!is.na(names(eQTLs_list))]
  SNP_summary <- mapply(AnalyseSNPs, eQTLs_list, tissues, clustering)

  message("Calculating OR for modules with canonical kegg pathways")

  #Calculate OR for modules with  canonical kegg pathways
  canOR <- lapply(X = SNP_summary["Module_Genes",],
                  FUN = getOdds,
                  pathways = KEGG_DATA$ENSG,
                  annotation = KEGG_DATA$PATHID2NAME)

  allOR <- lapply(colnames(SNP_summary), getEnrichREnrichment, SNPSummary=SNP_summary)
  names(allOR) <- colnames(SNP_summary)

  #Remove modules that are not enriched
  canOR <- lapply(canOR, function(x){Filter(Negate(function(x)nrow(x) == 0), x)})
  allOR <- lapply(allOR, function(x){Filter(Negate(function(x)nrow(x) == 0), x)})
  allOR <- updateallOR(allOR.in = allOR, eQTLs.in = eQTLs)

  SNP_summary <- rbind(SNP_summary, canOR, allOR)
  return(SNP_summary)
}
