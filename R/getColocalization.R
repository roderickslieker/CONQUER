#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom reshape2 colsplit
getColocalization <- function(rsID)
{
  cat(rsID,"\n")
  lead <- load(sprintf("./%s.RData", rsID))
  tempdata <- get(lead)
  rm(list = lead)
  
  LD <- tempdata$topHits
  eQTLs <- rbind(tempdata$eQTLs)
  eQTLs <- eQTLs[eQTLs$pValue <= eQTLs$pValueThreshold,] %>% na.omit()
  
  all.gencode <- unique(eQTLs$gencodeId)
  lg <- length(all.gencode)
  
  if(lg >= 1){
    data.qtls <- lapply(all.gencode, function(g){
      #cat(g)
      # Get QTL data
      query <- curl::curl_fetch_memory(sprintf("https://gtexportal.org/rest/v1/association/metasoft?gencodeId=%s&datasetId=gtex_v8", g))
      chars_with_Nan <- rawToChar(query$content)
      chars_with_NA <- gsub("NaN",'"NA"',chars_with_Nan)
      data.qtls.temp <- jsonlite::fromJSON(chars_with_NA)[[1]]
      Position <- reshape2::colsplit(data.qtls.temp[,"variantId"], "_", LETTERS[1:5])[,2]
      data.qtls.temp <- data.qtls.temp[which(Position %in% LD$start),]
      rm(Position)
      rownames(data.qtls.temp) <- NULL
      return(data.qtls.temp)
    }) 
    names(data.qtls) <- all.gencode
    
    allColoc <- lapply(names(data.qtls), calculateColocalization, data.qtls=data.qtls, LD=LD) %>% do.call(what=rbind)
    allColoc$Symbol <- eQTLs[match(allColoc$genecodeId, eQTLs$gencodeId),"gene"]
    allColoc$Lead <- ifelse(allColoc$snp == tempdata$SNP$variation, 1,NA)
    ntiss <- length(unique(allColoc$Tissue))
    return(allColoc)
    
    
  }else{
  
  }
}

