#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom reshape2 colsplit
#' @import magrittr
getColocalization <- function(rsID)
{
  cat(rsID,"\n")
  lead <- load(sprintf("./%s.RData", rsID))
  tempdata <- get(lead)
  rm(list = lead)

  LD <- tempdata$topHits
  #Covert names

  converted <- lapply(LD$variation, function(i){
    url <- sprintf("https://gtexportal.org/rest/v1/dataset/variant?format=json&snpId=%s&datasetId=gtex_v8", i)
    temp <- readLines(url)
    temp <- gsub("NaN",'"NA"',temp)
    fromJSON(temp)[[1]]
  }) %>% do.call(what=rbind)

  LD$variantID <- converted[match(LD$variation, converted$snpId),"variantId"]

  eQTLs <- tempdata$eQTLs
  eQTLs <- eQTLs[eQTLs$pValue <= 0.001,]
  eQTLs <- eQTLs[!is.na(eQTLs$nes),]
  all.gencode <- unique(eQTLs$gencodeId)

  lg <- length(all.gencode)

  if(lg >= 1){
    data.qtls <- lapply(all.gencode, function(g){
      # Get QTL data
      query <- curl::curl_fetch_memory(sprintf("https://gtexportal.org/rest/v1/association/metasoft?gencodeId=%s&datasetId=gtex_v8", g))
      chars_with_Nan <- rawToChar(query$content)
      chars_with_NA <- gsub("NaN",'"NA"',chars_with_Nan)
      data.qtls.temp <- jsonlite::fromJSON(chars_with_NA)[[1]]

      data.qtls.temp <- data.qtls.temp[which(data.qtls.temp$variantId %in% LD$variantID),]
      data.qtls.temp$SNP <- LD[match(data.qtls.temp$variantId, LD$variantID),1]
      rownames(data.qtls.temp) <- NULL
      return(data.qtls.temp)
    })
    names(data.qtls) <- all.gencode

    allColoc <- lapply(names(data.qtls), calculateColocalization, data.qtls=data.qtls, LD=LD) %>% do.call(what=rbind)
    allColoc$Symbol <- eQTLs[match(allColoc$genecodeId, eQTLs$gencodeId),"gene"]


    allColoc$Lead <- ifelse(allColoc$snp == tempdata$SNP$variation, 1,NA)
    return(allColoc)


  }
}

