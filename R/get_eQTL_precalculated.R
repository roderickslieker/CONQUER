#' get_eQTL_bulk
#' @param gene.in
#' @param lead
#' @param tissues
#' @keywords internal
#' @return [[data.frame]]
get_eQTL_precalculated <- function(lead, tissues){
  query <- curl::curl_fetch_memory(sprintf("https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json&snpId=%s&datasetId=gtex_v8", lead))
  chars_with_Nan <- rawToChar(query$content)
  chars_with_NA <- gsub("NaN",'"NA"',chars_with_Nan)
  data.qtls.temp <- jsonlite::fromJSON(chars_with_NA)[[1]]

  data.qtls.temp <- data.qtls.temp[data.qtls.temp$tissueSiteDetailId %in% tissues,]

  if(length(data.qtls.temp) != 0)
  {


    outdata <- data.frame(datasetId = data.qtls.temp$datasetId[1],
                          error = NA,
                          gencodeId = data.qtls.temp$gencodeId,
                          geneSymbol = data.qtls.temp$geneSymbol,
                          message = NA,
                          nes = data.qtls.temp$nes,
                          pValue = data.qtls.temp$pValue,
                          pValueThreshold = 0.001,
                          snpId = data.qtls.temp$snpId,
                          tStatistic = NA,
                          tissueSiteDetailId = data.qtls.temp$tissueSiteDetailId,
                          variantId = data.qtls.temp$variantId)


    if(length(outdata) == 0){
    }else if(nrow(outdata) == 0){
    }else{
      outdata$pValue <- as.numeric(outdata$pValue)
      return(outdata)
    }


  }

}
