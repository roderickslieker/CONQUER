#' get_eQTL_bulk
#' @param genesx
#' @param lead
#' @param tissues
#' @keywords internal
#' @importFrom purrr cross_df
#' @return [[data.frame]]
get_eQTL_bulk <- function(genesx, lead, tissues){
  tmp_list <- list("gencodeId" = genesx, "tissueSiteDetailId" = tissues, "variantId" = lead, "datasetId" = "gtex_v8")
  all_combinations_df <- purrr::cross_df(tmp_list)
  all_combinations_list <- split(all_combinations_df, seq(nrow(all_combinations_df)))
  all.comb <- lapply(all_combinations_list, jsonlite::toJSON)


  outdata <- lapply(1:length(all.comb), function(i){
    raw <- httr::POST("https://gtexportal.org/rest/v1/association/dyneqtl",
                      encode = "json",
                      body = all.comb[[i]])

    if(raw$status_code != 200)
    {}else{
      temp <- httr::content(raw,"text", encoding = "UTF-8")
      temp <- gsub("NaN",'\\"NA\\"',temp)
      parsed <- fromJSON(temp)
      parsed <- parsed$result
      #res <- do.call(cbind, parsed)
      if(ncol(parsed) == 10)
      {
        parsed <- data.frame(parsed[,1:6,drop=F],pValueThreshold=NA, parsed[,7:ncol(parsed),drop=F]) %>% as.matrix()
      }
      return(parsed)
    }


  }) %>% do.call(what=rbind) %>% data.frame(stringsAsFactors=FALSE)
  outdata$pValue <- as.numeric(outdata$pValue)
  outdata$pValueThreshold <- as.numeric(outdata$pValueThreshold)
  outdata <- outdata[!is.na(outdata$pValueThreshold),]
  return(outdata)
}
