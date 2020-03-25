#' getLDLink
#'
#' @param SNP
#' @param token
#' @param popilation
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[data.frame]]
getLDLink <- function(SNP, token, population="CEU") {
  rsID <- SNP$variation
  #make query link
  url <- sprintf("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=%s&pop=%s&r2_d=r2&token=%s",rsID, population, token)
  #get data from API
  response <- httr::GET(url = url)
  #check for errors
  httr::stop_for_status(response$status_code, task = paste("get LD for", rsID))
  #make response into text format
  LD <- httr::content(response, encoding = "UTF-8")
  if(length(grep("error", names(LD))) == 0)
  {
    LD <- list(LD, error = NULL)
  }
  if(length(grep("monoallelic",LD$error)) >= 1 |
     length(grep("not in 1000G reference panel", LD$error)) >= 1 |
     length(grep("not a biallelic variant.", LD$error)) >= 1)
  {
    cat(sprintf("Reverting to ENSEMBL:%s",LD$error))
    LD <- as.data.frame(jsonlite::fromJSON(paste0("http://rest.ensembl.org/ld/human/",
                                                  SNP$variation, "/",
                                                  paste0("1000GENOMES:phase_3:",SNP$population_name), "?attribs=1"),
                                           flatten = T)) #bij rs56116432 kwam lege lijst terug
    #add query variant
    finalData <- rbind(LD, SNP)
    #make r2, end and start columns numeric
    finalData[,c("r2", "start", "end")] <- apply(finalData[,c("r2", "start", "end")],
                                                 2, function(x) as.numeric(as.character(x)))

  }else{
    LD <- httr::content(response,type = "text", encoding = "UTF-8")
    #make text into dataframe
    LD <- as.data.frame(readr::read_tsv(LD))
    #remove rows without rsID
    LD <- LD[LD$RS_Number!=".",]
    #take region around high LD block, with block in middle
    highLD <- LD[LD$R2>=0.8,]
    block.radius <- max(abs(c(min(highLD$Distance), max(highLD$Distance))))+10000
    block.center <- median(c(as.numeric(SNP$start)+min(highLD$Distance),
                             as.numeric(SNP$start)+max(highLD$Distance)))
    start.r <- floor(block.center - block.radius)
    end.r <- ceiling(block.center + block.radius)
    region <- sprintf("%s:%s-%s",SNP$chr, start.r, end.r)

    #get hg38 positions and consequence types
    url <- sprintf("http://rest.ensembl.org/overlap/region/human/%s?feature=variation", region)
    r <- httr::GET(url = url, httr::content_type("application/json"))
    httr::stop_for_status(r, task = "get consequence types of LD SNPs")


    allSNPs <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)), simplifyDataFrame = TRUE, flatten = TRUE)
    allSNPs <- allSNPs[c('id', 'start','consequence_type', "clinical_significance")]

    out <- as.list(rep(NA, 4))
    for(i in 1:ncol(allSNPs)){
      if(i != 4)
      {
        out[[i]] <- do.call(c,allSNPs[,i])
      }else{
        temp <- allSNPs[,i]
        for(k in 1:length(temp)){
          if(length(temp[[k]])==0){temp[[k]] <- NA
          }else if(nrow(temp[[k]]) ==1){
          }else if(nrow(temp[[k]]) ==2){
            temp[[k]] <- sprintf("%s/%s",temp[[k]][1,1],temp[[k]][2,1])
          }
        }
        out[[i]] <- as.character(temp)

        #        idx <- which(!isEmpty(temp))
        #       data.in <- do.call(c, temp) %>%
        #        do.call(what = c)
        #     data.na <- rep(NA, length(temp))
        #    data.na[idx] <- data.in
        #   allSNPs[[i]] <- data.na

      }

    }
    allSNPs <- do.call(cbind, out) %>% as.data.frame(stringsAsFactors=F)
    colnames(allSNPs) <- c('id', 'start','consequence_type', "clinical_significance")
    allSNPs[,2 ] <- as.numeric(allSNPs[,2])

    #keep only SNPs with LD data
    allSNPs <- allSNPs[is.element(allSNPs$id, LD$RS_Number),]
    #merge LD info and hg38+consequence types
    finalData <- merge.data.frame(allSNPs, LD, by.x = 'id', by.y = 'RS_Number')
    finalData <- finalData[,c("id", "start", "consequence_type",
                              "R2", "clinical_significance")]
    colnames(finalData) <- c("variation", "start", "consequence_type",
                             "r2", "clinical_significance")

    finalData$chr <- SNP$chr
  }

  return(finalData)
}
