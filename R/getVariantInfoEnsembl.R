#' getVariantInfoEnsembl
#'
#' @param rsID
#' @param population
#'
#' @keywords internal
#' @usage NULL
#'
#' @return [[list]]
getVariantInfoEnsembl <- function(rsID, population="1000GENOMES:phase_3:CEU"){
  URL <- paste0("http://rest.ensembl.org/variation/human/", rsID)
  SNP <- tryCatch({
    jsonlite::fromJSON(URL, flatten = T)},
    error = function(e) {
      errorNum <- gsub(".*?([0-9]+).*", "\\1", e)
      if(as.numeric(errorNum) == 503){
        warning(paste0("The Rest API of ensembl is currently unavailable.\
                       Will try to retrieve information about ",
                       rsID,
                       " later.")
                )
      }else{
        warning(paste0(rsID, " is skipped. It might not be a valid rs ID"))
      }

    }
  )
  variantInfo <- as.list(c(r2 = 1,
                           start = as.numeric(SNP$mappings$start[1]),
                           consequence_type = SNP$most_severe_consequence,
                           end = as.numeric(SNP$mappings$end[1]),
                           strand = as.numeric(SNP$mappings$strand[1]),
                           population_name = population,
                           variation = as.character(SNP$name),
                           d_prime = 1,
                           chr = SNP$mappings$seq_region_name[1],
                           clinical_significance = as.logical(NA)))

  variantInfo$clinical_significance <- as.logical(NA)

  return(variantInfo)
}
