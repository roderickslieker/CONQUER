#' Main function to visualize the investigated SNPs
#' @param directory [[character]] Character string of the directory in which the SNP data and optional summary files are stored
#' @param SNPs [[vector]] Vector of SNP names that correspond to the SNP names provided in summarize. Can also be a single SNP name, but should always be a rs* number.
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{summarize("somedirectory","rs1558902")}
visualize <- function(directory, SNPs){
  #List all files in directory
  allFiles <- list.files(directory)
  if(identical(allFiles, character(0))){
    stop(sprintf("%s is not a valid path",directory))
  }
  sumFiles <- allFiles[grepl("CONQUER_SummaryX",allFiles)]

  res <- character(length(SNPs))
  for(i in 1:length(SNPs)){
    SNP <- SNPs[i]
    notInDirectory <- tryCatch({
      load(sprintf("%s/%s.RData",directory,SNP))
    },
    warning = function(cond){
      message(sprintf("%s was not found in %s",SNP, directory))
      return(i)
    })
    res[i] <- notInDirectory
  }
  SNPs <- res[grepl("rs",x = res)]

  SNPSummary <- NULL
  if(!identical(character(0),sumFiles)){
    load(paste0(directory,"/",sumFiles[1]))
  }
  loadedSNPs <- sapply(SNPs,function(SNP){get(SNP)},simplify = F)
  app <- visualizeDashboard(SNPs = loadedSNPs,SNPSummary = SNPSummary)
  shiny::runApp(app,launch.browser = T)
}









