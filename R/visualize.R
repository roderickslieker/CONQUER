#' Main function to visualize the investigated SNPs
#' @param directory [[character]] Character string of the directory in which the SNP data and optional summary files are stored
#' @param SNPs [[vector]] Vector of SNP names that correspond to the SNP names provided in summarize. Can also be a single SNP name, but should always be a rs* number.
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{summarize("somedirectory","rs1558902")}
visualize <- function(directory, SNPs, tissues=NULL){
  #List all files in directory
  allFiles <- list.files(directory)
  if(identical(allFiles, character(0))){
    stop(sprintf("%s is not a valid path",directory))
  }

  sumFiles <- allFiles[grepl("CONQUER_SummaryX",allFiles)]
  colocFiles <- allFiles[grepl("Colocalization_Summary",allFiles)]

  cat("..Loading SNP data.....","\n")

  loadedSNPs <- lapply(SNPs, function(SNP){
    filedir <- sprintf("%s/%s.RData", directory, SNP)
    check <- !file.exists(filedir)
    if(check)
    {
      stop(sprintf("%s was not found in %s. Run summarize again for this SNP.",SNP, directory))
    }else{
      tmp <- load(filedir)
      out <- get(tmp)
      rm(list=tmp, envir = .GlobalEnv)
      return(out)
    }
  })

  names(loadedSNPs) <- SNPs


  SNPSummary <- NULL
  if(!identical(character(0),sumFiles)){
    load(paste0(directory,"/",sumFiles[1]))
  }

  if(!identical(character(0),colocFiles)){
    load(paste0(directory,"/",colocFiles[1]))
  }else{
    stop("The colocalization file is missing in the directory provided. Please run summarize again to generate this object or add this file to your directory with SNP files")
  }

  if(!is.null(tissues))
  {
    SNPSummary <- SNPSummary[,tissues]
  }


  app <- visualizeDashboard(loadedSNPs = loadedSNPs,SNPSummary = SNPSummary, ColocSummary = all.coloc)
  shiny::runApp(app,launch.browser = T)
}









