#' Collect data for one or more SNPs and perform an optional integrative analysis on multiple tissues.
#' @param variants vector Vector with one or more SNP names (rs*)
#' @param precalculated Default is TRUE. GTEx by default only takes along SNPs 1MB from the
#' start site. If TRUE, only the precalculated SNPs by GTEx are used. When FALSE, CONQUER will
#' calculate all pairwise comparisons which may take substantially longer (up to days for large number of SNPs).
#' Note that this may yield more interesting results as GTEx may miss interesting QTLs.
#' @param multiAnalyze Default is FALSE. If TRUE, an integrated analysis will be performed across the SNPs.
#' @param tissues Tissues of interest. Should be a vector of names based on teh
#' @param directory String of the directory where the files should be stored for the SNPs investigated
#' @param token Token required for LDlink, which can be obtained from the LDlink website
#' @param population Letter code of population of interest for example " CEU"
#' @export
#' @examples \dontrun{
#' library(CONQUER)
#' summarize(variants = c("rs1558902","rs601945"),
#' directory=getwd(),
#' multiAnalyze=FALSE,
#' token="sometoken",
#' tissues=NULL)}

summarize <- function(variants, precalculated=TRUE, multiAnalyze=FALSE,
                      tissues = NULL ,directory=NULL, token=NULL, population="CEU") {
  if(is.null(tissues)) tissues <- conquer.db::gtexTissuesV8
  # Skip if existent
  if(is.null(token)) stop("Please provide a LDlink token!")

  filenames <- sprintf("%s/%s.RData",directory,variants)
  SNPsRemain <- variants[!file.exists(filenames)]
  cat(sprintf("%s variants are already present, %s to be done \n",sum(file.exists(filenames)),length(SNPsRemain)))

  if(length(SNPsRemain) != 0)
  {
    Chromatin <- c(Chromatin1,Chromatin2,Chromatin3)

    stats <- lapply(X = SNPsRemain,
                    FUN = getDataforSingleSNP,
                    directory = directory,
                    token = token,
                    population = population,
                    Chromatin=Chromatin,
                    allTissues = tissues,
                    precalculated = precalculated)
  }


  filenames <- sprintf("%s/%s.RData", directory, variants)
  SNPsRemain <- variants[!file.exists(filenames)]

  if(length(SNPsRemain) == 0)
  {
    allFiles <- list.files(directory)
    colocFiles <- allFiles[grepl("Colocalization_Summary",allFiles)]

    if(length(colocFiles) == 0){
      message("Colocalization has not yet been performed. Note that this may take time some time depending on the number of SNPs. Running..")
      all.coloc <- lapply(variants, getColocalization)
      names(all.coloc) <- variants
      save(all.coloc, file=paste0(directory, "/", "Colocalization_Summary.RData"))
    }

    summFiles <- allFiles[grepl("CONQUER_Summary",allFiles)]

    if(length(summFiles) == 0 & multiAnalyze){
      message("MultiAnalyze is true. The SNPs will be analyzed for the following tissues:")
      lapply(tissues,message)
      SNPSummary <- abstractAnalyze(variants = variants, directory = directory, tissues = tissues, clustering = "agnes")
      filename <- sprintf("CONQUER_Summary%s.RData", gsub("[.]","", make.names(Sys.time())))
      save(SNPSummary, file = paste0(directory, "/", filename))
      message(sprintf("CONQUER SNP summary saved in %s, with the following name %s", directory, filename))
    }else{
      message("Completed, you can now run visualize.")
    }
  }
  #return(outputLog)
}
