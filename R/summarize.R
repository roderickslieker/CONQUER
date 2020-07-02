#' Collect data for one or more SNPs and perform an optional integrative analysis on multiple tissues.
#' @param variants vector Vector with one or more SNP names (rs*)
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

summarize <- function(variants, multiAnalyze=FALSE, tissues ,directory=NULL, token=NULL, population="CEU") {
  # Skip if existent
  filenames <- sprintf("%s/%s.RData",directory,variants)
  SNPsRemain <- variants[!file.exists(filenames)]
  cat(sprintf("%s variants are already present, %s to be done \n",sum(file.exists(filenames)),length(SNPsRemain)))
  #if(parallel){
  #  BiocParallel::register(BiocParallel::bpstart(BiocParallel::SnowParam(cores)))
  #}

  # Chromatin states

  if(length(SNPsRemain) != 0)
  {
    Chromatin <- c(conquer.db::Chromatin1,conquer.db::Chromatin2,conquer.db::Chromatin3)

    stats <- lapply(X = SNPsRemain,
                    FUN = getDataforSingleSNP,
                    directory = directory,
                    token = token,
                    population = population,
                    Chromatin=Chromatin)
  }


  if(multiAnalyze){
    message("MultiAnalyze is true. The SNPs will be analyzed for the following tissues:")
    lapply(tissues,message)
    SNPSummary <- abstractAnalyze(variants = variants, directory = directory, tissues = tissues, clustering = "agnes")
    filename <- sprintf("CONQUER_Summary%s.RData", gsub("[.]","", make.names(Sys.time())))
    save(SNPSummary, file = paste0(directory, "/", filename))
    message(sprintf("CONQUER SNP summary saved in %s, with the following name %s", directory, filename))
  }
  if(length(SNPsRemain) != 0)
  {
    if(exists(stats))
    {
      outputLog <- do.call(rbind,stats) %>% data.frame()
    }

  }else{
    message("Completed, you can now run visualize.")
  }
  return(outputLog)
}
