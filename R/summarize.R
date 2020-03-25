#' Collect data for one or more SNPs and perform an optional integrative analysis on multiple tissues. 
#' @param variants [[vector]] Vector with one or more SNP names (rs*)
#' @param multiAnalyze [[boolean]] Default is FALSE. If TRUE, an integrated analysis will be performed across the SNPs. 
#' @param tissues [[vector]] Tissues of interest. Should be a vector of names based on teh 
#' @param directory [[character]] String of the directory where the files should be stored for the SNPs investigated
#' @param token [[character]] Token required for LDlink, which can be obtained from the LDlink website
#' @param population [[character]] Letter code of population of interest for example " CEU"
#' @export
#' @examples
summarize <- function(variants, multiAnalyze=FALSE, tissues ,directory=NULL, token=NULL, population="CEU") {
  # Skip if existent
  filenames <- sprintf("%s/%s.RData",directory,variants)
  SNPsRemain <- variants[!file.exists(filenames)]
  cat(sprintf("%s variants are already present, %s to be done \n",sum(file.exists(filenames)),length(SNPsRemain)))
  #if(parallel){
  #  BiocParallel::register(BiocParallel::bpstart(BiocParallel::SnowParam(cores)))
  #}
  stats <- lapply(X = SNPsRemain,
         FUN = getDataforSingleSNP,
         directory = directory,
         token = token,
         population = population)

  if(multiAnalyze){
    message("MultiAnalyze is true. The SNPs will be analyzed for the following tissues:")
    lapply(tissues,message)
    SNPSummary <- abstractAnalyze(variants = variants, directory = directory, tissues = tissues, clustering = "agnes")
    filename <- sprintf("CONQUER_Summary%s.RData", gsub("[.]","", make.names(Sys.time())))
    save(SNPSummary, file = paste0(directory, "/", filename))
    message(sprintf("CONQUER SNP summary saved in %s, with the following name %s", directory, filename))
  }
  outputLog <- do.call(rbind,stats) %>% data.frame()
  return(outputLog)
}
