#' AnalyseSNPs
#' @param eQTLData [[data.frame]] A data.frame containing the eQTL data.
#' @param tissue [[character]] A character vector containing the tissues of interest.
#' @param clustering [[character]] A character string indicating which clustring method should be used: `agnes` or `PAM`
#' @keywords internal
#' @import dplyr
#' @return [[list]]
AnalyseSNPs <- function(eQTLData, tissue, clustering){
  assign("expressionData",get(tissue, pos=2))
  QTLs <- eQTLData$gencodeId %>% unique()
  expressedGenes <- expressionData %>% colnames()
  SNPsnr <- eQTLData$SNP %>% unique() %>% length()
  message(paste0("Counted ", SNPsnr, " SNPs and ", length(QTLs),
                 " eQTLs for ",
                 eQTLData[1,"tissue"], " of which ",
                 intersect(expressedGenes,QTLs) %>% length(),
                 " eQTLs are expressed in the data."))

  QTLs <- QTLs[QTLs %in% expressedGenes]

  #Get Co-expressed genes
  coQTLs <- getCoexpression(geneSet = QTLs,expressionData = expressionData ,threshold = 0.90)

  message(sprintf("Found %s co-expressed genes",length(coQTLs)))
  #Cluster QTLS + Co-expressed genes.
  message(paste0("Clustering with ", clustering))
  if(clustering == "PAM"){
    cluster_info <- ClusterPam(expressionData[,coQTLs])
  }else if(clustering == "agnes"){
    cluster_info <- ClusterAgnesSpearman(expressionData[,coQTLs])
  }

  #Extract cluster assignment from output
  modules <- data.frame("gene" = names(cluster_info$cluster),
                        "cluster" = cluster_info$cluster)

  #modules <- wgcna_clustering(expressionData[,coQTLs])

  Module_list <- split(x = modules, f = modules$cluster)

  Module_list <- sapply(X = Module_list,
                        FUN = function(mod){
                          if((mod[,1] %in% QTLs) %>% sum() == 0 ){
                            return(NULL)
                          }
                          return(mod[,1])
                        })
  #Removing modules with zero QTLs
  Module_list <- Filter(Negate(is.null),Module_list)


  #Determine QTLs in Modules
  Module_QTLs <- sapply(X = Module_list,
                        FUN = function(mod){
                          return(mod[mod %in% QTLs])
                        })
  #Retrieve statistics from original  data
  Module_SNPs <- lapply(X = Module_QTLs,
                        FUN = function(modQTLs){
                          SNPinfo <- eQTLData[eQTLData$gencodeId %in% modQTLs,] %>% dplyr::arrange(by=Pval.ratio)
                          return(SNPinfo)
                        })

  #Score each module based on the Pval.Ratio
  Module_ratio_score <- lapply(X = Module_SNPs,
                               FUN = function(modSNPs){
                                 modSNPs$Pval.ratio %>% log(base = 10) %>% sum()
                               })



  out <- list("Module_Score" = Module_ratio_score,
              "Module_SNPs_eQTLs" = Module_SNPs,
              "Module_Genes" = Module_list,
              "Expression" = expressionData[,coQTLs])
  return(out)
}
