getOdds <- function(genSets,pathways,annotation){
  background <- pathways %>% unlist() %>% unique() %>% length()
  ModuleOdds <- lapply(X = genSets,
                       FUN = function(moduleGenes){
                         Odds <- lapply(X = pathways,
                                        FUN = calculateOdds,
                                        moduleGenes,
                                        background)
                         Odds <- do.call(rbind,Odds) %>% as.data.frame()
                         colnames(Odds) <- c("OR","low_CI","up_CI","bg","P.Val","intersect")
                         Odds$PathwayName<- annotation[match(rownames(Odds),names(annotation))]
                         Odds <- Odds[order(Odds$P.Val),]
                         Odds <- Odds[Odds$P.Val <= 0.05 & Odds$intersect > 1,]
                         return(Odds)
                       })
  return(ModuleOdds)
}
