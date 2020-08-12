getGenes <- function(rsID, SNPsin=loadedSNPs)
{
  temp <- SNPsin[[rsID]]
  temp <- c(temp$eQTLs$gene, temp$eQTLsTrans$gene) %>% unique()
  temp <- sort(temp)
  return(temp)
}
