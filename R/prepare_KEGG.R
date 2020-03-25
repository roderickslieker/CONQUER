#' This function was obtained from clusterProfiler and was not exported. Therefore this function has been included here
#' @param species [[character]]
#' @param KEGG_Type [[character]]
#' @param keyType [[character]]
#' @import clusterProfiler
prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
    build_Anno(kegg$KEGGPATHID2EXTID,
               kegg$KEGGPATHID2NAME)
}
