#' @importFrom clusterProfiler download_KEGG
prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
    build_Anno(kegg$KEGGPATHID2EXTID,
               kegg$KEGGPATHID2NAME)
}
