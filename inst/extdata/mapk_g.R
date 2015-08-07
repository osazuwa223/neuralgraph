library(KEGGgraph)
library(igraph)
g_nell <- tempfile() %T>%
  {KEGGgraph::retrieveKGML("04010", organism="hsa", destfile=., method="curl", quiet=TRUE)} %>%
  {KEGGgraph::parseKGML2Graph(., expandGenes=FALSE)} 
vertex_list <- KEGGgraph::getKEGGnodeData(g_nell) %>%
  {data.frame(
    kegg = unlist(lapply(., function(item) item@name[1])),
    label = unlist(lapply(., function(item)
      strsplit(item@graphics@name, ",")[[1]][1])), stringsAsFactors = F)}
g_init <- igraph.from.graphNEL(g_nell) 
V(g_init)$name <- vertex_list$kegg 
vertex_list <- dplyr::filter(vertex_list, !duplicated(kegg))
edge_list <- KEGGgraph::getKEGGedgeData(g_nell) %>%
  lapply(function(item){
    if(length(item@subtype) > 0){
      subtype_info <- item@subtype
      # KEGG uses a hierarchy of term for describing terms
      # for example, the first edge type is "activation", the second is "phosphorylation"
      # where phosphorylation is a type of activation.  The second term is more specific than
      # the first, so when it is provided, use it in lieu of the first type.
      if(length(subtype_info) > 1) {
        return(subtype_info[[2]]@name)
      } else {
        return(subtype_info$subtype@name)
      }
    } 
    NA
  }) %>%
  unlist %>%
  {cbind(get.edgelist(g_init), type = .)} %>%
    data.frame %>%
  {dplyr::filter(.,type == "phosphorylation")}
edge_list <- edge_list %>%
  as.data.frame %>%
  unique
vertex_list <- vertex_list %>%
  unique %>%
  {dplyr::filter(., !duplicated(kegg))}
g <- graph.data.frame(edge_list, directed = TRUE, vertices = vertex_list)
V(g)$kid <- V(g)$name
V(g)$name <- V(g)$label
g <- g - V(g)[igraph::degree(g) == 0]
mapk_g <- g
devtools::use_data(mapk_g, overwrite = TRUE)
