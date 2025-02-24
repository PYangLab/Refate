#' @title Build Trans-regulatory Networks (TRNs)
#'
#' @description This function builds TRNs and partitions them into communities using Louvain clustering.
#'
#' @param CCS A named numeric vector of CCS values for genes.
#' @param TFtarget A data frame of TF-target interactions with columns for source, target, and interaction type.
#' @param consistent A boolean indicating whether to use consistent filtering for all genes. Default is FALSE.
#'
#' @return A list of community subgraphs (igraph objects).
#'
community_graphs_build <- function(CCS, TFtarget, consistent = FALSE) {
  library(igraph)

  if (consistent) {
    TFtarget_filtered <- TFtarget[TFtarget[,1] %in% allGenes | TFtarget[,2] %in% allGenes, ]
  } else {
    TFtarget_filtered <- TFtarget[TFtarget[,1] %in% names(CCS) | TFtarget[,2] %in% names(CCS), ]
  }

  g <- graph_from_data_frame(TFtarget_filtered, directed = FALSE)
  set.seed(123)
  community <- cluster_louvain(g)
  community_list <- split(V(g), membership(community))

  community_graphs <- lapply(unique(membership(community)), function(i) {
    induced_subgraph(g, which(membership(community) == i))
  })

  return(community_graphs)
}

