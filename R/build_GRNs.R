#' @title Build GRNs for Candidate Drug
#'
#' @description This function builds a trans-regulatory network (GRN) for a candidate drug by connecting the drug to kinases, transcription factors (TFs), and target genes through PPI and TF-target interactions.
#'
#' @param drug A string representing the name of the drug.
#' @param fc A named numeric vector of fold-change values for genes.
#' @param kinome A vector of gene names representing kinases.
#' @param TFs A vector of gene names representing transcription factors.
#' @param String A data frame containing PPI interactions with columns "gene1" and "gene2".
#' @param TFtarget A data frame containing TF-target interactions with columns "gene1" (TF) and "gene2" (target gene).
#' @param CCS_threshold A numeric value specifying the CCS threshold for filtering target genes. Default is 0.02.
#'
#' @return A force-directed network visualization (HTML widget) and an igraph object representing the network.
#'
#' @examples
#' build_GRNs("quercetin", fc, kinome, TFs, String, TFtarget)
#'
build_GRNs <- function(drug, CCS, kinome, TFs, StringDB, TFtarget, CCS_threshold = 0.02) {
  library(igraph)
  library(networkD3)
  library(dplyr)

  targets <- na.omit(CCS[names(drugDB.weight[[drug]])] * drugDB.weight[[drug]])

  # Identify kinases that are drug targets
  kinases <- targets[which(names(targets) %in% kinome)]
  links1 <- data.frame(from = drug, to = names(kinases))
  links1 <- links1 %>% filter(to %in% names(CCS)[CCS > 0])

  # Kinase to TF via PPI
  links2 <- StringDB %>%
    filter(gene1 %in% links1$to & gene2 %in% TFs) %>%
    select(gene1, gene2)
  colnames(links2) <- c("from", "to")
  links1 <- links1 %>% filter(to %in% links2$from)
  links2 <- links2 %>% filter(to %in% names(CCS)[CCS > 0])

  # TF to target via TF-target interactions
  links3 <- TFtarget %>%
    filter(gene1 %in% links2$to & gene2 %in% names(CCS)) %>%
    select(gene1, gene2)
  colnames(links3) <- c("from", "to")
  links3 <- links3 %>% filter(to %in% names(CCS)[CCS > CCS_threshold])

  # Combine links
  links <- rbind(links1, links2, links3)

  # Assign node types
  node_types <- rep(NA, length(unique(c(links$from, links$to))))
  names(node_types) <- unique(c(links$from, links$to))
  node_types[links1$from] <- "drug"
  node_types[links1$to] <- "kinase"
  node_types[links2$to] <- "TF"
  node_types[links3$to] <- "target"

  # Convert to igraph object
  g <- graph_from_data_frame(links, directed = FALSE)

  # Group nodes for visualization
  groups <- ifelse(V(g)$name %in% names(node_types), node_types[V(g)$name], "unknown")

  # Convert to networkD3 format
  graph_d3 <- igraph_to_networkD3(g, group = groups)

  # Define color scale for node types
  ColourScale <- 'd3.scaleOrdinal()
                  .domain(["drug", "kinase", "TF", "target"])
                  .range(["black", "orange", "blue", "grey"]);'

  # Visualize the network
  p <- forceNetwork(
    Links = graph_d3$links, Nodes = graph_d3$nodes, NodeID = 'name',
    Source = 'source', Target = 'target', Group = 'group', opacityNoHover = 1,
    opacity = 1, zoom = TRUE, fontSize = 12, colourScale = ColourScale
  )
  return(p)
}
