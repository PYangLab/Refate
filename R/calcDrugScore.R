#' @title Build Backbone Network Using PPI and TF-Target Interactions
#'
#' @description This function builds a backbone network using PPI and TF-target interactions.
#' The network is partitioned into community subgraphs using the Louvain clustering algorithm.
#'
#' @param CCS A named numeric vector of Cellular Conversion Scores (CCS) for genes. Used to filter relevant genes in the network.
#' @param TFtarget A data frame of transcription factor (TF) target interactions with columns 'source' (TF), 'target' (gene), and 'type' (interaction type).
#' @param consistent A boolean indicating whether to use consistent filtering for all genes. Default is \code{FALSE}.
#'
#' @return A list of community subgraphs represented as igraph objects. Each subgraph corresponds to a densely connected community in the backbone network.
#'
#' @examples
#' # Build backbone network and partition into community subgraphs
#' community_graphs <- community_graphs_build(CCS, TFtarget)
#'
#' # View the number of community subgraphs
#' length(community_graphs)
#'
community_graphs_build <- function(CCS, TFtarget) {
  library(igraph)

  TFtarget_filtered <- TFtarget[TFtarget[,1] %in% names(CCS) | TFtarget[,2] %in% names(CCS), ]

  g <- graph_from_data_frame(TFtarget_filtered, directed = FALSE)
  set.seed(123)
  community <- cluster_louvain(g)
  community_list <- split(V(g), membership(community))

  community_graphs <- lapply(unique(membership(community)), function(i) {
    induced_subgraph(g, which(membership(community) == i))
  })

  return(community_graphs)
}


#' @title Predict Drug Impact on Trans-regulatory Networks (TRNs)
#'
#' @description This function predicts the impact of a given drug on TRNs by identifying the targeted subgraphs
#' and calculating cumulative scores based on the interaction of drug targets with genes in the community subgraphs.
#'
#' @param drug A string representing the name of the drug to be scored.
#' @param community_graphs A list of community subgraphs (igraph objects) representing partitioned TRNs.
#' @param TFWeight A numeric value for weighting transcription factor (TF) influence in the scoring calculation.
#' @param Scoring A string specifying the scoring method. Options include "scoring1", "scoring2", "scoring3", etc.
#' @param CCS A named numeric vector representing Cellular Conversion Scores (CCS) for genes.
#' @param drugdb A list containing drug-target relationships, where each drug is associated with a named vector
#' of target genes and their corresponding weights.
#'
#' @return A named numeric vector representing the cumulative drug score for the given drug.
#'
#' @examples
#' community_graphs <- community_graphs_build(CCS, TFtarget)
#' drugdb <- list(Drug1 = c(GeneA = 0.8, GeneB = 0.5), Drug2 = c(GeneC = 0.9))
#' CCS <- c(GeneA = 1.2, GeneB = 0.8, GeneC = 1.5)
#' drugPrediction("Drug1", community_graphs, TFWeight = 1, Scoring = "scoring1", CCS = CCS, drugdb = drugdb)
#'

drugPrediction <- function(drug,
                           community_graphs,
                           TFWeight,
                           CCS,
                           drugdb) {

  drug_target_genes <- names(drugdb[[drug]])
  # Initialize an empty list to store the results
  targeted_subgraphs <- list()
  targeted_subgraph_indices <- c()
  # Loop over the subgraphs
  for (i in 1:length(community_graphs)) {
    # Get the gene names for this subgraph
    subgraph_genes <- V(community_graphs[[i]])$name
    # Find which of these genes are targeted by the drug
    targeted_genes <- intersect(subgraph_genes, drug_target_genes)
    # If any genes are targeted, store the result
    if (length(targeted_genes) > 0) {
      TFweight = ifelse(subgraph_genes %in% TFs, as.numeric(TFWeight), 1) # Weighting for TFs
      CCS.sum <- sum(TFweight * CCS[subgraph_genes], na.rm = T)
      drug.effect <- drugdb[[drug]][targeted_genes]
      weighted.effect <- sum(drugdb[[drug]][targeted_genes], na.rm = T)
      weighted.effect.nodir <- sum(abs(drugdb[[drug]][targeted_genes]), na.rm = T)

      targeted_subgraphs[[length(targeted_subgraphs) + 1]] <- list(
        "Subgraph" = i,
        "TargetedGenes" = targeted_genes,
        "NumTargetedGenes" = length(targeted_genes),
        "TotalGenes" = length(subgraph_genes),
        "ProportionTargeted" = length(targeted_genes) / length(subgraph_genes),
        "CCS.sum" = CCS.sum,
        "drug.effect" = drug.effect,
        "weighted.effect" = weighted.effect,
        "weighted.effect.nodir" = weighted.effect.nodir
      )
      targeted_subgraph_indices <- c(targeted_subgraph_indices, i)
    }
  }

  scoring <- c()
  if (length(targeted_subgraphs) > 0) {
    idx <- sapply(targeted_subgraphs, function(x)x$CCS.sum > 0)

    scores <- sapply(targeted_subgraphs, function(x) {
      tmp <- sapply(x$drug.effect, function(y) y * x$CCS.sum)
      sum(tmp[which(tmp > 0)])})
    names(scores) <- targeted_subgraph_indices
    scoring[drug] <- sum(unlist(scores), na.rm = T)
  }
  return(scoring)
}


#' @title Prioritize Drugs Based on TRNs
#'
#' @description This function scores and ranks drugs based on their impact on TRNs using CCS and community graphs.
#'
#' @param CCS A named numeric vector of CCS values for genes.
#' @param community_graphs A list of community subgraphs (igraph objects).
#' @param community_cutoff An integer specifying the number of genes in a community subgraph.
#' @param drugdb A list of drug-target relationships with gene weights.
#' @param TFWeight A numeric value for weighting transcription factors.
#' @param Scoring A string specifying the scoring method. Options: "scoring1", "scoring2", "scoring3", etc.
#'
#' @return A named numeric vector of drug scores.
#'
calcDrugScore <- function (CCS, 
                           drugDB = drugDB.weight,
                           community_cutoff = 10, 
                           community_cutoff2 = 600, 
                           filtering = TRUE, 
                           drugTarget_cutoff = 500, 
                           TFWeight = 2,
                           direction = TRUE)
{
  if (direction) {
    drugDB <- drugDB
  } else {
    drugDB <- lapply(drugDB, function(x) {
      x <- abs(x)
      return(x)
    })
  }
  
  drugdb <- drugDB[sapply(drugDB, length) < drugTarget_cutoff]
  allDrugs <- names(drugdb)
  community_graphs <- community_graphs_build(CCS, TFtarget)
  community_graphs <- community_graphs[sapply(community_graphs, length) > community_cutoff]
  community_graphs <- community_graphs[sapply(community_graphs, length) < community_cutoff2]
  print(paste("The number of sub-graphs is:", length(community_graphs)))
  CC.prediction <- mclapply(allDrugs, function(drug) {
    drugPrediction(drug, community_graphs, TFWeight, CCS, 
                   drugdb)
  }, mc.cores = parallel::detectCores() - 2)
  CC.prediction <- do.call(c, CC.prediction)
  if (filtering) {
    CC.prediction <- CC.prediction[which(CC.prediction != 0)]
  }
  scaled.rank <- rank(CC.prediction)/(length(CC.prediction) + 1)
  scaled.rank <- sort(scaled.rank, decreasing = T)
  return(scaled.rank)
}


