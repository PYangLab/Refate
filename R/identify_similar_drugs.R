
drugScoringPerSubgraph <- function(drug,
                                   community_graphs,
                                   TFWeight,
                                   CCS,
                                   drugdb,
                                   Scoring) {

  drug_target_genes <- names(drugdb[[drug]])
  # Initialize an empty list to store the results
  targeted_subgraphs <- list()
  targeted_subgraph_indices <- c()
  # Loop over the subgraphs
  for (i in 1:length(community_graphs)) {
    print(i)
    # Get the gene names for this subgraph
    subgraph_genes <- V(community_graphs[[i]])$name
    # Find which of these genes are targeted by the drug
    targeted_genes <- intersect(subgraph_genes, drug_target_genes)
    # If any genes are targeted, store the result
    if (length(targeted_genes) > 0) {
      TFweight = ifelse(subgraph_genes %in% TFs, as.numeric(TFWeight), 1)
      CCS.sum <- sum(TFweight * CCS[subgraph_genes], na.rm = T)
      weighted.effect <- sum(drugdb[[drug]][targeted_genes], na.rm = T)
      drug.effect <- drugdb[[drug]][targeted_genes]

      targeted_subgraphs[[length(targeted_subgraphs) + 1]] <- list(
        "Subgraph" = i,
        "TargetedGenes" = targeted_genes,
        "NumTargetedGenes" = length(targeted_genes),
        "TotalGenes" = length(subgraph_genes),
        "ProportionTargeted" = length(targeted_genes) / length(subgraph_genes),
        "CCS.sum" = CCS.sum,
        "drug.effect" = drug.effect,
        "weighted.effect" = weighted.effect
      )
    } else {
      targeted_subgraphs[[length(targeted_subgraphs) + 1]] <- list(
        "Subgraph" = i,
        "TargetedGenes" = targeted_genes,
        "NumTargetedGenes" = length(targeted_genes),
        "TotalGenes" = length(subgraph_genes),
        "ProportionTargeted" = length(targeted_genes) / length(subgraph_genes),
        "CCS.sum" = 0,
        "drug.effect" = 0,
        "weighted.effect" = 0
      )
    }
  }

  scores <- sapply(targeted_subgraphs, function(x) {
    tmp <- sapply(x$drug.effect, function(y) y * x$CCS.sum)
    sum(tmp[which(tmp > 0)])})
  names(scores) <- targeted_subgraph_indices

  return(scores)
}


#' @title Identify Drugs Similar to Known Compounds
#'
#' @description This function identifies drugs that are most similar to known chemical compounds based on how they target subgraphs.
#' Similarity is measured by computing the correlation between subgraph-targeting scores.
#'
#' @param top_drugs A vector of drug names for which subgraph-level scores are calculated.
#' @param community_graphs A list of community subgraphs (igraph objects).
#' @param CCS A named numeric vector representing Cellular Conversion Scores (CCS) for genes.
#' @param drugDB.weight A list of drug-target relationships with corresponding weights.
#' @param known_drug A string representing the name of the known chemical compound.
#' @param method A string specifying the correlation method. Options are \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}. Default is \code{"pearson"}.
#' @param top Integer specifying how many top similar drugs to return.
#'
#' @return A named numeric vector of correlation values for the top similar drugs.
#'
identify_similar_drugs <- function(top_drugs,
                                   CCS,
                                   drugDB.weight,
                                   known_drug,
                                   method = "pearson",
                                   top = 10) {
  library(parallel)

  community_graphs <- community_graphs_build(CCS, TFtarget)
  community_graphs <- community_graphs[sapply(community_graphs, length) > 10]
  community_graphs <- community_graphs[sapply(community_graphs, length) < 600]

  # Calculate subgraph-level scores for all drugs in parallel
  drugByconversion_subgraph <- mclapply(top_drugs, function(drug) {
    scores <- drugScoringPerSubgraph(
      drug = drug,
      community_graphs = community_graphs,
      TFWeight = 2,
      drugdb = drugDB.weight,
      CCS = CCS
    )
    names(scores) <- 1:length(community_graphs)
    return(unlist(scores))
  }, mc.cores = parallel::detectCores() - 1)

  # Convert to a data frame with each drug as a column
  names(drugByconversion_subgraph) <- top_drugs
  drugByconversion_subgraph.re <- do.call(cbind, drugByconversion_subgraph)

  # Calculate correlation of each drug with the known drug
  correlations <- sapply(colnames(drugByconversion_subgraph.re), function(drug) {
    cor(drugByconversion_subgraph.re[, drug], drugByconversion_subgraph.re[, known_drug], method = method)
  })

  # Sort by similarity and return the top results
  sorted_correlations <- sort(correlations, decreasing = TRUE)
  top.similarity.drugs <- sorted_correlations[2:(top + 1)]  # Exclude the known drug itself

  return(top.similarity.drugs)
}
