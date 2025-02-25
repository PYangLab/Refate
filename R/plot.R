#' @title Plot Gene Ranks with Highlights for TFs and Specific Genes
#'
#' @description This function generates a bar plot to visualize the top-ranked genes based on their ranking scores, with transcription factors (TFs) and specific highlighted genes (`gt`) shown in different colors.
#'
#' @param gene.rank A named numeric vector representing the ranking scores for genes.
#' @param gt A character vector of gene names to be highlighted in the plot.
#'
#' @return A ggplot object displaying the top 20 ranked genes that are either transcription factors (TFs) or in the specified `gt` set.
#'
#' @examples
#' # Example data
#' gene.rank <- setNames(runif(100), paste0("Gene", 1:100))
#' TFs <- c("Gene5", "Gene10", "Gene20")  # Example TFs
#' gt <- c("Gene15", "Gene50")  # Genes to highlight
#'
#' # Generate the plot
#' rank.plot(gene.rank, gt)
#'
rank.plot <- function(gene.rank, gt) {
  library(ggplot2)

  dftoplot <- data.frame(gene.rank)
  dftoplot$type <- "Others"
  dftoplot$type[which(names(gene.rank) %in% TFs)] <- "TF"
  dftoplot$type[which(names(gene.rank) %in% gt)] <- "Highlight"
  dftoplot$gene <- rownames(dftoplot)
  dftoplot <- dftoplot[order(dftoplot$gene.rank, decreasing = TRUE), ]

  dftoplot %>%
    dplyr::filter(type %in% c("TF", "Highlight")) %>%
    head(20) %>%
    ggplot(aes(y = gene.rank, x = reorder(gene, gene.rank), fill = type)) +
    scale_fill_manual(values = c("#E22926", "#595A5A")) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    xlab("TFs") +
    ylab("CCS")
}

