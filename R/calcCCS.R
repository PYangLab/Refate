#' @title Calculate Cellular Conversion Score (CCS)
#'
#' @description This function calculates CCS for genes based on fold change (fc) and precomputed CPS values.
#'
#' @param fc A named numeric vector of fold-change values for genes.
#' @param specie A string specifying the species ("human" or "mouse"). Default is "human".
#'
#' @return A named numeric vector of CCS values.
#'
calcCCS <- function(fc, specie = "human") {
  if (specie == "human") {
    o <- intersect(names(fc), names(CPSs$human.combined))
    CCS <- fc[o] * CPSs$human.combined[o]
  } else if (specie == "mouse") {
    o <- intersect(names(fc), names(CPSs$mouse.combined))
    CCS <- fc[o] * CPSs$mouse.combined[o]
  }
  return(CCS)
}



