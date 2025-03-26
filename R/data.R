#' @title Cell Propensity Scores (CPS) for Human and Mouse Genes
#' @description Precomputed Cell Propensity Scores (CPS) for human and mouse genes, derived from multimodal single-cell atlas data.
#' @usage data(CPSs)
#' @format An object of class \code{data.frame} with 20,000 rows and 3 columns:
#' \describe{
#'   \item{\code{gene}}{Gene symbol}
#'   \item{\code{species}}{Species, either "human" or "mouse"}
#'   \item{\code{CPS}}{Cell Propensity Score}
#' }
#'
#' @details This dataset contains the precomputed CPS for human and mouse genes, which can be used directly in the analysis pipeline.
#'
#' @examples
#' data(CPSs)
#' head(CPSs)
#'
#' @name CPSs
#' @docType data
"CPSs"

#' @title Fold change examples
#' @description Example fold change values for genes.
#' @usage data(fc.fibroblast2hiPSC)
#' @format An object of class \code{vector}
#' \describe{
#'  \item{\code{name}}{Gene symbol}
#'  \item{\code{value}}{Fold change value}
#'  }
#'
#'  @examples
#'  data(fc.fibroblast2hiPSC)
#'  head(fc.fibroblast2hiPSC)
#'
#'  @name fc.fibroblast2hiPSC
#'  @docType data
"fc.fibroblast2hiPSC"

#' @title TF-target interactions
#' @description Curated TF-target interactions.
#' @usage data(TFtarget)
#' @format An object of class \code{data.frame} with 5 rows and 3 columns:
#' \describe{
#' \item{\code{gene1}}{TF gene symbol}
#' \item{\code{gene2}}{Target gene symbol}
#' }
#'
#' @examples
#' data(TFtarget)
#' head(TFtarget)
#'
#' @name TFtarget
#' @docType data
"TFtarget"

#' @title Drug-Target Interaction Data with Weights
#'
#' @description A list of drug-target interactions, where each element represents a drug and its target genes,
#' along with corresponding interaction weights (e.g., 1 for activation, -1 for repression).
#'
#' @format A list of named numeric vectors:
#' \describe{
#'   \item{Each element}{A named numeric vector representing target genes and their interaction weights for a specific drug.}
#'   \item{Gene names}{The names of genes targeted by the drug.}
#'   \item{Interaction weights}{Numeric values indicating the type of interaction: \code{1} for activation, \code{-1} for repression.}
#' }
#'
#' @examples
#' # Access the first drug in the list
#' drugDB.weight[[1]]
#'
#' # View target genes and interaction weights for the second drug
#' drugDB.weight[[2]]
#'
#' @name drugDB.weight
#' @docType data
"drugDB.weight"

#' @title Transcription factor (TF) and kinases
#' @description A list of transcription factors (TFs) and kinases.
#' @usage data(TFs_kinases)
#' @format A list of named character vectors
#' \describe{
#'  \item{Each element}{A named character vector representing TFs or kinases.}
#'  \item{Gene names}{The names of TFs or kinases.}
#'  }
#'
#' @examples
#' data(TFs_kinases)
#' head(TFs_kinases)
#'
#' @name TFs_kinases
#' @docType data
"TFs_kinases"

#' @title Protein-protein interaction (PPI) data
#' @description Protein-protein interaction data.
#' @usage data(StringDB)
#' @format An object of class \code{data.frame} with 45142 rows and 2 columns:
#' \describe{
#' \item{\code{gene1}}{Gene symbol}
#' \item{\code{gene2}}{Gene symbol}
#' }
#'
#' @examples
#' data(StringDB)
#' head(StringDB)
#'
#' @name StringDB
#' @docType data
#'
"StringDB"

#' @title meta data of drug databases
#' @description Meta data of drug databases.
#' @usage data(drugbases_metas)
#'
#' @name drugbases_metas
#' @docType data
#'
"drugbases_metas"

#' @title TTD Drug-target interaction data
#' @usage data(TTD)
#'
#' @name TTD
#' @docType data
"TTD"

#' @title CTD Drug-target interaction data
#' @usage data(CTD)
#'
#' @name CTD
#' @docType data
"CTD"

#' @title drugBank Drug-target interaction data
#' @usage data(drugBank)
#'
#' @name drugBank
#' @docType data
"drugBank"

#' @title STITCH Drug-target interaction data
#' @usage data(STITCH)
#'
#' @name STITCH
#' @docType data
"STITCH"

#' @title DrugRepurpose Drug-target interaction data
#' @usage data(DrugRepurpose)
#'
#' @name DrugRepurpose
#' @docType data
"DrugRepurpose"

#' @title DGIdb Drug-target interaction data
#' @usage data(DGIdb)
#'
#' @name DGIdb
#' @docType data
"DGIdb"
