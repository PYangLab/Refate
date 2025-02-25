#' @title Get Drug Annotations from Multiple Databases
#' @description This function retrieves all available annotations for a given drug from multiple drug databases.
#' @param drug_name A string representing the name of the drug.
#' @return A list of data frames containing drug annotations from different databases.
#'
#'

get_drug_annotations <- function(drug_name) {
  annotations <- list()
  idx <- integer()

  for (i in seq_along(drugbases_metas)) {
    df <- drugbases_metas[[i]]
    if (drug_name %in% df$drug.re) {
      drug_info <- df[df$drug.re == drug_name, setdiff(colnames(df), "drug.re")]
      if (nrow(drug_info) > 0) {
        idx <- c(idx, i)
        annotations[[length(annotations) + 1]] <- drug_info
      }
    }
  }
  names(annotations) <- names(drugbases_metas)[idx]
  return(annotations)
}

#' @title Get Standardized Drug Name
#' @description This function retrieves the standardized name of a given drug from multiple databases.
#' @param drug A string representing the drug identifier.
#' @return A string representing the standardized drug name.
#'
get_drug_name <- function(drug) {
  annotation <- get_drug_annotations(drug)
  ns <- names(annotation)
  drug.name <- NA

  if (any(grepl("CTD", ns))) {
    drug.name <- annotation[["CTD.meta"]]$X..ChemicalName[1]
  } else if (any(grepl("STITCH", ns))) {
    drug.name <- annotation[["STITCH.meta"]]$drug[1]
  } else if (any(grepl("drugBank", ns))) {
    drug.name <- annotation[["drugBank.meta"]]$name[1]
  } else if (any(grepl("DrugRepurpose", ns))) {
    drug.name <- annotation[["DrugRepurpose.meta"]]$pert_iname[1]
  } else if (any(grepl("TTD", ns))) {
    drug.name <- annotation[["TTD.meta"]]$drug[1]
  } else if (any(grepl("DGIdb", ns))) {
    drug.name <- annotation[["DGIdb.meta"]]$drug_claim_name[1]
  }

  return(ifelse(is.na(drug.name), "Unknown Drug Name", drug.name))
}

#' @title Get Drug Description or Mechanism of Action (MOA)
#' @description This function retrieves the description or mechanism of action (MOA) of a given drug from multiple databases.
#' @param drug A string representing the drug identifier.
#' @return A string containing the drug description or MOA.
#'
get_drug_description <- function(drug) {
  annotation <- get_drug_annotations(drug)
  ns <- names(annotation)
  drug.des <- NA

  if (any(grepl("drugBank", ns))) {
    drug.des <- annotation[["drugBank.meta"]]$description[1]
  } else if (any(grepl("DrugRepurpose", ns))) {
    drug.des <- annotation[["DrugRepurpose.meta"]]$moa[1]
  } else if (any(grepl("TTD", ns))) {
    drug.des <- annotation[["TTD.meta"]]$MOA[1]
  }

  return(ifelse(is.na(drug.des), "No Description Available", drug.des))
}
