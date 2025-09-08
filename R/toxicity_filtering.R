
#' @title Toxicity Filtering via PubChem and ChEBI
#' @description Filter drug candidates based on toxicity signals from PubChem PUG-View and ChEBI roles.
#' Resolves names -> CIDs (via your all_drug_cids and optional PubChem fallback), downloads/reads JSON,
#' extracts safety/role features, and applies conservative filters.
#'
#' @param drug_names Character vector of drug names.
#' @param all_drug_cids Data frame with columns `cid` and `drug_name` (name->CID map).
#' @param json_dir Directory to store/download PubChem JSON (default: user cache).
#' @param download_missing Logical, download missing PubChem JSON files.
#' @param id_col Preferred name column for output labels:
#'   one of "InputName","DrugName","Display","CID".
#' @param use_pubchem_fallback Logical, try PubChem name->CID for unmatched names.
#'
#' @return A list with:
#'   - `resolution`: name->CID resolution table
#'   - `full_table`: scored table for all CIDs
#'   - `filtered_table`: table with Keep flag and preferred names
#'   - `filtered_names`: character vector of kept names (per `id_col`)
#'
#' @examples
#' \dontrun{
#' drug_names <- c("Aspirin", "Paracetamol", "UnknownDrug")
#' result <- filter_toxicity(drug_names, all_drug_cids)
#' result$filtered_names
#' }
#'
suppressPackageStartupMessages({
  library(jsonlite)
  library(httr)
  library(xml2)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tools)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --------------------------- config (defaults) --------------------------------
DEFAULT_H_EXCLUDE <- c(
  "H300","H301","H310","H311","H330","H331","H340","H341","H350","H351","H360","H400"
)

DEFAULT_ICON_EXCLUDE <- c(
  "Acute Toxic","Explosive","Compressed Gas","Flammable","Oxidizer"
)

DEFAULT_EXCLUDE_ROLES <- setdiff(c(
  # Toxicity / harmful classes
  "mutagen", "neurotoxin", "explosive", "hepatotoxic agent", "carcinogenic agent",
  "teratogenic agent", "genotoxin", "poison", "persistent organic pollutant", "environmental contaminant",

  # Pesticides / agrochemicals
  "insecticide", "pyrethroid ester insecticide", "pyrethroid ester acaricide", "herbicide",
  "agrochemical", "pesticide synergist", "antifungal agrochemical", "antibiotic insecticide", "neonicotinoid insectide",

  # Irrelevant / unspecific
  "fragrance", "flavouring agent", "non-polar solvent", "solvent", "polar aprotic solvent", "fuel additive", "volatile oil component", "protic solvent",
  "food acidity regulator", "sweetening agent", "food additive", "NMR chemical shift reference compound",
  "refrigerant", "bleaching agent",

  # Non-human/environmental metabolites
  "algal metabolite", "fungal metabolite", "Penicillium metabolite", "Aspergillus metabolite", "Escherichia coli metabolite",
  "Saccharomyces cerevisiae metabolite", "Daphnia magna metabolite", "Daphnia galeata metabolite", "Mycoplasma genitalium metabolite",
  "bacterial metabolite", "plant metabolite", "marine metabolite",

  # Detrimental mechanisms
  "alkylating agent", "hapten", "oneirogen", "sensitiser", "emetic", "glutathione depleting agent",
  "electrophilic reagent",

  # Agents targeting unrelated processes
  "antinematodal drug", "antirheumatic drug", "schistosomicide drug", "antimalarial",
  "antibacterial drug", "antitubercular agent", "antiamoebic agent", "antifungal agent",
  "antileishmanial agent", "antiprotozoal drug", "antiviral agent", "anti-HIV agent",
  "anti-HBV agent", "anticoronaviral agent", "leprostatic drug", "analgesic", "anaesthetic",
  "local anaesthetic", "inhalation anaesthetic", "non-narcotic analgesic", "opioid analgesic",
  "antipyretic", "antiepileptic", "sedative", "anxiolytic drug", "antipsychotic",
  "juvenile hormone mimic", "animal growth promotant",
  "oxytocic", "female contraceptive drug", "synthetic oral contraceptive", "anticholesteremic drug",
  "cardiovascular drug", "antiarrhythmia drug", "cardiotonic drug", "bronchodilator agent",
  "anti-asthmatic drug", "gastrointestinal drug", "laxative", "antiemetic", "antipruritic drug",
  "anti-allergic agent", "H1-receptor antagonist", "hypoglycemic agent", "thyroid hormone",
  "antithyroid drug", "antiglaucoma drug", "mydriatic agent", "antiparkinson drug", "antidyskinesia agent"
), NA)

# name-keyed manual role patches (applied AFTER name join)
MANUAL_ROLE_PATCH <- list(
  "s12dichlorovinylcysteine"="nephrotoxin",
  "bisphenolb"="endocrine disruptor",
  "23dimethoxy14naphthoquinone"=c("cytotoxin","ROS generators"),
  "nsc228155"=c("cytotoxicity","cytotoxic"),
  "abt348"="cytotoxicity",
  "insm18"="Hepatotoxicity",
  "3nitrobenzanthrone"=c("mutagen","ROS generators"),
  "mg132"="proteasome inhibitor",
  "naphthoquinones"="ROS generators",
  "idarubicin"="cytotoxic",
  "fludarabine"="cytotoxic",
  "acidblue129"="dye",
  "oxazolone"="toxic",
  "emilace"="uncharacterized",
  "didecyldimethylammonium"="biocide",
  "pmid19788238c66"="uncharacterized"
)

# --------------------------- utilities ----------------------------------------
download_pubchem_json <- function(cids, json_dir, polite_delay = 0.25) {
  dir.create(json_dir, showWarnings = FALSE, recursive = TRUE)
  ua <- httr::user_agent("yourpkg/0.1 (contact: you@example.com)")
  for (cid in cids) {
    fp <- file.path(json_dir, sprintf("CID_%s.json", cid))
    if (file.exists(fp)) next
    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/", cid)

    resp <- try(
      httr::RETRY("GET", url, ua, httr::timeout(20),
                  terminate_on = c(400L, 401L, 403L, 404L),
                  times = 4, pause_base = 0.5, pause_cap = 4),
      silent = TRUE
    )

    code <- suppressWarnings(try(httr::status_code(resp), silent = TRUE))
    if (inherits(resp, "try-error") || is.na(code) || code != 200) {
      warning("Failed to retrieve CID: ", cid, " (HTTP ", ifelse(is.na(code), "NA", code), ")")
      next
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    tf <- tempfile("pubchem_", tmpdir = json_dir, fileext = ".json.tmp")
    writeLines(txt, tf, useBytes = TRUE)
    file.rename(tf, fp)
    Sys.sleep(polite_delay)
  }
  invisible(TRUE)
}

.find_section <- function(sections, heading) {
  for (section in sections) {
    if (!is.null(section$TOCHeading) && identical(section$TOCHeading, heading)) return(section)
    if (!is.null(section$Section)) {
      res <- .find_section(section$Section, heading)
      if (!is.null(res)) return(res)
    }
  }
  NULL
}

extract_pubchem_info <- function(json_path) {
  json_data <- try(jsonlite::fromJSON(json_path, simplifyVector = FALSE), silent = TRUE)
  if (inherits(json_data, "try-error") || is.null(json_data$Record$Section)) return(NULL)
  find_section <- function(heading) .find_section(json_data$Record$Section, heading)

  extract_computed_properties <- function() {
    computed <- find_section("Computed Properties"); if (is.null(computed)) return(NULL)
    props <- list()
    for (sub in computed$Section) {
      name <- sub$TOCHeading
      val <- tryCatch({
        if (!is.null(sub$Information[[1]]$Value$Number)) {
          sub$Information[[1]]$Value$Number[[1]]
        } else if (!is.null(sub$Information[[1]]$Value$StringWithMarkup)) {
          sub$Information[[1]]$Value$StringWithMarkup[[1]]$String
        } else NA
      }, error = function(e) NA)
      props[[name]] <- val
    }
    props
  }

  extract_ghs_h <- function() {
    ghs <- find_section("GHS Classification"); if (is.null(ghs)) return(NA_character_)
    tryCatch({
      info <- unlist(lapply(ghs$Information, function(x) lapply(x$Value$StringWithMarkup, \(y) y$String)))
      codes <- unique(unlist(regmatches(info, gregexpr("H[0-9]{3}", info))))
      if (length(codes)==0) NA_character_ else codes
    }, error = function(e) NA_character_)
  }

  extract_ghs_icons <- function() {
    ghs <- find_section("Chemical Safety"); if (is.null(ghs)) return(NA_character_)
    tryCatch({
      icons <- lapply(ghs$Information[[1]]$Value$StringWithMarkup[[1]]$Markup, \(x) x$Extra)
      if (length(icons)==0) NA_character_ else unlist(icons)
    }, error = function(e) NA_character_)
  }

  extract_other_hazards <- function() {
    hz <- find_section("Hazard Classes and Categories"); if (is.null(hz)) return(NA_character_)
    tryCatch({
      vals <- unlist(lapply(hz$Information, \(x) lapply(x$Value$StringWithMarkup, \(y) y$String)))
      if (length(vals)==0) NA_character_ else vals
    }, error = function(e) NA_character_)
  }

  extract_mesh <- function() {
    s <- find_section("MeSH Pharmacological Classification"); if (is.null(s)) return(NA_character_)
    tryCatch({
      v <- sapply(s$Information, \(x) x$Value$StringWithMarkup[[1]]$String)
      if (length(v)==0) NA_character_ else v
    }, error = function(e) NA_character_)
  }

  extract_cas <- function() {
    # Note: CAS sometimes lives under "Other Identifiers" → "CAS"; this will catch top-level only.
    s <- find_section("CAS"); if (is.null(s)) return(NA_character_)
    tryCatch(s$Information[[1]]$Value$StringWithMarkup[[1]]$String, error = function(e) NA_character_)
  }

  extract_bioassays <- function() {
    s <- find_section("BioAssay Results"); if (is.null(s)) return(NA_character_)
    tryCatch({
      v <- sapply(s$Information, \(x) x$Value$ExternalTableName)
      if (length(v)==0) NA_character_ else v
    }, error = function(e) NA_character_)
  }

  extract_chebi_refs <- function() {
    s <- find_section("ChEBI Ontology"); if (is.null(s)) return(NA_character_)
    tryCatch({
      v <- sapply(s$Information, \(x) x$ReferenceNumber)
      if (length(v)==0) NA_character_ else v
    }, error = function(e) NA_character_)
  }

  extract_chem_classes <- function() {
    s <- find_section("Drugs"); if (is.null(s)) return(NA_character_)
    tryCatch({
      v <- sapply(s$Information, \(x) x$Value$StringWithMarkup[[1]]$String)
      if (length(v)==0) NA_character_ else v
    }, error = function(e) NA_character_)
  }

  extract_other_ids <- function() {
    oi <- find_section("Other Identifiers")
    sections <- if (is.null(oi)) list() else (oi$Section %||% list())
    other  <- .find_section(sections, "ChEBI ID")
    dsstox <- .find_section(sections, "DSSTox Substance ID")
    chebi_id  <- tryCatch(other$Information[[1]]$Value$StringWithMarkup[[1]]$String, error=function(e) NA_character_)
    dsstox_id <- tryCatch(dsstox$Information[[1]]$Value$StringWithMarkup[[1]]$String, error=function(e) NA_character_)
    list(CHEBI_ID = chebi_id, DSSTox_ID = dsstox_id)
  }

  other_ids <- extract_other_ids()

  list(
    ComputedProperties = extract_computed_properties(),
    GHS_Hazards        = extract_ghs_h(),
    GHS_Hazard_Icons   = extract_ghs_icons(),
    Other_Hazards      = extract_other_hazards(),
    MeSH_Classification= extract_mesh(),
    CAS_Number         = extract_cas(),
    CHEBI_ID           = other_ids$CHEBI_ID,
    DSSTox_ID          = other_ids$DSSTox_ID,
    BioAssays          = extract_bioassays(),
    ChEBI_Ontology     = extract_chebi_refs(),
    Chemical_Classes   = extract_chem_classes()
  )
}

get_chebi_roles <- function(chebi_id, cache_dir = NULL) {
  if (is.null(chebi_id) || length(chebi_id)==0 || is.na(chebi_id) || chebi_id=="") return(NA_character_)
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  fp <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(chebi_id, ".rds")) else NULL
  if (!is.null(fp) && file.exists(fp)) {
    roles <- readRDS(fp); return(if (length(roles)==0) NA_character_ else roles)
  }

  url <- paste0("https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId=", chebi_id)
  response <- try(httr::GET(url), silent = TRUE)
  if (inherits(response,"try-error") || httr::status_code(response) != 200) return(NA_character_)

  xml_txt <- httr::content(response, as = "text", encoding = "UTF-8")
  xml_doc <- try(xml2::read_xml(xml_txt), silent = TRUE)
  if (inherits(xml_doc, "try-error")) return(NA_character_)

  ns <- c(S="http://schemas.xmlsoap.org/soap/envelope/",
          chebi="https://www.ebi.ac.uk/webservices/chebi")
  roles <- xml2::xml_find_all(xml_doc, ".//chebi:OntologyParents[chebi:type='has role']/chebi:chebiName", ns) |>
    xml2::xml_text() |> unique()

  if (!is.null(fp)) saveRDS(roles, fp)
  if (length(roles)==0) return(NA_character_)
  roles
}

# ----------------------------- core builder -----------------------------------
build_pubchem_chebi_table <- function(cids,
                                      json_dir,
                                      download_missing = TRUE,
                                      chebi_cache_dir = file.path(json_dir, "_chebi_cache"),
                                      manual_role_patch = MANUAL_ROLE_PATCH,
                                      .progress = TRUE) {
  if (download_missing) download_pubchem_json(cids, json_dir)

  files <- file.path(json_dir, paste0("CID_", cids, ".json"))
  names(files) <- as.character(cids)

  if (.progress) message("Parsing PubChem JSON ...")
  parsed <- purrr::map(files, extract_pubchem_info)

  if (.progress) message("Fetching ChEBI roles ...")
  chebi_ids <- purrr::map_chr(parsed, \(x) if (is.null(x)) NA_character_ else x$CHEBI_ID %||% NA_character_)
  roles <- purrr::map(chebi_ids, get_chebi_roles, cache_dir = chebi_cache_dir)
  names(roles) <- names(parsed)  # ensure roles are named by CID

  # NOTE: MANUAL_ROLE_PATCH is name-keyed; applied later after join (see filter_toxicity()).
  list(pubchem = parsed, chebi_roles = roles)
}

# ----------------------------- scoring & flags --------------------------------
score_toxicity_table <- function(pubchem_list,
                                 chebi_roles,
                                 h_exclude = DEFAULT_H_EXCLUDE,
                                 icon_exclude = DEFAULT_ICON_EXCLUDE,
                                 exclude_roles = DEFAULT_EXCLUDE_ROLES) {

  safe_get <- function(x, what) if (is.null(x[[what]])) NA else x[[what]]
  get_prop <- function(x, key) {
    val <- tryCatch(safe_get(x$ComputedProperties, key), error=function(e) NA)
    if (length(val)>1) val <- val[[1]]
    suppressWarnings(as.numeric(val))
  }

  df <- tibble::tibble(
    CID = names(pubchem_list),
    MW  = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Molecular Weight")),
    TPSA = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Topological Polar Surface Area")),
    RotatableBonds = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Rotatable Bond Count")),
    XLogP3 = purrr::map_chr(pubchem_list, ~ {
      x <- safe_get(.x$ComputedProperties, "XLogP3"); if (is.null(x)) NA_character_ else as.character(x)
    }),
    HB_Donor = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Hydrogen Bond Donor Count")),
    HB_Acceptor = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Hydrogen Bond Acceptor Count")),
    Formal_charge = purrr::map_dbl(pubchem_list, ~ get_prop(.x, "Formal Charge")),
    CAS_Number = purrr::map_chr(pubchem_list, ~ safe_get(.x, "CAS_Number") %||% NA_character_),
    has_BioAssay = purrr::map_lgl(pubchem_list, ~ !all(is.na(safe_get(.x, "BioAssays")))),
    has_CHEBI_ID = purrr::map_lgl(pubchem_list, ~ !is.na(safe_get(.x, "CHEBI_ID"))),
    has_DSSTox = purrr::map_lgl(pubchem_list, ~ !is.na(safe_get(.x, "DSSTox_ID"))),
    has_MeSH = purrr::map_lgl(pubchem_list, ~ !all(is.na(safe_get(.x, "MeSH_Classification")))),
    other_hazards_n = purrr::map_int(pubchem_list, ~ {
      oh <- safe_get(.x, "Other_Hazards")
      if (all(is.na(oh))) 0L else length(oh)
    }),
    GHS_H = purrr::map(pubchem_list, ~ safe_get(.x, "GHS_Hazards")),
    GHS_Icons = purrr::map(pubchem_list, ~ safe_get(.x, "GHS_Hazard_Icons")),
    Chebi_Roles = chebi_roles[names(pubchem_list)]  # align by CID
  )

  df <- df %>%
    dplyr::mutate(
      Is_Safe_H = purrr::map(GHS_H, ~ { if (all(is.na(.x))) NA else length(intersect(.x, h_exclude)) == 0 }) %>% unlist(),
      Is_Safe_Icons = purrr::map(GHS_Icons, ~ { if (all(is.na(.x))) NA else length(intersect(.x, icon_exclude)) == 0 }) %>% unlist(),
      Is_Safe_Chebi = purrr::map(Chebi_Roles, ~ { if (all(is.na(.x))) NA else length(intersect(.x, exclude_roles)) == 0 }) %>% unlist()
    )

  df
}

# ----------------------------- candidate filter -------------------------------
filter_candidates <- function(tox_df,
                              mw_range = c(200, 500),
                              max_other_hazards = 5) {
  safe_logic <- function(h, i, c) {
    (is.na(h) || h) & (is.na(i) || i) & (is.na(c) || c)
  }
  flags <- mapply(safe_logic, tox_df$Is_Safe_H, tox_df$Is_Safe_Icons, tox_df$Is_Safe_Chebi)
  in_mw <- is.na(tox_df$MW) | (tox_df$MW > mw_range[1] & tox_df$MW < mw_range[2])
  keep <- in_mw & flags & (tox_df$other_hazards_n <= max_other_hazards)
  dplyr::mutate(tox_df, Keep = keep)
}

# ----------------------------- CID-based wrapper ------------------------------
filter_toxicity_by_cids <- function(cids,
                                    all_drug_cids,
                                    json_dir = file.path(tools::R_user_dir("yourpkg","cache"), "pubchem_json"),
                                    download_missing = TRUE,
                                    h_exclude = DEFAULT_H_EXCLUDE,
                                    icon_exclude = DEFAULT_ICON_EXCLUDE,
                                    exclude_roles = DEFAULT_EXCLUDE_ROLES,
                                    mw_range = c(150, 600),
                                    max_other_hazards = 5,
                                    manual_role_patch = MANUAL_ROLE_PATCH,
                                    id_col = c("InputName","DrugName","Display","CID"),
                                    .progress = TRUE) {

  built <- build_pubchem_chebi_table(
    cids             = cids,
    json_dir         = json_dir,
    download_missing = download_missing,
    chebi_cache_dir  = file.path(json_dir, "_chebi_cache"),
    manual_role_patch = manual_role_patch,
    .progress        = .progress
  )

  full <- score_toxicity_table(
    pubchem_list  = built$pubchem,
    chebi_roles   = built$chebi_roles,
    h_exclude     = h_exclude,
    icon_exclude  = icon_exclude,
    exclude_roles = exclude_roles
  )

  full <- full %>%
    dplyr::mutate(CID = as.character(CID)) %>%
    dplyr::left_join(
      all_drug_cids %>% dplyr::mutate(cid = as.character(cid)),
      by = c("CID" = "cid")
    ) %>%
    dplyr::rename(DrugName = drug_name) %>%
    dplyr::mutate(Display = dplyr::coalesce(DrugName, CID))

  out <- filter_candidates(
    tox_df            = full,
    mw_range          = mw_range,
    max_other_hazards = max_other_hazards
  )

  pref <- match.arg(id_col)
  if (!pref %in% names(out)) pref <- "Display"
  filtered_names <- out %>% dplyr::filter(Keep) %>% dplyr::pull(pref) %>% unique()

  list(full_table = full, filtered_table = out, filtered_names = filtered_names)
}

# -------------- Name -> CID resolution (using your map; PubChem fallback) -----
normalize_name <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_replace_all("[[:punct:]]", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()
}

build_name_cid_index <- function(all_drug_cids) {
  stopifnot(all(c("cid","drug_name") %in% names(all_drug_cids)))
  all_drug_cids %>%
    dplyr::transmute(CID = as.character(cid),
                     drug_name = as.character(drug_name),
                     name_norm = normalize_name(drug_name)) %>%
    dplyr::distinct(name_norm, .keep_all = TRUE)
}

pubchem_cid_from_name <- function(name) {
  url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/JSON",
                 URLencode(name, reserved = TRUE))
  resp <- try(httr::GET(url), silent = TRUE)
  if (inherits(resp, "try-error") || httr::status_code(resp) != 200) return(NA_character_)
  j <- try(jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8")), silent = TRUE)
  if (inherits(j, "try-error")) return(NA_character_)
  cid <- try(as.character(j$IdentifierList$CID[[1]]), silent = TRUE)
  if (inherits(cid, "try-error") || length(cid)==0) return(NA_character_)
  cid
}

names_to_cids <- function(drug_names, name_cid_index, use_pubchem_fallback = TRUE, polite_delay = 0.15) {
  df_in <- tibble::tibble(InputName = as.character(drug_names),
                          name_norm = normalize_name(drug_names))
  m <- df_in %>%
    dplyr::left_join(name_cid_index %>% dplyr::select(name_norm, CID, drug_name), by = "name_norm")

  if (use_pubchem_fallback) {
    need <- which(is.na(m$CID))
    if (length(need)) {
      for (i in need) {
        nm <- df_in$InputName[i]
        cid <- pubchem_cid_from_name(nm)
        m$CID[i] <- cid
        Sys.sleep(polite_delay)
      }
    }
  }

  m %>% dplyr::mutate(CID = as.character(CID)) %>% dplyr::distinct(InputName, .keep_all = TRUE)
}

# ------------------------- Names-based wrapper (entrypoint) -------------------
#' @export
filter_toxicity <- function(drug_names,
                            all_drug_cids,
                            json_dir = file.path(tools::R_user_dir("yourpkg","cache"), "pubchem_json"),
                            download_missing = TRUE,
                            id_col = c("InputName","DrugName","Display","CID"),
                            use_pubchem_fallback = TRUE,
                            ...) {
  idx <- build_name_cid_index(all_drug_cids)
  resolved <- names_to_cids(drug_names, idx, use_pubchem_fallback = use_pubchem_fallback)

  name_map <- resolved %>% dplyr::transmute(InputName, CID)
  cids <- resolved$CID[!is.na(resolved$CID)]
  if (!length(cids)) {
    warning("No CIDs resolved from provided names.")
    return(list(resolution = resolved))
  }

  res <- filter_toxicity_by_cids(
    cids             = cids,
    all_drug_cids    = all_drug_cids,
    json_dir         = json_dir,
    download_missing = download_missing,
    id_col           = id_col
  )

  # Reattach original input names & apply name-keyed manual patches
  filtered_tbl <- res$filtered_table %>%
    dplyr::left_join(name_map, by = "CID") %>%
    dplyr::mutate(PreferredName = dplyr::coalesce(InputName, DrugName, Display, CID))

  # Apply MANUAL_ROLE_PATCH by normalized PreferredName
  apply_manual_patch <- function(roles, key, patch_map) {
    add <- patch_map[[key]]
    if (is.null(add)) return(roles)
    unique(c(roles %||% character(), add))
  }
  if (length(MANUAL_ROLE_PATCH)) {
    norm_names <- normalize_name(filtered_tbl$PreferredName)
    filtered_tbl$Chebi_Roles <- mapply(
      function(r, n) apply_manual_patch(r, n, MANUAL_ROLE_PATCH),
      filtered_tbl$Chebi_Roles, norm_names,
      SIMPLIFY = FALSE
    )
    # Recompute Is_Safe_Chebi after patch
    filtered_tbl$Is_Safe_Chebi <- purrr::map_lgl(
      filtered_tbl$Chebi_Roles,
      ~ { if (all(is.na(.x))) NA else length(intersect(.x, DEFAULT_EXCLUDE_ROLES)) == 0 }
    ) %>% unlist()
    # Recompute Keep with the same gate
    safe_logic <- function(h, i, c) (is.na(h)||h) & (is.na(i)||i) & (is.na(c)||c)
    flags <- mapply(safe_logic, filtered_tbl$Is_Safe_H, filtered_tbl$Is_Safe_Icons, filtered_tbl$Is_Safe_Chebi)
    in_mw <- is.na(filtered_tbl$MW) | (filtered_tbl$MW > 200 & filtered_tbl$MW < 500)
    filtered_tbl$Keep <- in_mw & flags & (filtered_tbl$other_hazards_n <= 5)
  }

  preferred_col <- match.arg(id_col)
  filtered_names <- filtered_tbl %>% dplyr::filter(Keep) %>% dplyr::pull(preferred_col) %>% unique()

  list(
    resolution      = resolved,
    full_table      = res$full_table,
    filtered_table  = filtered_tbl,
    filtered_names  = filtered_names
  )
}

#' @export
filter_toxicity_by_cids <- function(cids,
                                    all_drug_cids,
                                    json_dir = file.path(tools::R_user_dir("yourpkg","cache"), "pubchem_json"),
                                    download_missing = TRUE,
                                    h_exclude = DEFAULT_H_EXCLUDE,
                                    icon_exclude = DEFAULT_ICON_EXCLUDE,
                                    exclude_roles = DEFAULT_EXCLUDE_ROLES,
                                    mw_range = c(150, 600),
                                    max_other_hazards = 5,
                                    manual_role_patch = MANUAL_ROLE_PATCH,
                                    id_col = c("InputName","DrugName","Display","CID"),
                                    .progress = TRUE) {

  built <- build_pubchem_chebi_table(
    cids             = cids,
    json_dir         = json_dir,
    download_missing = download_missing,
    chebi_cache_dir  = file.path(json_dir, "_chebi_cache"),
    manual_role_patch = manual_role_patch,
    .progress        = .progress
  )

  full <- score_toxicity_table(
    pubchem_list  = built$pubchem,
    chebi_roles   = built$chebi_roles,
    h_exclude     = h_exclude,
    icon_exclude  = icon_exclude,
    exclude_roles = exclude_roles
  )

  full <- full %>%
    dplyr::mutate(CID = as.character(CID)) %>%
    dplyr::left_join(all_drug_cids %>% dplyr::mutate(cid = as.character(cid)),
                     by = c("CID" = "cid")) %>%
    dplyr::rename(DrugName = drug_name) %>%
    dplyr::mutate(Display = dplyr::coalesce(DrugName, CID))

  out <- filter_candidates(full, mw_range = mw_range, max_other_hazards = max_other_hazards)

  pref <- match.arg(id_col)
  if (!pref %in% names(out)) pref <- "Display"
  filtered_names <- out %>% dplyr::filter(Keep) %>% dplyr::pull(pref) %>% unique()

  list(full_table = full, filtered_table = out, filtered_names = filtered_names)
}

