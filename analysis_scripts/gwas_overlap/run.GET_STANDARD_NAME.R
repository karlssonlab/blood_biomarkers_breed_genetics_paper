# Load required libraries used in this script
library(tidyverse)      # Data manipulation (dplyr) and tibbles
library(httr2)          # HTTP requests to OpenAI API
library(jsonlite)       # Parse JSON returned by the model
library(stringr)        # Extract JSON block from model output

## -----------------------------------
## Config: model, API key, and system prompt
## -----------------------------------

model <- "gpt-5.1"  # OpenAI chat-completions model name

api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY is not set. Put it in ~/.Renviron and restart R.")
}

# System prompt instructing the model to conservatively standardize biochemical names
system_prompt_standardize <- "
You are a biochemical nomenclature assistant.

Your task is to return a *standardized metabolite or protein name* ONLY when there is a widely accepted, unambiguous name used in:
- HMDB
- ChEBI
- KEGG
- UniProt (for proteins)
- Common metabolomics shorthand (e.g., TMAO, GSH, ATP)

You must be **conservative**.  
If you are not completely certain that a single standardized name is appropriate, return the **original name unchanged**.

Rules:

1. If the input name already matches a standard biochemical name, return it unchanged.
2. If the input name is a clear synonym of a *widely recognized molecule*, return the standard name.  
   Examples:  
     - \"TMAO\" → \"Trimethylamine N-oxide\"  
     - \"D-Leucic-Acid\" → \"2-hydroxy-4-methylpentanoic acid (D-enantiomer)\"  
3. If the name could refer to multiple compounds, **do not standardize** — return the original input.
4. Do **not infer** stereochemistry, adducts, or structural details unless explicitly present.
5. Do **not expand abbreviations unless they have a globally accepted expansion** (e.g., ATP, GABA).
6. Do **not** rename blood cell types, proteins, clinical labs, or other non-metabolite traits.
7. Keep output extremely short: **return ONLY the standardized name** (or the original name if not safely standardizable).
8. Never explain your reasoning. Never return JSON. Only return the final name string.

Your output must be:
- a single string
- containing only the standardized name or the unchanged original name
"

## -----------------------------------
## Constants and input locations
## -----------------------------------

# P-value thresholds (defined here for consistency with other analyses; not used below)
sigP <- 5e-8
sugP <- 1e-6 

# Google Sheet URL (kept for reference; script uses local RDS below)
supp_data_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/"

## =====================================
## Load all data (from RDS) into the global environment
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # loads objects like raw_S2_ALL_PHENOS into .GlobalEnv
}

# Restrict to clinical analytes and plasma metabolites, keeping unique phenotype labels to standardize
bloodPhenos <- raw_S2_ALL_PHENOS %>%
  filter(paper_phenotype_category %in% c("Clinical analytes", "Plasma metabolites")) %>%
  select(phenotype, plot_label) %>%
  distinct()

## ------------------------------------------------
## Low-level OpenAI request helper (builds the request; does not perform it)
## ------------------------------------------------

openai_standardize_req <- function(trait_name) {
  # User message includes the trait name; request asks for JSON for robust parsing downstream
  user_msg <- paste0(
    "Input trait/metabolite name:\n",
    "trait1: ", trait_name, "\n\n",
    "Return ONLY a JSON object with fields {\"original\", \"standard_name\"}."
  )
  
  request("https://api.openai.com/v1/chat/completions") |>
    req_auth_bearer_token(api_key) |>
    req_body_json(list(
      model = model,
      temperature = 0,
      messages = list(
        list(role = "system", content = system_prompt_standardize),
        list(role = "user",   content = user_msg)
      )
    ))
}

## ------------------------------------------------
## Parse JSON from model content into a single standardized name
## ------------------------------------------------

parse_standard_name <- function(content, trait_name) {
  # Extract the first {...} block in case the model returns extra text
  json_text <- str_extract(content, "\\{[\\s\\S]*\\}")
  if (is.na(json_text)) {
    stop("Model did not return JSON for trait: ", trait_name,
         "\nFull content:\n", content)
  }
  
  obj <- tryCatch(
    jsonlite::fromJSON(json_text, simplifyVector = FALSE),
    error = function(e) {
      stop("Could not parse JSON for trait: ", trait_name,
           "\nJSON text:\n", json_text,
           "\nError: ", conditionMessage(e))
    }
  )
  
  if (!is.list(obj)) {
    stop("Parsed JSON is not an object for trait: ", trait_name)
  }
  
  # Enforce expected schema from the model response
  if (!all(c("original", "standard_name") %in% names(obj))) {
    stop("JSON missing required fields for trait: ", trait_name,
         "\nNames found: ", paste(names(obj), collapse = ", "))
  }
  
  original      <- as.character(obj[["original"]])
  standard_name <- as.character(obj[["standard_name"]])
  
  # Defensive fallbacks: if missing/empty or not scalar, keep the original input
  if (is.na(standard_name) || standard_name == "") {
    standard_name <- trait_name
  }
  if (length(standard_name) != 1) {
    standard_name <- trait_name
  }
  
  standard_name
}

## ------------------------------------------------
## Standardize a single trait (request -> response -> parse)
## ------------------------------------------------

standardize_one_trait <- function(trait_name) {
  # Early return for missing/blank labels
  if (is.na(trait_name) || trimws(trait_name) == "") return(trait_name)
  
  resp <- openai_standardize_req(trait_name) |>
    req_perform()
  
  body <- resp_body_json(resp, simplifyVector = FALSE)
  
  # Surface API-side errors with the trait name for easier debugging
  if (!is.null(body$error)) {
    stop("OpenAI API error while standardizing '", trait_name, "': ",
         body$error$message)
  }
  
  content <- body[["choices"]][[1]][["message"]][["content"]]
  parse_standard_name(content, trait_name)
}

## ------------------------------------------------
## Vectorized wrapper: standardize many traits with progress messages
## ------------------------------------------------
standardize_traits <- function(trait_vec) {
  if (!is.character(trait_vec)) {
    stop("trait_vec must be a character vector.")
  }
  
  n <- length(trait_vec)
  message("Standardizing ", n, " traits...")
  
  purrr::map_chr(seq_along(trait_vec), function(i) {
    trait <- trait_vec[i]
    message("  [", i, "/", n, "] ", trait)
    standardize_one_trait(trait)
  })
}

# Extract the unique plot labels to send for standardization
traits <- bloodPhenos %>% select(plot_label) %>% distinct() %>% pull(plot_label)

## ------------------------------------------------
## Run standardization and join results back to the full phenotype table
## ------------------------------------------------

standard_names <- standardize_traits(traits)
results <- tibble(plot_label = traits,
                standard_name = standard_names)
new_bloodPhenos <- raw_S2_ALL_PHENOS %>% left_join(results)
