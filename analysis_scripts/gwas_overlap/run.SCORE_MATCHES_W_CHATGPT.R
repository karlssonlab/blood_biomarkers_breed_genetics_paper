# ------------------------------
# Libraries used in this script
# ------------------------------
library(tidyverse)      # dplyr/purrr/tibble utilities used throughout
library(googlesheets4)  # read_sheet() / sheet_write() for Google Sheets I/O
library(httr2)          # HTTP requests to the OpenAI API
library(jsonlite)       # JSON parsing of model output
library(stringr)        # str_extract() for pulling JSON out of model content
library(glue)           # glue() templating for prompts

## -----------------------------------
## Config: model + credentials + scoring rubric
## -----------------------------------

model <- "gpt-5.1"  # Model name for chat/completions

# Read API key from environment and hard-stop if missing
api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY is not set. Put it in ~/.Renviron and restart R.")
}

# System prompt defining the strict 1–5 scoring scheme and JSON-only output contract
system_prompt <- "You score pairs of **blood-measured biological traits** on a strict 1–5 scale.

1 = SAME MEASURE  
    • Exact same molecule/entity OR strict subset/superset  
    • Identical text (case/whitespace-insensitive) → 1  
    • Aggregate vs component of same measure (e.g., total WBC ↔ neutrophils; BUN ↔ urea)

2 = SAME SPECIFIC PATHWAY / DIRECT COMPONENT  
    Assign 2 when both traits belong to the same defined metabolic or blood-cell pathway, including:  
      • Substrate–product or ≤2 steps apart  
      • Same canonical pathway module: glycolysis, TCA, PPP, urea cycle, β-oxidation, carnitine shuttle,  
        bile acid synthesis/transport, amino acid catabolism (incl. BCAA), glutathione cycle,  
        tryptophan–kynurenine/indole, 1-carbon/methylation cycle, heme–bilirubin, purine/pyrimidine  
      • Enzyme ↔ substrate/cofactor  
      • Metabolite ↔ its conjugate  
      • Metabolite ↔ direct clinical chemistry readout  
      • Direct lineage composites (total WBC ↔ lymphocytes; RBC ↔ reticulocytes)  
      • Metabolite ↔ its primary carrier or transport protein in blood (e.g., heptadecanoic acid ↔ albumin / ALB)

3 = FUNCTIONALLY RELATED (NOT pathway neighbors)  
    • Same biological process: oxidative stress, immune activation, thrombosis, gut–liver axis,  
      amino acid homeostasis, RBC turnover  
    • Same biological pathway or same mechanistic pathway but not close enough to score as 2

4 = SYSTEMIC / STATE CORRELATION  
    • Co-vary due to shared physiological state: renal/hepatic function, inflammation,  
      cardiometabolic status, nutritional state, disease severity  
    • No mechanistic pathway connection

5 = NO CLEAR RELATIONSHIP

Quick checks:  
    1 → identical or strict component/superset  
    2 → same defined pathway module OR enzyme–substrate/cofactor OR direct clinical readout OR primary carrier–metabolite pair in blood  
    3 → same process but not pathway-linked  
    4 → systemic correlation only  
    5 → none of the above

Return ONLY JSON with: trait1, trait2, score, rationale."

## -----------------------------------
## Constants and source locations
## -----------------------------------

# Significance thresholds for filtering GWAS overlap regions / human trait associations
sigP <- 5e-8
sugP <- 1e-6 

# Google Sheets sources
#supp_data_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/"
alt_url <- "https://docs.google.com/spreadsheets/d/1JmhdoWwaedZUDA1azA7a7QQGUEC3pQ95s5XrxPM21UM/"


# Create empty table with desired columns
empty_df <- tibble(
  standard_name = character(),
  human_trait_name = character(),
  score = numeric(),
  rationale = character()
)

# Check if sheet exists
if (!("Data S_TRAIT_MATCH" %in% sheet_names(alt_url))) {
  
  sheet_write(
    ss = alt_url,
    data = empty_df,
    sheet = "Data S_TRAIT_MATCH"
  )
  
}

# ------------------------------------------------------------------------------
# Load data objects (only if not already present in the current R environment)
# ------------------------------------------------------------------------------

if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv) 
}

## -----------------------------------
## Load input tables (only if not already in memory)
## -----------------------------------

if (!exists("raw_gwas_overlap")){
raw_gwas_overlap <- as_tibble(read.delim("../../data/human_dog_GWAS_intersect/S_GWAS_RAW_OVERLAP.tsv",na.strings=c("NA",""),header=T))
#  raw_S2_ALL_PHENOS <- as_tibble(read_sheet(supp_data_url, sheet = "Table S2_ALL_PHENOS", na = "NA"))
}

# Restrict to blood-measured phenotypes of interest and keep label/name columns used downstream
bloodPhenos <- raw_S2_ALL_PHENOS %>% filter(paper_phenotype_category %in% c("Clinical analytes","Plasma metabolites")) %>% 
  select(phenotype,plot_label,standard_name) %>% distinct()

## -----------------------------------
## Pair selection: choose which trait pairs still need scoring
## -----------------------------------

make_pairs <- function() {
  # Use global `all` (prepared later) and keep a local copy for transformations
  all_F <- all
  matchingF <- as_tibble(read_sheet(alt_url, sheet = "Data S_TRAIT_MATCH", na = "NA"))
  matchingF <- matchingF %>% mutate(standard_name=as.character(standard_name),
                                    human_trait_name=as.character(human_trait_name),
                                    score=as.numeric(score),
                                    rationale=as.character(rationale))
  # Identify already-scored pairs (excluding generic human traits) with strong matches (score <= 2)
  foundF <- all_F %>% left_join(matchingF,join_by(standard_name, human_trait_name))
  foundF <- foundF %>% filter(!is.na(score)) %>% filter(human_trait_name!="metabolite"&human_trait_name!="cholesterol") %>% group_by(standard_name,region) %>% summarize(score=min(score))
  foundF <- foundF %>% filter(score<=2) 
  
  # Candidate pairs are those not yet present in the match sheet and not already "covered" by a strong match in the same region
  missingF <- all_F %>% anti_join(matchingF %>% select(standard_name,human_trait_name),join_by(standard_name, human_trait_name))
  missingF <- missingF %>% anti_join(foundF %>% select(standard_name,region),join_by(standard_name,region)) %>% distinct()
  
  # Collapse to min human P per (dog trait, human trait) to prioritize strongest associations
  missingF <- missingF %>% group_by(standard_name,human_trait_name) %>% summarize(human_P=min(human_P)) %>% ungroup() 
  print(paste0("Missing=",nrow(missingF)))
  
  # Build a set of candidate human-P cutoffs derived from rounded -log10(P)
  range <- missingF %>% mutate(human_logP=round(-log10(human_P),0)) %>% select(human_logP) %>% 
    arrange(-human_logP) %>% distinct() %>% mutate(cutoff=10**-human_logP) %>% select(-human_logP)
  
  # For each cutoff, count how many distinct pairs would be included
  counts <- crossing(range,missingF) %>% filter(human_P<cutoff) %>% 
    select(-human_P) %>% distinct() %>% group_by(cutoff) %>% count() %>% arrange(n)
  
  # Pick the most stringent cutoff (smallest P) that yields at least 1 row; fallback to least stringent if needed
  chosen <- counts %>%
    filter(n >= 1) %>%
    arrange(cutoff) %>%
    head(n=1)
  
  if (nrow(chosen) == 0) {
    chosen <- counts %>%
      arrange(-cutoff) %>%
      head(n=1)
  }  
  if (nrow(chosen) == 0) {
    stop("All needed pairs scored.")
  }
  
  print(paste0("p cutoff is ",chosen$cutoff))
  print(paste0("pairs for comparison = ",chosen$n))
  
  # Limit to up to 10 pairs per batch to keep API calls bounded per loop iteration
  torun <- missingF %>% filter(human_P < chosen$cutoff) %>% select(trait1=standard_name,trait2=human_trait_name) %>% distinct() %>% slice_head(n=10)
  return(torun)
}

## -----------------------------------
## OpenAI request builder (chat/completions)
## -----------------------------------

openai_chat_req <- function(user_content) {
  request("https://api.openai.com/v1/chat/completions") |>
    req_auth_bearer_token(api_key) |>
    req_timeout(20) |>  # Request timeout in seconds
    req_body_json(list(
      model = model,
      temperature = 0,
      messages = list(
        list(role = "system", content = system_prompt),
        list(role = "user",   content = user_content)
      )
    ))
}

## -----------------------------------
## Parse JSON from model output (defensive extraction + validation)
## -----------------------------------

parse_json_from_content <- function(content, trait1, trait2) {
  # Log raw content for debugging when JSON parsing fails
  cat("\n--- RAW MODEL CONTENT ---\n", content, "\n--------------------------\n")
  
  # Extract the first JSON object-like block from the text response
  json_text <- str_extract(content, "\\{[\\s\\S]*\\}")
  if (is.na(json_text)) {
    stop("Model did not return JSON for pair: ", trait1, " | ", trait2)
  }
  
  # Parse without simplification to keep strict control over types/fields
  out <- tryCatch(
    jsonlite::fromJSON(json_text, simplifyVector = FALSE),
    error = function(e) {
      stop("Could not parse JSON for pair: ", trait1, " | ", trait2,
           "\nJSON text:\n", json_text,
           "\nError: ", conditionMessage(e))
    }
  )
  
  if (!is.list(out)) {
    stop("Parsed JSON is not an object for pair: ", trait1, " | ", trait2,
         "\nstr(out):\n", paste(capture.output(str(out)), collapse = "\n"))
  }
  
  # Require all contract fields to be present
  required <- c("trait1", "trait2", "score", "rationale")
  missing <- required[!required %in% names(out)]
  if (length(missing) > 0) {
    stop("JSON missing fields: ", paste(missing, collapse = ", "),
         " for pair: ", trait1, " | ", trait2,
         "\nNames found: ", paste(names(out), collapse = ", "),
         "\nJSON text:\n", json_text)
  }
  
  # Normalize to a one-row tibble with expected types
  tibble(
    trait1    = as.character(out[["trait1"]]),
    trait2    = as.character(out[["trait2"]]),
    score     = as.integer(out[["score"]]),
    rationale = as.character(out[["rationale"]])
  )
}

## -----------------------------------
## One API call for a single pair
## -----------------------------------

score_pair_once <- function(trait1, trait2) {
  
  # Construct user prompt for this trait pair; system prompt enforces scoring rubric + JSON-only response
  user_prompt <- glue("
Score this biological trait pair using the strict 1–5 scheme.

trait1: {trait1}
trait2: {trait2}

Return ONLY a single JSON object with fields: trait1, trait2, score, rationale.
")
  
  resp <- openai_chat_req(user_prompt) |>
    req_perform()
  
  body <- resp_body_json(resp, simplifyVector = FALSE)
  
  # If the API returned a structured error payload, stop with its message
  if (!is.null(body$error)) {
    stop("OpenAI API error for pair ", trait1, " | ", trait2, ": ",
         body$error$message)
  }
  
  # Extract assistant content from the first choice
  content <- body[["choices"]][[1]][["message"]][["content"]]
  
  parse_json_from_content(content, trait1, trait2)
}

## -----------------------------------
## Retry wrapper around score_pair_once (with exponential backoff)
## -----------------------------------

score_pair <- function(trait1, trait2,
                       max_tries   = 5,
                       base_delay  = 2) {
  attempt <- 1
  
  repeat {
    res <- tryCatch(
      score_pair_once(trait1, trait2),
      error = identity
    )
    
    # Success: a tibble result was returned
    if (inherits(res, "tbl_df")) {
      return(res)
    }
    
    # Failure: evaluate whether this error is retryable
    err_msg <- conditionMessage(res)
    retryable <- FALSE
    status <- NA_integer_
    
    # HTTP-level errors from httr2
    if (inherits(res, "httr2_http_error")) {
      status <- res$response$status_code
      if (status %in% c(408, 429, 500, 502, 503, 504)) {
        retryable <- TRUE
      }
    } else {
      # Non-HTTP errors: common transient network/timeout patterns
      if (grepl("Timeout was reached", err_msg, ignore.case = TRUE) ||
          grepl("timed out after", err_msg, ignore.case = TRUE)   ||
          grepl("resolve host", err_msg, ignore.case = TRUE)) {
        retryable <- TRUE
      }
    }
    
    # Abort if not retryable or we have exhausted attempts
    if (!retryable || attempt >= max_tries) {
      stop(
        "Failed for pair ", trait1, " | ", trait2,
        " after ", attempt, " attempt(s).",
        if (!is.na(status)) paste0(" HTTP status: ", status, "."),
        "\nLast error: ", err_msg
      )
    }
    
    # Exponential backoff between retries
    delay <- base_delay * 2^(attempt - 1)
    message(
      "Attempt ", attempt, " failed for pair ", trait1, " | ", trait2,
      if (!is.na(status)) paste0(" (HTTP ", status, ")") else "",
      ". Retrying in ", delay, "s...\n  Error: ", err_msg
    )
    
    Sys.sleep(delay)
    attempt <- attempt + 1
  }
}

## -----------------------------------
## Batch scoring helper: score a data frame of (trait1, trait2) pairs
## -----------------------------------

score_pairs_df <- function(pairs_df) {
  stopifnot(all(c("trait1", "trait2") %in% names(pairs_df)))
  
  # Split into 1-row tibbles to keep each API call isolated and recoverable
  pairs_df |>
    mutate(row_id = row_number()) |>
    group_split(row_id) |>
    map_dfr(~ {
      t1_original <- .x$trait1[1]
      t2_original <- .x$trait2[1]
      score_pair(t1_original, t2_original)
    })
}

## -----------------------------------
## Build the full candidate set of dog-region vs human-trait pairs to potentially score
## -----------------------------------

# Filter overlap table to significant human associations and suggestive dog region P-values
all <- raw_gwas_overlap %>% filter(lifted_over&!is.na(human_P) & human_P <= sigP & region_P <= sugP) %>% select(phenotype,region,region_P,human_trait_name,human_P)
all <- all %>% group_by(phenotype,region,region_P,human_trait_name) %>% summarize(human_P=min(human_P))
all <- bloodPhenos %>% right_join(all) %>% 
  select(standard_name,region,region_P,human_trait_name,human_P) %>% distinct()
all <- all %>% arrange(-region_P,-human_P)

## -----------------------------------
## Main loop: repeatedly pick unscored pairs, score them, and write back to Google Sheet
## -----------------------------------

for (i in 1:5){
  # Refresh current matching table from the sheet each iteration to avoid overwriting remote edits
  matching <- as_tibble(read_sheet(alt_url, sheet = "Data S_TRAIT_MATCH", na = "NA"))
  
  # Select next batch of pairs to score and call the API
  pairs_df <- make_pairs()
  results <- score_pairs_df(pairs_df) 
  print(dim(matching))
  sheet_write(
    ss = alt_url,
    data = matching,
    sheet = "Data S_TRAIT_MATCH" ### this output is in S_TRAIT_MATCH.tsv in data Dryad
  )
  
  # Drop large objects between iterations
  rm(pairs_df)
  rm(matching)
}
