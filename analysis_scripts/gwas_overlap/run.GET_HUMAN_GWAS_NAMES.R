## ================================
##  Libraries (only those used below)
## ================================
library(tidyverse)   # dplyr/tibble/readr/stringr pipes and helpers
library(httr)        # OpenAI API POST requests
library(jsonlite)    # JSON encode/decode for OpenAI requests/responses
library(glue)        # Build LLM prompts with interpolation
library(googlesheets4)  # Append results to Google Sheets

## ================================
##  Inputs / constants
## ================================
alt_url <- "https://docs.google.com/spreadsheets/d/1JmhdoWwaedZUDA1azA7a7QQGUEC3pQ95s5XrxPM21UM/"

# Allowlist of PubMed IDs to keep from the GWAS Catalog
ok_pubmedIDs <- c(
  35213538, 38448586, 34503513, 33462484, 32888494, 39024449,
  34594039, 29084231, 32888493, 36357675, 35668104, 34648354,
  35347128, 36635386, 36764567, 37253714, 33031748, 27863252,
  35050183, 30595370, 35995766, 36168886, 38626723, 38116116, 38658550
)

# Genome-wide significance threshold
sigP <- 5e-8

# Local GWAS Catalog associations file (TSV)
in_gwas_catalog <- "../../data/human_dog_GWAS_intersect/gwas_catalog_v1.0-associations_e113_r2025-02-08.tsv"


## ================================
##  OpenAI config
## ================================
openai_model <- "gpt-4.1-mini"   # change if you want

# Read API key from environment; fail fast if missing
api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY is not set. Put it in ~/.Renviron and restart R.")
}

# System prompt: classify whether a trait is a direct BLOOD measurement and standardize the name
system_prompt_standardize <- "
You classify trait names measured in BLOOD and return a canonical name and include=TRUE/FALSE.

WORKFLOW:
1) First decide if INCLUDE should be TRUE or FALSE.
2) If INCLUDE = FALSE, do NOT spend effort refining the canonical name; set canonical_name equal to the input and give a short reason.
3) If INCLUDE = TRUE, then generate a cleaned, standardized canonical_name.

INCLUDE = TRUE only if the trait is a direct blood measurement:
  - Serum/plasma metabolites (small molecules, amino acids, lipids, acylcarnitines)
- Clinical chemistry analytes (creatinine, urea, glucose, calcium, bilirubin, ALT/AST, etc.)
- CBC traits (RBC, WBC, neutrophils, eosinophils, platelets, hematocrit, MCV, etc.)
- Serum/plasma proteins, cytokines, complement factors, hormones

INCLUDE = FALSE if:
  - Urine/CSF/tissue/fecal measurements
- Disease-specific phrases (\\\"in chronic kidney disease\\\", \\\"in elite athletes\\\")
  - Traits not measured in blood (height, diabetes, bone density, psychiatric traits)
  - Unknown feature IDs (X-11787, QI7389)
  - Environmental, behavioral, questionnaire-based traits

CANONICAL NAME RULES (only apply when INCLUDE = TRUE):
  - Normalize to the underlying biochemical/clinical entity.
  - Remove prefixes/suffixes: \\\"Serum\\\", \\\"Plasma\\\", \\\"Metabolite levels\\\", \\\"levels\\\", \\\"concentrations\\\".
  - Remove transformations: \\\"(minimum)\\\", \\\"(maximum)\\\", \\\"(mean)\\\", \\\"inv-norm transformed\\\".
  - Remove disease qualifiers (\\\"in CKD\\\" etc.).
  - For complex lipids: remove chain-length/unsaturation annotations such as (38:3), (18:1/20:4), O-(36:2); keep only the lipid class.
    Examples:
      \\\"Phosphatidylserine (38:3)\\\" → \\\"phosphatidylserine\\\"
      \\\"PC(20:4/18:1)\\\" → \\\"phosphatidylcholine\\\"
      \\\"TG(18:0/18:1/16:0)\\\" → \\\"triacylglycerol\\\"
      \\\"LPC 18:2\\\" → \\\"lysophosphatidylcholine\\\"

  - Protein naming rule:
    Collapse isoforms, subfamilies, and member numbers to the core family name.
    Examples:
      \\\"leukocyte immunoglobulin-like receptor subfamily A member 5\\\" → \\\"leukocyte immunoglobulin-like receptor\\\"
      \\\"interleukin-1 receptor type II precursor\\\" → \\\"interleukin-1 receptor\\\"
      \\\"apolipoprotein E isoform 2\\\" → \\\"apolipoprotein E\\\"

  - Keep chemical names unchanged.
  - If ambiguous but clearly a blood measurement, choose the simplest reasonable canonical name.

SPECIAL CASE WHEN INCLUDE = FALSE:
  - Set canonical_name = input (no additional simplification).
  - Reason should explain why it is excluded (e.g., non-blood trait, urine measurement, disease-specific cohort, unknown feature ID).

OUTPUT FORMAT:
Return ONE JSON object:
{
  \\\"results\\\": [
    {\\\"input\\\": \\\"...\\\", \\\"canonical_name\\\": \\\"...\\\", \\\"include\\\": true/false, \\\"reason\\\": \\\"...\\\"}
  ]
}
"



## ================================
##  OpenAI call helper
## ================================
# One-shot chat completion call; returns parsed JSON (R list)
call_openai_once <- function(user_prompt,
                             api_key = Sys.getenv("OPENAI_API_KEY"),
                             model = openai_model,
                             system_prompt = system_prompt_standardize) {
  
  resp <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    add_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ),
    body = jsonlite::toJSON(
      list(
        model = model,
        messages = list(
          list(role = "system", content = system_prompt),
          list(role = "user",   content = user_prompt)
        ),
        temperature = 0
      ),
      auto_unbox = TRUE
    )
  )
  
  # Hard fail on non-200 so we can see the server error payload
  if (httr::status_code(resp) != 200) {
    txt_err <- httr::content(resp, as = "text", encoding = "UTF-8")
    stop("OpenAI API error (status ", httr::status_code(resp), "):\n", txt_err)
  }
  
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  
  # Ensure expected structure is present before extracting the message content
  if (is.null(parsed$choices) || length(parsed$choices) < 1) {
    stop("No 'choices' field found in OpenAI response:\n", txt)
  }
  
  content_str <- parsed$choices[[1]]$message$content
  
  out <- jsonlite::fromJSON(content_str)
  out
}

## ================================
##  Canonicalization for a vector
## ================================
# Send a vector of traits and return a tibble with input/canonical/include/reason
canonicalize_blood_traits <- function(traits_vec) {
  
  trait_lines <- paste0("- ", traits_vec, collapse = "\n")
  
  user_prompt <- glue("
You are given a list of trait names (one per line) under 'TRAITS:'.
For EACH trait, return an entry in the `results` array as:
  {{\"input\": \"...\", \"canonical_name\": \"...\", \"include\": true/false, \"reason\": \"...\"}}.

TRAITS:
{trait_lines}
")
  
  res <- call_openai_once(user_prompt)
  
  if (!"results" %in% names(res)) {
    stop(
      "OpenAI response missing 'results' field. Got:\n",
      paste(capture.output(str(res)), collapse = "\n")
    )
  }
  
  as_tibble(res$results)
}

## ================================
##  Simple batch helper
## ================================
# Return index batches for chunked API calls
split_into_batches <- function(x, batch_size = 80) {
  n <- length(x)
  idx <- seq_len(n)
  split(idx, ceiling(idx / batch_size))
}

## ================================
##  Clean trait names before LLM
## ================================
# Pre-clean GWAS Catalog trait strings to reduce prompt noise (keep meaning intact)
clean_trait_names <- function(df) {
  df %>%
    mutate(
      trait_clean = gwas_catalog_trait_name,
      
      # Strip 'Metabolite levels (...)' wrapper and keep the inside text
      trait_clean = str_replace(
        trait_clean,
        "^Metabolite levels \\((.*)\\)$",
        "\\1"
      ),
      trait_clean = str_replace(
        trait_clean,
        "^Metabolite levels \\(([^)]*)\\)",
        "\\1"
      ),
      
      # Remove parenthetical notes about summary/transformations
      trait_clean = str_replace_all(
        trait_clean,
        "\\s*\\((?=[^)]*(maximum|minimum|mean|inv-?norm|transformed)[^)]*)[^)]*\\)",
        ""
      ),
      
      trait_clean = str_squish(trait_clean)
    ) %>% 
      mutate(
        # Normalize a few generic words and plurals that add noise for matching
        trait_clean = trait_clean %>%
          str_replace_all("\\bcounts?\\b", "") %>%   # remove "count" or "counts"
          str_replace_all("\\blevels?\\b", "") %>%   # remove "level" or "levels"
          str_replace_all("\\bcells\\b", "cell") %>% # plural → singular
          str_squish()                               # clean extra whitespace
      )
  
}

## ================================
##  Optional post-processing: collapse over-specific canonical names
## ================================
# Apply regex-based simplifications after LLM canonicalization (currently not called below)
simplify_canonical <- function(x) {
  x %>%
    # 1. Cholesteryl / cholesterol ratios & subfractions
    str_replace("^cholesteryl esters? to total lipids ratio.*$", "cholesteryl ester to total lipids ratio") %>%
    str_replace("^cholesteryl ester in .*", "cholesteryl ester") %>%
    str_replace("^cholesterol to total lipids ratio in .*", "cholesterol to total lipids ratio") %>%
    str_replace("^free cholesterol to total lipids ratio in .*", "free cholesterol to total lipids ratio") %>%
    
    # 2. Lipoprotein particle concentration / diameter
    str_replace("^.* HDL particle concentration$", "HDL particle concentration") %>%
    str_replace("^.* LDL particle concentration$", "LDL particle concentration") %>%
    str_replace("^.* VLDL particle concentration$", "VLDL particle concentration") %>%
    str_replace("^average diameter for HDL particles$", "HDL particle diameter") %>%
    str_replace("^average diameter for LDL particles$", "LDL particle diameter") %>%
    str_replace("^average diameter for VLDL particles$", "VLDL particle diameter") %>%
    
    # 3. Family-level protein collapses (optional, if you want pathway-level)
    str_replace("^cadherin-.*$", "cadherin") %>%
    str_replace("^coagulation factor [IVXLC]+$", "coagulation factor") %>%
    str_replace("^c-c motif chemokine .*", "CC chemokine") %>%
    str_replace("^c-x-c motif chemokine .*", "CXC chemokine") %>%
    str_replace("^interleukin-[0-9A-Z]+( receptor.*)?$", "interleukin / receptor") %>%
    str_replace("^leucine-rich repeat.*$", "leucine-rich repeat protein") %>%
    str_replace("^immunoglobulin .*$", "immunoglobulin")
}


## ----------------------------
## 1. Load and filter GWAS Catalog to a trait list
## ----------------------------

# Load the catalog once into the session if not already present
if (!exists("raw_gwas_human") | !exists("raw_human_trait_list")) {
  raw_gwas_human <- read_tsv(in_gwas_catalog, show_col_types = FALSE)
}

# Keep genome-wide significant hits with chromosome info and matching PubMed IDs
filtered_gwas_human <- raw_gwas_human %>%
  filter(`P-VALUE` <= sigP & !is.na(CHR_ID)) %>%
  select(-SNP_ID_CURRENT) %>%
  filter(PUBMEDID %in% ok_pubmedIDs)

dim(filtered_gwas_human)

# Reduce to unique GWAS Catalog trait strings
long_gwas_names <- filtered_gwas_human %>%
    select(gwas_catalog_trait_name = `DISEASE/TRAIT`) %>%
    distinct()

dim(long_gwas_names)

# Create cleaned trait strings used for LLM prompting/join-back
long_gwas_names <- clean_trait_names(long_gwas_names)

## ----------------------------
## 2. Canonicalize trait names via LLM and append results to Google Sheet
## ----------------------------

# Unique cleaned traits for batching
gwas_names <- long_gwas_names %>% select(trait_clean) %>% distinct() %>% arrange(trait_clean)

set.seed(2)
batch_size <- 5

traits <- gwas_names$trait_clean
idx_batches <- split_into_batches(traits, batch_size = batch_size)
  
all_results <- list()
counter <- 0

# Iterate through batches, classify/canonicalize, then append to the target sheet
for (batch_i in seq_along(idx_batches)) {
  idx <- idx_batches[[batch_i]]
  message("\n=== Batch ", batch_i,
          " (traits ", min(idx), "–", max(idx), ") ===")
  results_batch <- canonicalize_blood_traits(traits[idx])
  results_out <- results_batch %>% select(trait_clean=input,canonical_name,include,reason) %>% 
      inner_join(long_gwas_names) %>% select(-trait_clean)
  results_out <- results_out %>% select(gwas_catalog_trait_name,canonical_name,include,reason)
  results_out <- results_out %>%
    mutate(
      # If excluded, blank out canonical_name so downstream uses only included canonicals
      canonical_name = if_else(include == FALSE, "", canonical_name)
      )
  sheet_append(
    data = results_out,
    ss = alt_url,
    sheet = "Data S_HUMAN_TRAITS" ### this output is in DataDryad repository in S_HUMAN_TRAITS.tsv
  )
}
