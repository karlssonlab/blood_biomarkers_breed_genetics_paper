# ============================================================
# Libraries and helper functions
# ============================================================
library(tidyverse)      # Data manipulation, modeling workflows, and visualization (dplyr/tidyr/purrr/stringr/etc.)
library(rstatix)        # anova_test() and get_anova_table() helpers

# Project-specific data cleaning utilities (e.g., breed name standardization)
source("util.cleaning.R")

# ============================================================
# Output configuration (Google Sheet target identifiers)
# ============================================================
supp_sheet_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/edit?usp=sharing"
out_sheet <- "Data S9_ANOVA_ON_LMER"

# ============================================================
# Load all data (only if not already present in the environment)
# ============================================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # exports objects from the RDS list into the global environment (raw_pheno, raw_anc, etc.)
}

# ============================================================
# Prepare breed-level covariates (weight, lifespan) and clean fields
# ============================================================
breeds <- raw_S12_BREED_PHENOS %>% select(breed,breed.weight.kg,lifespan=`median_age_death.(McMillan.2024)`) %>% distinct()
breeds <- standardize_breed_names(breeds,breed)
breeds <- convert_yesno_to_logical(breeds)

# ============================================================
# Build ANOVA dataset: join breed covariates to per-phenotype LM/LMER t-values
# ============================================================
d_anova <- standardize_breed_names(raw_S8_LMER,breed) %>% select(phenotype,breed,REML.t.val)
d_anova <- breeds %>% inner_join(d_anova)
dim(d_anova)
dim(d_anova)
d_anova <- d_anova %>% drop_na()

# Model specification (kept as both formula object and string for reporting)
formula <- as.formula("REML.t.val ~ breed.weight.kg + lifespan")
formula_str <- deparse(formula)

# ============================================================
# Per-phenotype ANOVA, then FDR adjust per model term across all phenotypes
# ============================================================
anova_result <- d_anova %>% 
  left_join(raw_S2_ALL_PHENOS %>% select(phenotype) %>% distinct()) %>% 
  group_by(phenotype) %>%          
  anova_test(formula) %>%      
  get_anova_table() %>% ungroup() %>% group_by(Effect) %>%
  mutate(p_adj_fdr = p.adjust(p, method = "BH")) %>% ungroup()

# ============================================================
# Per-phenotype linear model with standardized variables (effect sizes comparable across traits)
# Extract coefficients and confidence intervals for each term
# ============================================================
lm_result <- d_anova %>%
  group_by(phenotype) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(
      scale(REML.t.val) ~ scale(breed.weight.kg) + scale(lifespan),
      data = .x
    )),
    coefs = map(fit, ~ {
      cf <- as.data.frame(summary(.x)$coefficients) %>%
        tibble::rownames_to_column("term")
      ci <- as.data.frame(confint(.x)) %>%
        tibble::rownames_to_column("term") %>%
        rename(conf.low = `2.5 %`, conf.high = `97.5 %`)
      left_join(cf, ci, by = "term")
    })
  ) %>%
  unnest(coefs) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Convert term names like "scale(breed.weight.kg)" back to "breed.weight.kg"
    term = term %>%
      stringr::str_remove("^scale\\(") %>%
      stringr::str_remove("\\)$") %>%
      stringr::str_remove("TRUE$")
  ) %>%
  transmute(
    phenotype,
    term,
    beta_std   = Estimate,
    std.error  = `Std. Error`,
    statistic  = `t value`,
    p.value    = `Pr(>|t|)`,
    conf.low,
    conf.high
  ) %>%
  ungroup()
# Align LM term names with ANOVA 'Effect' field for merging
lm_result <- lm_result %>% mutate(Effect = str_remove(term, "TRUE$"))
by_cols <- c("phenotype", "Effect")

# ============================================================
# Merge ANOVA + LM outputs into a single table with suffixes to disambiguate columns
# ============================================================
merged_results <- full_join(
  anova_result %>%
    mutate(model_used = formula_str) %>%
    rename_with(~ paste0(.x, ".anova"), -all_of(c(by_cols, "model_used"))),
  lm_result %>%
    mutate(model_used = formula_str) %>%
    rename_with(~ paste0(.x, ".lm"), -all_of(c(by_cols, "model_used"))),
  by = c(by_cols, "model_used")
) %>%
  select(-`p<.05.anova`) %>%
  select(-matches("\\.x$|\\.y$")) %>%
  relocate(model_used, .before = 1) %>%
  arrange(across(all_of(c("phenotype", "Effect"))))
