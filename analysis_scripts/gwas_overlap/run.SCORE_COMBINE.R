# =========================
# Libraries used in this script
# =========================
library(tidyverse)      # dplyr/tibble piping and data manipulation used throughout
library(googlesheets4)  # (optional) write results to Google Sheets

# =========================
# Helper functions
# =========================
source("util.cleaning.R")

# =========================
# I/O parameters (Google Sheet destination)
# =========================
supp_sheet_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/edit?usp=sharing"
new_sheet_name <- "Data S11_COMBINED"


## =====================================
## Load all data (only if not already in the environment)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_anc, etc. as separate objects
}

# =========================
# Build combined evidence table:
# - start from LM/ANOVA results for a specific model
# - join Cox model hazard ratios and p-values
# - compute signed z-scores, support metrics, and combined p-values (with FDR)
# =========================
model_to_use <- "REML.t.val ~ breed.weight.kg + lifespan"
d <- raw_S9_ANOVA_ON_LMER %>% filter(model_used==model_to_use) %>% 
  select(phenotype,model_used,Effect,beta_std.lm,p.value.lm) %>% drop_na()

# Add Cox outputs (exp(coef) and p-value), aligned by phenotype
d <- d %>% left_join(raw_S10_COX_HAZ %>% select(phenotype,exp_coef.cox = `exp(coef)`,p.cox = `Pr(>|z|)`))

# Compute signed z-scores (direction from effect signs), plus combined statistics per Effect
d <- d %>% group_by(Effect) %>% mutate(logHR = log(exp_coef.cox),
    z_lm  = sign(beta_std.lm) * qnorm(p.value.lm / 2, lower.tail = FALSE),
    z_cox = sign(logHR)       * qnorm(p.cox       / 2, lower.tail = FALSE),
    support_score = -z_lm * z_cox,
    z_comb = (-z_lm + z_cox) / sqrt(2),
    p_comb = 2 * pnorm(-abs(z_comb)),
    p_comb_fdr = p.adjust(p_comb, "BH")
  )

# Flag whether each phenotype had Cox results available after the join
d <- d %>% mutate(in_cox_analysis=if_else(is.na(exp_coef.cox),FALSE,TRUE))

# =========================
# Correlation (per Effect) between category-residualized z-scores:
# - adjust z-scores for phenotype category
# - compute Spearman correlation of residuals within each Effect
# =========================
df <- d %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category)) 

df <- df %>% mutate(r_cox=resid(lm(z_cox ~ paper_phenotype_category)),
                             r_lm=resid(lm(z_lm  ~ paper_phenotype_category)))
cor_by_effect <- df %>%
  filter(!is.na(r_cox), !is.na(r_lm)) %>%
  group_by(Effect) %>%
  group_modify(~{
    ct <- cor.test(.x$r_cox, .x$r_lm, method = "spearman")
    tibble(
      n = nrow(.x),
      r = unname(ct$estimate),
      p = ct$p.value,
      cor_str = paste0(
        "spearman rho = ",
        round(unname(ct$estimate), 2),
        ", p = ",
        signif(ct$p.value, 2)
      )
    )
  }) %>%
  ungroup()
cor_by_effect

# =========================
# Final output table: keep key columns and sort by LM signal strength
# =========================
d <- d %>% select(phenotype,model_used,Effect,in_cox_analysis,z_lm,z_cox,support_score,z_comb,p_comb,p_comb_fdr)
d <- d %>% arrange(-abs(z_lm))
