library(tidyverse)
library(survival)
library(coxme)
library(lubridate)
library(survminer)
library(googlesheets4)
library(psych)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(ggforestplot)
library(broom)

# 937 dogs with all metabolite and CBC/Chem data analyzed (159 traits)
##### read in intermediate phenotypes analyzed in study #####
if (!exists("dat")) {
  dat <- readRDS("DAP_supp_data.rds") # supplemental material tables and datafiles located in DataDryad
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_S15_BREED_ANC, etc.
}

intermed_phenos <- raw_S2_ALL_PHENOS %>% select(phenotype) %>% pull()

intermed_phenos <- gsub("_ln_transformed", "", intermed_phenos)
intermed_phenos <- gsub("_sqrt_transformed", "", intermed_phenos)

##### read in Ben's RData file with longitudinal data ####
#file can be found in DataDryad within longitudinal_Rdata_files_for_survival_analysis
load("../tables/survival_analysis_input_files/addingCBC_chem_forVista_data5")

# 957 dogs with CBC/Chem longitudinal data
# 109 dogs have passed away
Dead = ldat %>% filter(death == "Dead")
nrow(Dead)
table(Dead$sex)
table(Dead$sterile)
summary(Dead$recruitmentAge)


##### simple cox model without any intermediate phenotypes as covariates (just effect of recruitment age, sex, sterilization, predicted height, best weight) #####
# scale the weight variables and recruitment age 
ldat$recruitmentAge <- scale(ldat$recruitmentAge)


multivar_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + Weight_Class_10KGBin_at_HLES + sex + sterile, ldat)

summary(multivar_coxmodel)$coefficients[, "Pr(>|z|)"]
summary(multivar_coxmodel)
multivar_ggforest <- ggforest(multivar_coxmodel)
ggsave("../vista_plots/cbc_chem_cox_multivariate_hazard_ratio.pdf", multivar_ggforest, width = 8, height = 8)


##### blood (CBC/Chem) cox proportional hazard model #####

# combine all the CBC and chem variables together 
ldat_add <- ldat %>% select(dog_id, recruitmentAge, sex, Weight_Class_10KGBin_at_HLES, tstart, tstop, death)
ldat_blood_scaled_df <- as.data.frame(scale(ldat[, VistasVars]))
ldat_blood_scaled_df <- cbind(ldat_add, ldat_blood_scaled_df)


# for each intermediate blood-based phenotype, test the effects of recruitment age, sex, sterilization, predicted height and best weight and the intermediate blood-based phenotype on survival status
# results are saved in a list of dataframes named by the blood-based phenotype (coefficients and summary of cox-ph model saved in that dataframe)
blood_univar_cox_results <- setNames(lapply(VistasVars, function(cov) {
  print(cov)
  formula <- as.formula(paste("Surv(tstart, tstop, death == 'Dead') ~ recruitmentAge + sex + Weight_Class_10KGBin_at_HLES +", cov))
  model <- coxph(formula, data = ldat_blood_scaled_df)
  as.data.frame(summary(model)$coefficients)
}), VistasVars)

# make a dataframe for blood phenotypes with log hazard coefficient, exponentiated coefficient and p-value
blood_coxph_df <- bind_rows(lapply(blood_univar_cox_results, function(df) df %>% select(coef, `exp(coef)`, `Pr(>|z|)`)))

# create a column containing  FDR/Benjamini-Hochberg adjusted p-values
blood_coxph_df$p_adj_fdr <- p.adjust(blood_coxph_df$`Pr(>|z|)`, method = "BH")

# create a column pertaining to covariates analyzed in survival analysis and how they related to life/death from the row names
blood_coxph_df$covariates <- rownames(blood_coxph_df) 

#set the row names to null
rownames(blood_coxph_df) <- NULL


blood_coxph_df <- blood_coxph_df %>% filter(covariates %in% VistasVars)

# filter dataframe to include p-values < 0.05 and further FDR adjusted p-values < 0.05
blood_coxph_sig_pval_df <- blood_coxph_df %>% filter(`Pr(>|z|)` < 0.05)
blood_coxph_sig_padj_df <- blood_coxph_df %>% filter(p_adj_fdr < 0.05)

# arrange by p_adj_fdr and only include intermediate phenotypes
blood_coxph_sig_padj_df %>% arrange(p_adj_fdr)

blood_coxph_sig_pval_df <- blood_coxph_sig_padj_df %>% filter(covariates %in% VistasVars) %>% arrange(p_adj_fdr)

bloodVars_included_mortality_ordered <- blood_coxph_df %>% arrange(abs(`Pr(>|z|)`)) %>% pull(covariates)

pdf("../vista_plots/plot_hazard_ratio_blood_phenotypes_scaled_weights_N-36_20251208.pdf")
lapply(bloodVars_included_mortality_ordered, function(cov) {
  formula <- as.formula(paste0("Surv(tstart, tstop, death == 'Dead') ~ recruitmentAge + sex + Weight_Class_10KGBin_at_HLES + `", cov, "`"))
  model <- coxph(formula, data = ldat_blood_scaled_df)
  ggforest(model)
})
dev.off()

# remove ldat, since the metabolite dataset R dataframe is also called ldat
cbc_chem_ldat <- ldat
rm(ldat)

##### metabolite cox proportional hazard model ######
##### read in Ben's metabolite RData file with longitudinal data ####
#file can be found in DataDryad within longitudinal_Rdata_files_for_survival_analysis
load("../tables/survival_analysis_input_files/longitudinalData2")
# clean up metabolite names in ldat and dogmzs
names(ldat) <- gsub('\\"', "", names(ldat))
names(ldat) <- gsub('/', "-", names(ldat))
names(ldat) <- gsub(',', "-", names(ldat))
names(ldat) <- gsub(' ', "-", names(ldat))
names(ldat) <- gsub(')', "", names(ldat))
names(ldat) <- gsub('\\(', "", names(ldat))

dogmzs <- gsub('\\"', "", dogmzs)
dogmzs <- gsub('/', "-", dogmzs)
dogmzs <- gsub(',', "-", dogmzs)
dogmzs <- gsub(' ', "-", dogmzs)
dogmzs <- gsub(')', "", dogmzs)
dogmzs <- gsub('\\(', "", dogmzs)

dogmzs_included <- dogmzs[dogmzs %in% intermed_phenos]

# scale metabolites
ldat_mzs_scaled_df <- as.data.frame(scale(ldat[, dogmzs_included]))
ldat_add <- ldat %>% select(dog_id, recruitmentAge, sex, bestWt, tstart, tstop, death)
ldat_mzs_scaled_df <- cbind(ldat_add, ldat_mzs_scaled_df)

# perform individual level survival analysis on 1 metabolite at a time 
univar_cox_metabolites_results <- setNames(lapply(dogmzs_included, function(cov){
  print(cov)
  formula <- as.formula(paste0("Surv(tstart, tstop, death == 'Dead') ~ recruitmentAge + sex + bestWt + `", cov, "`"))
  model <- coxph(formula, data = ldat_mzs_scaled_df)
  as.data.frame(summary(model)$coefficients)
}), dogmzs_included)

# make a dataframe for metabolites with log hazard coefficient, exponentiated coefficient and p-value
metabolites_coxph_df <- bind_rows(lapply(univar_cox_metabolites_results, function(df) df %>% select(coef, `exp(coef)`, `Pr(>|z|)`)))

metabolites_coxph_df$p_adj_fdr <- p.adjust(metabolites_coxph_df$`Pr(>|z|)`, method = "BH")

metabolites_coxph_df$covariates <- rownames(metabolites_coxph_df) 

rownames(metabolites_coxph_df) <- NULL

metabolites_coxph_df$covariates <- gsub("`", "", metabolites_coxph_df$covariates)

metabolites_coxph_df <- metabolites_coxph_df %>% filter(covariates %in% dogmzs_included)

metabolites_coxph_sig_pval_df <- metabolites_coxph_df %>% filter(`Pr(>|z|)` < 0.05)
metabolites_coxph_sig_padj_df <- metabolites_coxph_df %>% filter(p_adj_fdr < 0.05)

metabolites_coxph_sig_padj_df  %>% arrange(p_adj_fdr)


mzs_included_mortality_ordered <- metabolites_coxph_df %>% arrange(abs(`Pr(>|z|)`)) %>% pull(covariates)

pdf("../vista_plots/plot_hazard_ratio_metabolites_N-123_scaled_weights_recruitmentage_20251208.pdf")
lapply(mzs_included_mortality_ordered, function(cov) {
  print(cov)
  formula <- as.formula(paste0("Surv(tstart, tstop, death == 'Dead') ~ recruitmentAge + sex + bestWt + `", cov, "`"))
  model <- coxph(formula, data = ldat_mzs_scaled_df)
  ggforest(model)
})
dev.off()


# combine into a single dataframe metabolites_coxph_sig_padj_df and blood_coxph_df
endophenotype_coxph_df <- rbind(blood_coxph_df, metabolites_coxph_df)

# change cbc and chem phenotype to phenotype_name in phenotype column of intermediate_phenotype
endophenotype_coxph_df <- endophenotype_coxph_df %>% mutate(phenotype = case_when(
  str_detect(covariates, "_cbc_") ~ paste0(covariates, "_sqrt_transformed"),
  str_detect(covariates, "_cp_") ~ paste0(covariates, "_ln_transformed"),
  TRUE ~ covariates ))


#### potassium and globulins - high risk breeds composite variable in cox ph model ####
lmer_results <- read_sheet(dap_supp_url, sheet = "Data S8_LMER")

potassium_sig_results <- lmer_results %>% filter(phenotype == "krt_cp_potassium_value_ln_transformed") %>% arrange(desc(abs(REML.t.val))) %>% filter(ML.anova.p < 0.05)

# take top 9 breeds since potassium has 9 significant effects
globulin_sig_results <- lmer_results %>% filter(phenotype == "krt_cp_globulins_value_ln_transformed") %>% arrange(desc(abs(REML.t.val))) %>% filter(ML.anova.p < 0.05)

# albumin globulin ratio top breeds
alb_glob_sig_results <- lmer_results %>% filter(phenotype == "krt_cp_alb_glob_ratio_value_ln_transformed") %>% arrange(desc(abs(REML.t.val))) %>% filter(ML.anova.p < 0.05)

# sodium potassium ratio top breeds
sodium_potassium_sig_results <- lmer_results %>% filter(phenotype == "krt_cp_sp_ratio_value_ln_transformed") %>% arrange(desc(abs(REML.t.val))) %>% filter(ML.anova.p < 0.05)


potassium_top_breeds <- potassium_sig_results %>% pull(breed)
globulin_top_breeds <- globulin_sig_results %>% head(n=9) %>% pull(breed)
alb_glob_top_breeds <- alb_glob_sig_results %>% head(n=9) %>% pull(breed)
sodium_potassium_top_breeds <- sodium_potassium_sig_results %>% head(n = 9) %>% pull(breed)

ancestry <- read_sheet("https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/edit?gid=2137210148#gid=2137210148", sheet = "Data S15_BREED_ANC")
ancestry$dog_id <- as.character(ancestry$dog_id)

# ancestries added in 
anc_ldat_blood_scaled_df <- ldat_blood_scaled_df %>% left_join(ancestry, by = "dog_id")
anc_ldat_blood_scaled_df <- anc_ldat_blood_scaled_df %>% pivot_wider(names_from= "pop", values_from = "pct", values_fill = 0)

# create composite breed scores for the top 9 significant breeds in LMER analysis
anc_ldat_blood_scaled_df <- anc_ldat_blood_scaled_df %>% mutate(potassium_risk_breeds_composite_ancestry = rowSums(across(all_of(potassium_top_breeds))))
anc_ldat_blood_scaled_df <- anc_ldat_blood_scaled_df %>% mutate(globulin_risk_breeds_composite_ancestry = rowSums(across(all_of(globulin_top_breeds))))
anc_ldat_blood_scaled_df <- anc_ldat_blood_scaled_df %>% mutate(alb_glob_risk_breeds_composite_ancestry = rowSums(across(all_of(alb_glob_top_breeds))))
anc_ldat_blood_scaled_df <- anc_ldat_blood_scaled_df %>% mutate(sodium_potassium_risk_breeds_composite_ancestry = rowSums(across(all_of(sodium_potassium_top_breeds))))
summary(anc_ldat_blood_scaled_df$potassium_risk_breeds_composite_ancestry)

# run cox ph model with composite breed score as fixed effect 
potassium_anc_composite_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + Weight_Class_10KGBin_at_HLES + sex + krt_cp_potassium_value + potassium_risk_breeds_composite_ancestry, data = anc_ldat_blood_scaled_df)

summary(potassium_anc_composite_coxmodel)

potassium_tidy_cox <- as.data.frame(tidy(potassium_anc_composite_coxmodel, exponentiate = TRUE, conf.int = TRUE))

# forest plot with ggplot and the tidy dataframe above
potassium_anc_high_risk_breeds_plot <- ggplot(potassium_tidy_cox, aes(x = estimate, y = term)) +
  geom_point(size = 3) +                                   # points for estimates
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),   # horizontal error bars
                 height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +        # reference line at 1
  scale_x_continuous(name = "hazard ratio") +
  geom_text(aes(label = ifelse(p.value < 0.001, "<0.001", sprintf("p=%.3f", p.value))),
            hjust = -2.6, vjust = -1, size = 4)+
  geom_text(aes(label = round(estimate,2)), vjust = -1, size = 4) +
  ylab("") +
  theme_cowplot()

ggsave("../vista_plots/potassium_forestplot_w_ancestry_high_risk_breeds_20251208.pdf", potassium_anc_high_risk_breeds_plot, width = 12, height = 6)

# globulin
globulin_anc_composite_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + Weight_Class_10KGBin_at_HLES + sex + krt_cp_globulins_value + globulin_risk_breeds_composite_ancestry, data = anc_ldat_blood_scaled_df)

summary(globulin_anc_composite_coxmodel)

globulin_tidy_cox <- as.data.frame(tidy(globulin_anc_composite_coxmodel, exponentiate = TRUE, conf.int = TRUE))

# forest plot with ggplot and the tidy dataframe above
globulin_anc_high_risk_breeds_plot <- ggplot(globulin_tidy_cox, aes(x = estimate, y = term)) +
  geom_point(size = 3) +                                   # points for estimates
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),   # horizontal error bars
                 height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +        # reference line at 1
  scale_x_continuous(name = "hazard ratio") +
  geom_text(aes(label = ifelse(p.value < 0.001, "<0.001", sprintf("p=%.3f", p.value))),
            hjust = -2.6, vjust = -1, size = 4)+
  geom_text(aes(label = round(estimate,2)), vjust = -1, size = 4) +
  ylab("") +
  theme_cowplot()

ggsave("../vista_plots/globulin_forestplot_w_ancestry_high_risk_breeds_20251208.pdf", globulin_anc_high_risk_breeds_plot, width = 12, height = 6)


# albumin globulin ratio
# run cox ph model with composite breed score as fixed effect 
alb_glob_composite_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + Weight_Class_10KGBin_at_HLES + sex + krt_cp_alb_glob_ratio_value + alb_glob_risk_breeds_composite_ancestry, data = anc_ldat_blood_scaled_df)

summary(alb_glob_composite_coxmodel)

alb_glob_tidy_cox <- as.data.frame(tidy(alb_glob_composite_coxmodel, exponentiate = TRUE, conf.int = TRUE))


# sodium potassium ratio 
sodium_potassium_composite_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + Weight_Class_10KGBin_at_HLES + sex + krt_cp_sp_ratio_value + sodium_potassium_risk_breeds_composite_ancestry, data = anc_ldat_blood_scaled_df)

summary(sodium_potassium_composite_coxmodel)

sodium_potassium_tidy_cox <- as.data.frame(tidy(sodium_potassium_composite_coxmodel, exponentiate = TRUE, conf.int = TRUE))

# composite ancestry for hydroxyproline and effect of top 9 risk breeds impacting hydroxyproline
# ancestries added in 
anc_ldat_mzs_scaled_df <- ldat_mzs_scaled_df %>% left_join(ancestry, by = "dog_id")
anc_ldat_mzs_scaled_df <- anc_ldat_mzs_scaled_df %>% pivot_wider(names_from= "pop", values_from = "pct", values_fill = 0)

hydroxyp_sig_results <- lmer_results %>% filter(phenotype == "Hydroxyproline") %>% arrange(desc(abs(REML.t.val))) %>% filter(ML.anova.p < 0.05)
hydroxyp_top_breeds <- hydroxyp_sig_results %>% pull(breed)

# create composite breed scores for the top 9 significant breeds in LMER analysis
anc_ldat_mzs_scaled_df <- anc_ldat_mzs_scaled_df %>% mutate(hydroxyproline_risk_breeds_composite_ancestry = rowSums(across(all_of(hydroxyp_top_breeds))))

# hydroxyproline multivariate cox model with composite breed ancestry as predictor
hydroxyproline_composite_coxmodel <- coxph(Surv(tstart, tstop, death == "Dead") ~ recruitmentAge + bestWt + sex + Hydroxyproline + hydroxyproline_risk_breeds_composite_ancestry, data = anc_ldat_mzs_scaled_df)

summary(hydroxyproline_composite_coxmodel)

hydroxyp_tidy_cox <- as.data.frame(tidy(hydroxyproline_composite_coxmodel, exponentiate = TRUE, conf.int = TRUE))

# combine the 5 traits with composite breed ancestry results 
alb_glob_tidy_cox$phenotype = "albumin_globulin_ratio"
globulin_tidy_cox$phenotype = "globulin"
potassium_tidy_cox$phenotype = "potassium"
sodium_potassium_tidy_cox$phenotype = "sodium_potassium_ratio"
hydroxyp_tidy_cox$phenotype = "Hydroxyproline"

five_aging_based_traits_tidy_cox = rbind(alb_glob_tidy_cox, globulin_tidy_cox, potassium_tidy_cox, sodium_potassium_tidy_cox, hydroxyp_tidy_cox)

five_aging_based_traits_tidy_cox$term <- gsub("bestWt", "best_weight", five_aging_based_traits_tidy_cox$term)
write_sheet(five_aging_based_traits_tidy_cox, ss = dap_supp_url, sheet = "Data S_COX_HAZ_W_BREED_ANC")

