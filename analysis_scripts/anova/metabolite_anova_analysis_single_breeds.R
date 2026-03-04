# Goal: Conducting ANOVA analysis on metabolites in single breed dogs 
#       Dependent variable: metabolite values
#       Independent variable: age, weight, sex status, and breed
#       Testing the effects of age, weight, and breed (for single breed dogs) on metabolite values

library(tidyverse)
library(rstatix)


### input files ###
# load metabolite data
# these files can be found in DataDryad within dap_gwas_input_files directory
metabolites_120 <- read_tsv("../blood_data/blood_phenotype_datasets/precision2024_120metabolites_937dogs_phenotypes.tsv")
metabolites_3 <- read_tsv("../blood_data/blood_phenotype_datasets/precision2024_3metabolites_937dogs_phenotypes.tsv")
metabolites_all <- metabolites_120 %>% left_join(metabolites_3, by = "dog_id")


# more efficient to read from supplemental table
if (!exists("dat")) {
  dat <- readRDS("DAP_supp_data.rds") # supplemental material tables and datafiles located in DataDryad
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_S15_BREED_ANC, etc.
}


# read in age, weight from quantitative covariate file from GWAS found in DataDryad within dap_gwas_input_files directory
age_weight_df <- read_tsv("../blood_data/metabolite_gwas_input/age_weight_hoursFasted_937_precision_dogs.tsv")
age_weight_df <- age_weight_df %>% select(-c("Sample_Dog_Hours_Fasted", "FID"))

sex_status_df <- read_tsv("../blood_data/metabolite_gwas_input/Sex_Class_at_HLES_937_precision_dogs.tsv")
sex_status_df <- sex_status_df %>% select(-c("FID"))
covar_df <- age_weight_df %>% left_join(sex_status_df)

# read in kurtosis file to obtain metabolite results that will be included
# set threshold to including phenotypes with kurtosis < 10
# kurtosis file found in DataDryad located in baseline_phenotype_kurtosis_survey_response_distribution_files
metabolite_kurtosis_df <- read_tsv("../blood_data/kurtosis/metabolites_kurtosis.tsv")

metabolites_included <- metabolite_kurtosis_df %>% filter(kurtosis_result < 10) %>% pull(pheno)
metabolites_all  <- metabolites_all %>% select(dog_id, all_of(metabolites_included))


metabolite_df <- metabolites_all %>% left_join(covar_df)

# replace the ',' in single variable called Sex_Class_at_HLES with "_" and rename to sex_status
metabolite_df$sex_status <- gsub(", ", "_", metabolite_df$Sex_Class_at_HLES)
metabolite_df <- metabolite_df %>% select(-c(Sex_Class_at_HLES))

names(metabolite_df)[which(names(metabolite_df) == "age_yrs_at_sample_collection")] <- "age" 
names(metabolite_df)[which(names(metabolite_df) == "weight_kg")] <- "weight" 



# convert metabolite_df from wide to long format
metabolite_long_df <- metabolite_df %>% pivot_longer(cols = !c("dog_id", "age", "weight", "sex_status"), names_to = "metabolite_name", values_to = "metabolite_value")


#### for dogs that are single breed given ancestry, see the effect of breed on metabolite levels in addition to age, weight, sex, and sterilization status#####

# single breed dogs with metabolite data
singlebreed_cohort <- raw_S14_SEQ_DOGS %>% filter(single_breed == TRUE)

singlebreed_metabolite_cohort <- singlebreed_cohort %>% filter(dog_id %in% metabolite_long_df$dog_id)

# to get a sense for number of dogs in each breed category (use standardized breed that has standardized breed names for single-breed dogs)
breed_count_anova <- singlebreed_metabolite_cohort %>% group_by(standardized_breed) %>% count()
breed_count_anova %>% filter(n > 1) %>% arrange(desc(n)) %>% tail()

# we retain 81 breeds (more than 1 dog in each breed)
breed_count_anova_kept <- breed_count_anova %>% filter(n > 1) 
nrow(breed_count_anova_kept)

singlebreed_metabolite_cohort <- singlebreed_metabolite_cohort %>% select(c(dog_id, standardized_breed))

singlebreed_metabolite_long_df <- singlebreed_metabolite_cohort %>% inner_join(metabolite_long_df, by = "dog_id")

sum(is.na(singlebreed_metabolite_long_df))


# keep breeds that are in breed_count_anova_kept (having more than 1 dog)
singlebreed_metabolite_long_df <- singlebreed_metabolite_long_df %>% filter(standardized_breed %in% breed_count_anova_kept$standardized_breed)

# only retain metabolites that will be included in genetic analysis (kurtosis < 10)
singlebreed_metabolite_long_df <- singlebreed_metabolite_long_df %>% filter(metabolite_name %in% metabolites_included)


# number of single breed dogs in ANOVA analysis (474)
num_dogs_anova <- length(unique(singlebreed_metabolite_long_df$dog_id))

sample_size <- singlebreed_metabolite_long_df %>% group_by(metabolite_name) %>% summarize(num_dogs = n_distinct(dog_id))
# number of single breeds represented in ANOVA analysis (81)
num_breeds_represented = length(unique(singlebreed_metabolite_long_df$standardized_breed))

singlebreed_metabolite_df <- metabolite_df %>% inner_join(singlebreed_metabolite_cohort)


# trying to figure out effect of breed after including weight in model
# variance inflation factor (VIF) of weight when breed is included in model for creatinine
mod1 <- lm(Creatinine ~ age + weight, data = singlebreed_metabolite_df)


# Model 2: age + weight + breed
mod2 <- lm(Creatinine ~ age + weight + standardized_breed, data = singlebreed_metabolite_df)

anova(mod1, mod2)  # Tests whether adding breed significantly improves model fit

library(car)
Anova(mod2, type = 2)  # Type II SS

vif(mod2)

# performing ANOVA (including breed as independent variable as well)
names(singlebreed_metabolite_long_df)[which(names(singlebreed_metabolite_long_df) == "standardized_breed")] <- "Breed"

anova_metabolite_breed_summary <- as_tibble(singlebreed_metabolite_long_df %>% group_by(metabolite_name) %>% anova_test(metabolite_value ~ age + Breed + weight + sex_status))

# replace the word "Breed" for population with "breed" to make result easier to analyze
anova_metabolite_breed_summary$Effect <- gsub("Breed", "breed", anova_metabolite_breed_summary$Effect)

anova_metabolite_breed_summary$model <- "age + breed + weight + sex status"
anova_metabolite_breed_summary$Effect <- gsub("sex_status", "sex status", anova_metabolite_breed_summary$Effect)

# add number of breeds and number of dogs in analysis to dataframe and type of blood-based phenotype
anova_metabolite_breed_summary$num_dogs <- num_dogs_anova
anova_metabolite_breed_summary$num_breeds <- num_breeds_represented
anova_metabolite_breed_summary$paper_phenotype_category <- "Plasma metabolites"

#anova_metabolite_breed_summary_arranged <- anova_metabolite_breed_summary %>% arrange(desc(ges))

write_tsv(anova_metabolite_breed_summary, paste0("./anova_results/metabolite_singlebreed_cohort_N-", num_dogs_anova, "_", num_breeds_represented, "breeds_anova-age_breed_weight_sexStatus.tsv"), col_names = TRUE)


#### single breed ANOVA excluding breed effect on 456 dogs ####
anova_metabolite_singlebreed_nobreed_summary <- as_tibble(singlebreed_metabolite_long_df %>% group_by(metabolite_name) %>% anova_test(metabolite_value ~ age + weight + sex_status))

anova_metabolite_singlebreed_nobreed_summary$model <- "age + weight + sex status"


anova_metabolite_singlebreed_nobreed_summary$Effect <- gsub("sex_status", "sex status", anova_metabolite_singlebreed_nobreed_summary$Effect)

# add number of breeds and number of dogs in analysis to dataframe and type of blood-based phenotype
anova_metabolite_singlebreed_nobreed_summary$num_dogs <- num_dogs_anova
anova_metabolite_singlebreed_nobreed_summary$num_breeds <- num_breeds_represented
anova_metabolite_singlebreed_nobreed_summary$paper_phenotype_category <- "Plasma metabolites"

anova_metabolite_singlebreed_nobreed_summary_arranged <- anova_metabolite_singlebreed_nobreed_summary %>% arrange(desc(ges))

write_tsv(anova_metabolite_singlebreed_nobreed_summary, paste0("./anova_results/metabolite_singlebreed_cohort_N-", num_dogs_anova, "_", num_breeds_represented, "_anova-age_weight_sexStatus.tsv"), col_names = TRUE)


### summarize analysis for paper ###
breed_effects <- anova_metabolite_breed_summary %>% filter(Effect == "breed")
significant_breed_effects <- breed_effects %>% filter(p < 0.05)
nrow(significant_breed_effects)
mean(significant_breed_effects$ges)
sd(significant_breed_effects$ges)

age_effects_with_breed <- anova_metabolite_breed_summary %>% filter(Effect == "age")
significant_age_effects_with_breed <- age_effects_with_breed %>% filter(p < 0.05)
nrow(significant_age_effects_with_breed)
mean(significant_age_effects_with_breed$ges)
sd(significant_age_effects_with_breed$ges)

