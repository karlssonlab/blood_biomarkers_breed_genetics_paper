# Goal: Conducting ANOVA analysis on blood chemistry values 
#       Dependent variable: blood chemistry values
#       Independent variable: age, weight, and breed
#       Testing the effects of age, weight,  and breed (for single breed dogs) on blood chemistry values


library(tidyverse)
library(rstatix)
library(googlesheets4)

options(scipen = 999)

### input files ###
# read in blood chemistry phenotype file
# these files can be found in DataDryad within dap_gwas_input_files directory (specifically in blood_cbc_chem subdirectory)
blood_chemistry <- read_tsv("../blood_data/blood_phenotype_datasets/precision2024_bloodchem_16phenotypes.tsv") %>% select(-c(FID))

blood_chemistry_qcovars <- read_tsv("../blood_data/cbc_chem_gwas_input/age_at_sample_collection_weight_kg_hours_fasted_chem_966_precision_dogs.tsv", col_names = TRUE) %>% select(-c(FID, Sample_Dog_Hours_Fasted))
blood_chemistry_dcovars <- read_tsv("../blood_data/cbc_chem_gwas_input/Sex_Class_at_HLES_972_precision_dogs.tsv", col_names = TRUE) %>% select(-c(FID))


# for now, I am using sex_class_at_HLES from DogOverview table that is used as GWAS covariate for blood chemistry and CBC

blood_chemistry_covars <- blood_chemistry_qcovars %>% inner_join(blood_chemistry_dcovars, by = "dog_id")

blood_chem_df <- blood_chemistry_covars %>% inner_join(blood_chemistry, by = "dog_id")
names(blood_chem_df)[names(blood_chem_df) == "age_yrs_at_sample_collection"] <- "age"
names(blood_chem_df)[names(blood_chem_df) == "weight_kg"] <- "weight"
names(blood_chem_df)[names(blood_chem_df) == "Sex_Class_at_HLES"] <- "sex_class_at_HLES"


# read in sequenced dog demography containing information on single breed dog status 
if (!exists("dat")) {
  dat <- readRDS("DAP_supp_data.rds") # supplemental material tables and datafiles located in DataDryad
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_S15_BREED_ANC, etc.
}



sum(is.na(blood_chem_df))
# remove any dogs that have NAs for a particular blood chemistry value
blood_chem_df <- na.omit(blood_chem_df)

# single breed dogs with metabolite data
singlebreed_cohort <- raw_S14_SEQ_DOGS %>% filter(single_breed == TRUE)



# convert blood_chem_df from wide to long format
#blood_chem_long_df <- blood_chem_df %>% pivot_longer(cols = !c("dog_id", "age", "weight", "sex_class_at_HLES", "genomic_predicted_height"), names_to = "blood_chem_name", values_to = "blood_chem_value")
blood_chem_long_df <- blood_chem_df %>% pivot_longer(cols = !c("dog_id", "age", "weight", "sex_class_at_HLES"), names_to = "blood_chem_name", values_to = "blood_chem_value")

singlebreed_bchem_cohort <- singlebreed_cohort %>% filter(dog_id %in% blood_chem_long_df$dog_id)
singlebreed_bchem_cohort <- singlebreed_bchem_cohort %>% select(c(dog_id, standardized_breed))



breed_count_anova <- singlebreed_bchem_cohort %>% group_by(standardized_breed) %>% count()

breed_count_anova_kept <- breed_count_anova %>% filter(n > 1)

singlebreed_blood_chem_long_df <- blood_chem_long_df %>% filter(dog_id %in% singlebreed_bchem_cohort$dog_id)

# get breed information
singlebreed_blood_chem_long_df <- singlebreed_blood_chem_long_df %>% inner_join(singlebreed_bchem_cohort, by = "dog_id")


# retain breeds that have greater than 1 dog represented
singlebreed_blood_chem_long_df <- singlebreed_blood_chem_long_df %>% filter(standardized_breed %in% breed_count_anova_kept$standardized_breed)

# single breed blood chem with more than 1 dog in each breed represented
lm_df <- blood_chem_df %>% inner_join(singlebreed_bchem_cohort, by = "dog_id") %>% filter(standardized_breed %in% breed_count_anova_kept$standardized_breed)
age_weight_mod <- lm(krt_cp_creatinine_value_ln_transformed ~ age + weight, data = lm_df )
breed_age_weight_mod <- lm(krt_cp_creatinine_value_ln_transformed ~ age + weight + standardized_breed, data = lm_df)
anova(age_weight_mod, breed_age_weight_mod)

library(car)
Anova(breed_age_weight_mod, type = 2)  # Type II SS
# 542 single breed dogs 
length(unique(singlebreed_blood_chem_long_df$dog_id))

# sample sizes for anova (475 single breed dogs)
anova_bloodchem_sample_sizes <- singlebreed_blood_chem_long_df %>% group_by(blood_chem_name) %>% summarize(num_dogs = n_distinct(dog_id))
num_dogs <- unique(anova_bloodchem_sample_sizes$num_dogs)
# number of breeds in analysis
num_breeds_represented <- length(unique(singlebreed_blood_chem_long_df$standardized_breed))

anova_breed_bloodchem_summary <- as_tibble(singlebreed_blood_chem_long_df %>% group_by(blood_chem_name) %>% anova_test(blood_chem_value ~ age + Breed + weight + sex_class_at_HLES))



# covariates used in anova
anova_breed_bloodchem_summary$model <- "age + breed + weight + sex_class_at_HLES"

# replace the word "Breed" for population with "breed" to make result easier to analyze
anova_breed_bloodchem_summary$Effect <- gsub("standardized_breed", "breed", anova_breed_bloodchem_summary$Effect)
anova_breed_bloodchem_summary$model <- gsub("sex_class_at_HLES", "sex status", anova_breed_bloodchem_summary$model)
anova_breed_bloodchem_summary$Effect <- gsub("sex_class_at_HLES", "sex status", anova_breed_bloodchem_summary$Effect)

# include full blood phenotype names in addition to abbreviated phenotype names (not for supplemental table version)
bchem_key <- raw_S2_ALL_PHENOS %>% select(phenotype, paper_phenotype_category)
anova_breed_bloodchem_summary <- anova_breed_bloodchem_summary%>% left_join(bchem_key, by = c("blood_chem_name" = "phenotype"))

# add number of dogs and breeds analyzed and the phenotype category for inclusion in supplemental table
anova_breed_bloodchem_summary$num_dogs <- num_dogs
anova_breed_bloodchem_summary$num_breeds <- num_breeds_represented

anova_breed_bloodchem_summary$paper_phenotype_category <- "Clinical analytes"
anova_breed_bloodchem_summary_ges_sort <- anova_breed_bloodchem_summary %>% arrange(desc(ges))


write_tsv(anova_breed_bloodchem_summary, paste0("./anova_results/updated_precision_set_2025/blood_chemistry_singlebreed_cohort_N-", length(unique(singlebreed_blood_chem_long_df$dog_id)), "_", num_breeds_represented, "_anova-age_breed_weight_sexClassatHLES.tsv"), col_names = TRUE)


#### for dogs that are single breed given ancestry file, see the effect of absence of breed on blood chemistry (covariates: age, weight, sex, and sterilization status) #####


anova_bloodchem_nobreed_summary <- as_tibble(singlebreed_blood_chem_long_df  %>% group_by(blood_chem_name) %>% anova_test(blood_chem_value ~ age + weight + sex_class_at_HLES))


anova_bloodchem_nobreed_summary$model <- "age + weight + sex_class_at_HLES"

# replace terms with consistent labels used throughout metabolite and CBC ANOVA data tables
anova_bloodchem_nobreed_summary$model <- gsub("sex_class_at_HLES", "sex status", anova_bloodchem_nobreed_summary$model)
anova_bloodchem_nobreed_summary$Effect <- gsub("sex_class_at_HLES", "sex status", anova_bloodchem_nobreed_summary$Effect)




anova_bloodchem_nobreed_summary_arranged <- anova_bloodchem_nobreed_summary %>% arrange(desc(ges))

# include full blood phenotype names in addition to abbreviated phenotype names (not for supplemental table version)
anova_bloodchem_nobreed_summary <- anova_bloodchem_nobreed_summary %>% left_join(bchem_key, by = c("blood_chem_name" = "phenotype"))

# add number of dogs and breeds analyzed and the phenotype category for inclusion in supplemental table
anova_bloodchem_nobreed_summary$num_dogs <- num_dogs
anova_bloodchem_nobreed_summary$num_breeds <- num_breeds_represented

anova_bloodchem_nobreed_summary$paper_phenotype_category <- "Clinical analytes"


write_tsv(anova_bloodchem_nobreed_summary, paste0("./anova_results/updated_precision_set_2025/blood_chemistry_singlebreed_cohort_N-", length(unique(singlebreed_blood_chem_long_df$dog_id)), "_", num_breeds_represented, "_anova-age_weight_sexClassatHLES.tsv"), col_names = TRUE)

