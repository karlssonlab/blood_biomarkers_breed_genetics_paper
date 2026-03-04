# Goal: Conducting ANOVA analysis on complete blood count values 
#       Dependent variable: blood count values
#       Independent variable: age, weight, sex status and breed
#       Testing the effects of age, weight, sex status and breed (for single breed dogs) on complete blood count values


library(tidyverse)
library(rstatix)



### input files ###
# read in blood chemistry phenotype file 
# these files can be found in DataDryad within dap_gwas_input_files directory (specifically in blood_cbc_chem subdirectory)
cbc <- read_tsv("../blood_data/blood_phenotype_datasets/precision2024_20cbc_phenotypes.tsv") %>% select(-c(FID))

cbc_qcovars <- read_tsv("../blood_data/cbc_chem_gwas_input/age_at_sample_collection_weight_kg_hours_fasted_cbc_935_precision_dogs.tsv", col_names = TRUE) %>% select(-c(FID, Sample_Dog_Hours_Fasted))
cbc_dcovars <- read_tsv("../blood_data/cbc_chem_gwas_input/Sex_Class_at_HLES_972_precision_dogs.tsv", col_names = TRUE) %>% select(-c(FID))


# I am using sex_class_at_HLES from DogOverview table that is used as GWAS covariate for blood chemistry and CBC

cbc_covars <- cbc_qcovars %>% inner_join(cbc_dcovars, by = "dog_id")

cbc_df <- cbc_covars %>% inner_join(cbc, by = "dog_id")

# read in sequenced dog demography containing information on single breed dog status
if (!exists("dat")) {
  dat <- readRDS("DAP_supp_data.rds") # supplemental material tables and datafiles located in DataDryad
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_S15_BREED_ANC, etc.
}


# select dog_id, age, breed, weight, sex and sterilization status with all cbc results
names(cbc_df)[names(cbc_df) == "age_yrs_at_sample_collection"] <- "age"
names(cbc_df)[names(cbc_df) == "weight_kg"] <- "weight"
names(cbc_df)[names(cbc_df) == "Sex_Class_at_HLES"] <- "sex_class_at_HLES"

singlebreed_cohort <- raw_S14_SEQ_DOGS %>% filter(single_breed == TRUE)
singlebreed_cbc_cohort <- singlebreed_cohort %>% filter(dog_id %in% cbc_df$dog_id)
singlebreed_cbc_cohort <- singlebreed_cbc_cohort %>% select(c(dog_id, standardized_breed))



breed_count_anova <- singlebreed_cbc_cohort %>% group_by(standardized_breed) %>% count()



breed_count_anova_kept <- breed_count_anova %>% filter(n > 1)


cbc_df <- cbc_df %>% select(c(dog_id, age, weight, sex_class_at_HLES, any_of(names(cbc))))

cbc_df <- cbc_df %>% inner_join(singlebreed_cbc_cohort, by = "dog_id")

cbc_df <- cbc_df %>% filter(standardized_breed %in% breed_count_anova_kept$standardized_breed)


# get number of NAs per column
cbc_na_counts <- cbc_df %>%
  summarise_all(~ sum(is.na(.)))


# split platelet from the rest since it has the most NAs and perform ANOVA on this dataset separately
cbc_platelet_df <- cbc_df %>% select(dog_id, age, standardized_breed, weight, sex_class_at_HLES, krt_cbc_plt_sqrt_transformed)

# all other cbc columns present in this dataframe for anova
cbc_noplt_df <- cbc_df %>% select(-c(krt_cbc_plt_sqrt_transformed))

# convert cbc_df from wide to long format
cbc_long_df <- cbc_noplt_df %>% pivot_longer(cols = !c("dog_id", "age", "standardized_breed", "weight", "sex_class_at_HLES"), names_to = "cbc_name", values_to = "cbc_value")


# drop any NAs excluding ones from platelet since cbc_long_df will run ANOVA on all other blood count phenotypes except for platelets
cbc_long_df <- cbc_long_df %>% na.omit()


# 468 dogs max 
length(unique(cbc_long_df$dog_id))

# 83 breeds
num_breeds_represented = length(unique(cbc_long_df$standardized_breed))

anova_cbc_notplt_sample_sizes <- cbc_long_df %>% group_by(cbc_name) %>% summarize(num_dogs = n_distinct(dog_id))


anova_cbc_notplt_summary <- as_tibble(cbc_long_df %>% group_by(cbc_name) %>% anova_test(cbc_value ~ age + standardized_breed + weight + sex_class_at_HLES))

anova_cbc_notplt_summary$model <- "age + breed + weight + sex_class_at_HLES"

anova_cbc_notplt_summary <- anova_cbc_notplt_summary %>% inner_join(anova_cbc_notplt_sample_sizes, by = "cbc_name")

anova_cbc_notplt_summary$num_breed <- num_breeds_represented

# include full blood phenotype names in addition to abbreviated phenotype names
cbc_key <- raw_S2_ALL_PHENOS %>% select(phenotype, paper_phenotype_category)
anova_cbc_notplt_summary <- anova_cbc_notplt_summary %>% left_join(cbc_key, by = c("cbc_name" = "phenotype"))

anova_cbc_notplt_summary_ges_sort <- anova_cbc_notplt_summary %>% arrange(desc(ges))

# platelet anova 
# remove any NAs in cbc_platelet dataframe
cbc_platelet_df <- cbc_platelet_df %>% na.omit()

# 81 breeds
plt_num_breeds = length(unique(cbc_platelet_df$standardized_breed))

# run ANOVA on platelet phenotype
anova_cbc_plt_summary <- cbc_platelet_df %>% anova_test(krt_cbc_plt_sqrt_transformed ~ age + standardized_breed + weight + sex_class_at_HLES)

# define model covariates
anova_cbc_plt_summary$model <- "age + breed + weight + sex_class_at_HLES"

# add number of dogs for cbc_platelet anova
anova_cbc_plt_summary$num_dogs <- c(length(unique(cbc_platelet_df$dog_id)))
anova_cbc_plt_summary$num_breed <- plt_num_breeds
# include full blood phenotype names in addition to abbreviated phenotype names
anova_cbc_plt_summary$cbc_name <- "krt_cbc_plt_sqrt_transformed"
anova_cbc_plt_summary <- anova_cbc_plt_summary %>% left_join(cbc_key, by = c("cbc_name" = "phenotype"))




anova_cbc_summary <- rbind(anova_cbc_notplt_summary, anova_cbc_plt_summary)
# replace the word "standardized_breed" for population with "breed" to make result easier to analyze
anova_cbc_summary$Effect <- gsub("standardized_breed", "breed", anova_cbc_summary$Effect)
anova_cbc_summary$Effect <- gsub("sex_class_at_HLES", "sex status", anova_cbc_summary$Effect)


anova_cbc_summary$model<- gsub("sex_class_at_HLES", "sex status", anova_cbc_summary$model)



write_tsv(anova_cbc_summary, "./anova_results/cbc_singlebreed_cohort_anova-age_breed_weight_sexClassatHLES.tsv", col_names = TRUE)


#### excluding breed as predictor ####
anova_cbc_nobreed_notplt_summary <- as_tibble(cbc_long_df %>% group_by(cbc_name) %>% anova_test(cbc_value ~ age + weight + sex_class_at_HLES))


anova_cbc_nobreed_notplt_summary$model <- "age + weight + sex_class_at_HLES"

anova_cbc_nobreed_notplt_summary <- anova_cbc_nobreed_notplt_summary %>% inner_join(anova_cbc_notplt_sample_sizes, by = "cbc_name")

anova_cbc_nobreed_notplt_summary$num_breed <- num_breeds_represented
# include full blood phenotype names in addition to abbreviated phenotype names
anova_cbc_nobreed_notplt_summary <- anova_cbc_nobreed_notplt_summary %>% left_join(cbc_key, by = c("cbc_name" = "phenotype"))

# for platelets

anova_cbc_nobreed_plt_summary <- cbc_platelet_df %>% anova_test(krt_cbc_plt_sqrt_transformed ~ age + weight + sex_class_at_HLES)

# define model covariates
anova_cbc_nobreed_plt_summary$model <- "age + weight + sex_class_at_HLES"

# add number of dogs for cbc_platelet anova
anova_cbc_nobreed_plt_summary$num_dogs <- c(length(unique(cbc_platelet_df$dog_id)))
anova_cbc_nobreed_plt_summary$num_breed <- plt_num_breeds
# include full blood phenotype names in addition to abbreviated phenotype names
anova_cbc_nobreed_plt_summary$cbc_name <- "krt_cbc_plt_sqrt_transformed"
anova_cbc_nobreed_plt_summary <- anova_cbc_nobreed_plt_summary %>% left_join(cbc_key, by = c("cbc_name" = "phenotype"))




anova_cbc_nobreed_summary <- rbind(anova_cbc_nobreed_notplt_summary, anova_cbc_nobreed_plt_summary)


anova_cbc_nobreed_summary$Effect <- gsub("sex_class_at_HLES", "sex status", anova_cbc_nobreed_summary$Effect)
anova_cbc_nobreed_summary$model <- gsub("sex_class_at_HLES", "sex status", anova_cbc_nobreed_summary$model)

write_tsv(anova_cbc_nobreed_summary, "./anova_results/cbc_singlebreed_cohort_anova-age_weight_sexClassatHLES.tsv", col_names = TRUE)


