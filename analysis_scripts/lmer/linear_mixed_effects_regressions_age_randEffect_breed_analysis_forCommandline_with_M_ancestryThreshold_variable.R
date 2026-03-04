# Linear mixed effects regression models for any phenotype (ie Darwin's Ark OCD traits here)

# Conditions:
# 1. breed ancestry up to 100% from a breed
# 2. only dogs with <100% ancestry from any one breed
# 3. only breeds with a SD >5% ancestry across dogs carrying that breed

# Variables:
# independent variable: phenotype of interest
# fixed effects: ancestry from each breed
# random effects: covariance matrix of genetic relatedness, age bracket (3 equal age brackets - tertiles)

# Method:
# 1. build models with restricted maximum likelihood (REML) to obtain unbiased estimates, standard deviations, and Wald statistics (t.val) for the fixed effects of breed on factor scores
# 2. perform analysis-of-variance (ANOVA) on each factor model to obtain the breed F-statistics
# 3. contribute models with maximum likelihood (ML) - breed.i
# 4. perform ANOVA to obtain likelihood ratio for each breed and report p-values for likelihood tests between the model with ML (all breeds) and models with ML -breed.i
# NOTE! We report p-values for the likelihood tests between the ML models but not for REML models. Why? See: https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html
# 5. obtain the proportion of variance explained by the fixed effects (breed ancestries) as the marginal Nakagawa's R2 for each factor modeled with and without restricted maximum likelihood
# 6. obtain the conditional R2, which is the variance explained by both the random effects, which are the age bracket and kinship covariance matrix, and fixed effects, which are the breed ancestries
# NOTE! We could not always ascertain the conditional R2 due to singularity in some models from the random effects structure being too complex.

# submit by phenotype (each phenotype is its own job)

.libPaths(c("/home/vsohrab/R/x86_64-pc-linux-gnu-library/4.4",.libPaths()))

# Load Libraries
# plotting and visualization
require(ggplot2)
require(ggpubr)
require(tidyverse)

# data management
require(tidyr)
require(dplyr)
require(plyr)
require(readr)
require(data.table)
require(reshape2)
require(zoo)

# modeling
require(devtools)
require(lme4)
# remotes::install_github("palday/coefplot2", subdir = "pkg")
require(coefplot2)
# devtools::install_github("variani/lme4qtl")
require(lme4qtl)
require(psych)
require(insight)
require(performance)

options(scipen = 999)

# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

#### input files ####
args <- commandArgs(trailingOnly=T)
grm <- args[1]
ancestry <- args[2]
pheno_data <- args[3]
random_effects_data <- args[4]
outdir <- args[5]
outname <- args[6]
pheno_outname <- args[7]
ancestry_cutoff <- args[8]

setwd(outdir)

#### Load Data ####
GRM = ReadGRMBin(grm)

breedcalls =  read_tsv(ancestry) %>% as.data.table()

names(breedcalls) <- c("dog", "breed", "pct")

breedcalls$breed <- tolower(breedcalls$breed)


phenotypes <- read_tsv(pheno_data)
colnames(phenotypes)[which(colnames(phenotypes) == "IID")] <- "dog_id"
randeffects <- read_tsv(random_effects_data, col_names = TRUE)
names(randeffects) <- c("FID", "dog_id", "age")

phenotypes$dog_id <- as.numeric(phenotypes$dog_id)
randeffects$dog_id <- as.numeric(randeffects$dog_id)

phenotypes$FID <- as.numeric(phenotypes$FID)
randeffects$FID <- as.numeric(randeffects$FID)

head(phenotypes)

phenotypes <- phenotypes %>% select(dog_id, FID, pheno_outname)
head(phenotypes)

head(randeffects)
pheno_df <- randeffects %>% inner_join(phenotypes, by = c("dog_id", "FID")) %>% select(-c("FID"))

head(pheno_df)

pheno_long_df <- pheno_df %>% pivot_longer(cols = !c("dog_id", "age"), names_to = "pheno_name", values_to = "pheno_value")


#### Mixed Linear Effects Models with Covariance Matrix ####
# get the covariance matrix from genetic relationship matrix (GRM)
GRMCorr <- matrix(NA, ncol = length(GRM$diag), nrow = length(GRM$diag))
GRMCorr[lower.tri(GRMCorr)] <- GRM$off
GRMCorr[upper.tri(GRMCorr)] <- t(GRMCorr)[upper.tri(t(GRMCorr))]
diag(GRMCorr) <- 1 # GRM not exact 1 on diagonal due to some nonsense I spoke with Diane about though can't recall exactly why
rownames(GRMCorr) = GRM$id$V2 # change from V1 to V2 as FID is now 0 instead of IID (we want to likely get the IIDs here, which is column 2)
colnames(GRMCorr) = GRM$id$V2

# as this is not a positive definite matrix, smooth correlations
GRMCorr = cor.smooth(GRMCorr)

# build model for different ancestry cutoffs:
# 55%, 45%, 35%, 25%, 15%, 5%
# assess at each cutoff, how many breeds remain
# this was for when I submitted jobs with M as a variable
# for metabolites, we would like to keep all dogs
M = as.numeric(ancestry_cutoff)

# settings for model
min = 0 # minimum pct considered
max = M # maximum pct considered

# get breed calls for dogs listed in the GRM
BreedLoad = copy(breedcalls[, list(pct = mean(pct, na.rm = T)), by = c("dog","breed")])

# remove spaces from breed
BreedLoad[, breed := gsub(" ","_",breed)]

# eliminate calls below min%, set to 0%
BreedLoad[pct < min, pct := 0]

# keep only dogs that have metabolite values
BreedLoad = BreedLoad[dog %in% unique(pheno_df$dog_id)]

# keep only dogs in the GRM
BreedLoad = BreedLoad[dog %in% GRM$id$V2]

# set 0% to 0.01% for purpose of correlation matrix
BreedLoad[pct == 0, pct := 0.01]

# remove dogs with any call not within [min,max] %
BreedLoad = BreedLoad[!dog %in% unique(BreedLoad[pct > max]$dog)]
print("Dogs kept in model:")
print(length(unique(BreedLoad$dog)))

# retain calls from breeds with SD >= 5% of ancestry between 5% and 100%
## 3. only breeds with a SD >5% ancestry across dogs carrying that breed (ensure that we do not include breeds that have only 100% ancestry but capture variation)
BreedLoad = BreedLoad[breed %in% BreedLoad[pct >= 0.05, sd(pct, na.rm = T), by = breed][V1 >= 0.05]$breed]
print("Breeds kept in model:")
print(length(unique(BreedLoad$breed)))

# get correlation
BreedLoadCorr = pivot_wider(data = BreedLoad,
                            id_cols = c("dog"),
                            names_from = "breed",
                            values_from = "pct",
                            values_fill = 0) %>% as.data.table()
BreedLoad_dogs = BreedLoadCorr$dog
BreedLoadCorr[, dog := NULL]
BreedLoadCorr = round(cor(t(BreedLoadCorr)),3)


BreedLoadCorr = cor.smooth(BreedLoadCorr)
rownames(BreedLoadCorr) = BreedLoad_dogs
colnames(BreedLoadCorr) = BreedLoad_dogs

# reset 0.01% to 0%
BreedLoad[pct < min, pct := 0]

# list breeds col names
formula_paste = paste(unique(BreedLoad$breed), collapse = " + ")

# breed info
BreedLoad_breeds = data.table(breed = unique(BreedLoad$breed),
                              Nbreed = BreedLoad[pct > 0, .N, by = "breed"]$N)
breeds_run = unique(BreedLoad$breed)

# scale percent ancestry by breed
BreedLoad[, pct := scale(pct), by = "breed"]

# pivot by breed
BreedLoad = pivot_wider(data = BreedLoad,
                        id_cols = "dog",
                        names_from = "breed",
                        values_from = "pct",
                        values_fill = 0) %>% as.data.table()

# merge with blood phenotype info
BreedLoad = merge(pheno_long_df,
                  BreedLoad,
                  by.x = "dog_id", by.y = "dog", allow.cartesian = T) %>% as.data.table()

# dog info
BreedLoad_dogs = unique(BreedLoad$dog_id)

# dog as factor
BreedLoad$dog_id = as.factor(BreedLoad$dog_id)


print("Dogs kept in model:")
print(length(unique(BreedLoad$dog)))

# age brackets

# put age into 3 equal categories (tertiles/terciles)
BreedLoad$age_bracket <- ntile(BreedLoad$age, 3)


# 
# BreedLoad[, age_bracket := "3yr_to_10yr"]
# BreedLoad[age>=10, age_bracket := "10yr_min"]
# BreedLoad[age<=3, age_bracket := "3yr_max"]

age_bracket_df <- BreedLoad %>% group_by(dog_id) %>% select(dog_id, age, age_bracket) %>% distinct()

table(age_bracket_df$age_bracket)




BreedLoad_PhenoValue_LMER_CovMat_byPhenotype = list()
BreedLoad_PhenoValue_LMER_CovMat_byBreed = list()
modelType = "LMER_CovMat-GRM"


mi = 1 # blood pheno index always 1 since doing this one blood phenotype at a time on server
bi = 1 # breed index for inner loop

formula_allDogs = as.formula(paste("pheno_value ~ (1 | dog_id) + (1 | age_bracket) +", paste(BreedLoad_breeds$breed, collapse = " + ")))
  
# model with REML, all breeds:
model_withBreed_REML = relmatLmer(data = BreedLoad[pheno_name == pheno_outname],
                                    formula = formula_allDogs,
                                    relmat = list(dog_id = GRMCorr),
                                    na.action = "na.exclude",
                                    REML = T)
model_summ = data.table(coef(summary(model_withBreed_REML)), keep.rownames = T)
  
# model with ML, all breeds:
model_withBreed_ML = relmatLmer(data = BreedLoad[pheno_name == pheno_outname],
                                  formula = formula_allDogs,
                                  relmat = list(dog_id = GRMCorr),
                                  na.action = "na.exclude",
                                  REML = F)
  
BreedLoad_PhenoValue_LMER_CovMat_byPhenotype[[mi]] = compare_performance(model_withBreed_ML, model_withBreed_REML) %>% as.data.table()
  
BreedLoad_PhenoValue_LMER_CovMat_byPhenotype[[mi]]$pheno = pheno_outname
  
# perform ANOVA on model with REML to obtain F-stats for each breed:
model_anova = as.data.table(anova(model_withBreed_REML), keep.rownames = T)
  
# perform ANOVA of model with ML +/- breed to get likelihood ratio:
for (B in BreedLoad_breeds$breed){
    print(B)
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]] = data.table(pheno = pheno_outname,
                                                                 breed = B,
                                                                 N = BreedLoad_breeds[breed == B]$Nbreed,
                                                                 modelType = modelType,
                                                                 modelMinPctCutoff = min,
                                                                 modelMaxPctCutoff = max,
								 num_dogs = sum(!is.na(BreedLoad[pheno_name == pheno_outname]$pheno_value)))
    
    model_withoutBreed_ML = relmatLmer(data = BreedLoad[pheno_name == pheno_outname],
                                       formula = as.formula(paste("pheno_value ~ (1 | dog_id) + (1 | age_bracket) +", paste(BreedLoad_breeds[breed!=B]$breed, collapse = " + "))),
                                       relmat = list(dog_id = GRMCorr),
                                       na.action = "na.exclude",
                                       REML = F)
    
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.effect = model_summ[rn==B]$Estimate
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.SD = model_summ[rn==B]$`Std. Error`
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.t.val = model_summ[rn==B]$`t value`
    
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.anova.sumSq = model_anova[rn==B]$`Sum Sq`
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.anova.meanSq = model_anova[rn==B]$`Mean Sq`
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$REML.anova.F.val = model_anova[rn==B]$`F value`
    
    BreedLoad_anova = anova(model_withBreed_ML,model_withoutBreed_ML)
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$ML.anova.chisq = BreedLoad_anova$Chisq[2]
    BreedLoad_PhenoValue_LMER_CovMat_byBreed[[bi]]$ML.anova.p = BreedLoad_anova$`Pr(>Chisq)`[2]
    
    bi = bi+1 
  }


# collapse all breed-factor results:
BreedLoad_PhenoValue_LMER_CovMat_byBreed = rbindlist(BreedLoad_PhenoValue_LMER_CovMat_byBreed) %>% as.data.table()

# adjust p-values from model with ML +/- breed ANOVA:
BreedLoad_PhenoValue_LMER_CovMat_byBreed[, ML.anova.p.adj_bonferroni := p.adjust(ML.anova.p, method = "bonferroni"), by = "pheno"]
BreedLoad_PhenoValue_LMER_CovMat_byBreed[, ML.anova.p.adj_benjhoch_FDR := p.adjust(ML.anova.p, method = "fdr"), by = "pheno"]

BreedLoad_PhenoValue_LMER_CovMat_byBreed$ancestryCutoff = M

# collapse all metabolite results:
BreedLoad_PhenoValue_LMER_CovMat_byPhenotype = rbindlist(BreedLoad_PhenoValue_LMER_CovMat_byPhenotype, fill=TRUE) %>% as.data.table()
BreedLoad_PhenoValue_LMER_CovMat_byPhenotype$ancestryCutoff = M

write.table(x = BreedLoad_PhenoValue_LMER_CovMat_byBreed,
            file = paste(outname, "_", pheno_outname, "_LMER_CovMat_byBreed_ancCutoff-",M,".csv",sep=""),
            sep = ",",
            row.names = F)

write.table(x = BreedLoad_PhenoValue_LMER_CovMat_byPhenotype,
            file = paste(outname,"_", pheno_outname, "_LMER_CovMat_byMetabolite_ancCutoff-",M,".csv",sep=""),
            sep = ",",
            row.names = F)

#save.image(file = paste("DAP_BloodPhenoAnalysis_LMER_ancCutoff-",M,".Rdata", sep = ""))
