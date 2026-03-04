library(tidyverse)
library(ggpubr)
library(cowplot)
library(scales)
source("util.cleaning.R")

# ------------------------------------------------------------------------------
# Script: plot_fig_S1_DATASETS.R
# Purpose: Load cohort metadata and generate a multi-panel supplemental figure.
# Output: Saves a multi-panel PDF to `outpdf`.
# ------------------------------------------------------------------------------

# Output PDF path
outpdf <- "../../../figures_pdf/fig_S1_DATASETS.pdf"

# ------------------------------------------------------------------------------
# Load data objects (only if not already present in the current R environment)
# ------------------------------------------------------------------------------

if (!exists("dat")) {
dat <- readRDS("../../data/DAP_supp_data.rds")
list2env(dat, .GlobalEnv) 
}

# ------------------------------------------------------------------------------
# Standardize breed naming across input tables
# ------------------------------------------------------------------------------

indata <- standardize_breed_names(raw_S14_SEQ_DOGS, standardized_breed)
alldogs <- standardize_breed_names(raw_S13_DAP_ALL, Breed)

# ------------------------------------------------------------------------------
# Clean/derive metadata fields used in plotting
# - Convert "yes"/"no" columns to logical
# - Create combined clinical analytes flag
# ------------------------------------------------------------------------------

indata <- indata %>%
  mutate(across(
    where(~ all(.x %in% c("yes", "no"))),
    ~ .x == "yes"
  ))
indata <- indata %>% mutate(clinical_analytes=if_else(serum_chemistry_panel|complete_blood_counts,TRUE,FALSE)) %>% select(-serum_chemistry_panel,-complete_blood_counts)

# ------------------------------------------------------------------------------
# Harmonize variable names/types and define dataset labels
# ------------------------------------------------------------------------------

alldogs <- alldogs %>% rename(age_years=Estimated_Age_Years_at_HLES)
alldogs <- alldogs %>% mutate(dog_id=as.integer(dog_id),age_years=as.numeric(age_years))
alldogs <- alldogs %>% mutate(set="all dogs")

# Keep only sequenced dogs with required identifiers/measurements for plotting
indata <- indata %>%  filter(!is.na(dog_id) & !is.na(weight_kg))
indata <- indata %>% mutate(source="seqDAP")
indata <- indata %>% rename(age_years=Estimated_Age_Years_at_HLES)

# ------------------------------------------------------------------------------
# Define and build "set" membership table for sequenced dogs
# Each spec yields a label for a subset used as a y-axis group in plots.
# ------------------------------------------------------------------------------

label_specs <- list(
  list(filter = NULL, set = "sequenced dogs", order = 0),
  list(filter = quote(survey_data), set = "seq. dogs w/ survey data", order = 2),
  list(filter = quote(clinical_analytes), set = "seq. dogs w/ clinical analytes", order = 2),
  list(filter = quote(plasma_metabolites), set = "seq. dogs w/ plasma metabolites", order = 2)
)

# Apply each specification to indata and stack results
labels <- map_dfr(label_specs, function(spec) {
  df <- if (is.null(spec$filter)) indata else indata %>% filter(!!spec$filter)
  df %>%
    select(dog_id, source) %>%
    mutate(set = spec$set, order = spec$order)
})

# Add the all-dogs cohort as an additional "set"
labels <- bind_rows(
  labels,
  alldogs %>%
    mutate(source = "allDAP") %>%
    select(dog_id, source) %>%
    mutate(set = "all dogs", order = -1)
)

# Keep only columns needed downstream and define single vs mixed-breed label
indata <- indata %>% select(dog_id,source,age_years,Breed=standardized_breed,single_breed,sex,sterilization_status,weight_kg,top_breed,max_percent_ancestry,num_breeds_greater_than_five_percent)
indata <- indata %>% mutate(breed_type=if_else(single_breed,"single-breed","mixed breed"))

# ------------------------------------------------------------------------------
# Plot theme helper (shared across box/violin/bar plots)
# ------------------------------------------------------------------------------

theme_histo <- function(){
  theme_cowplot(12) %+replace%
    theme(axis.title = element_text(size=8),axis.text = element_text(size=7),
          legend.position = "none",
          legend.text = element_text(size=5.5),legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.justification.inside = c(1.05, 1.25),
          legend.box.just = "right",
          axis.ticks.y=element_blank(),
          strip.text.x=element_text(size=7,hjust=0.5,face="bold",lineheight=1), margin(5,3,5,3)
    )
}

# ------------------------------------------------------------------------------
# Build plotting tables for comparisons across dataset "sets"
# ------------------------------------------------------------------------------

# Determine order and display labels (including N) for each set on the y-axis
counts <- labels %>% group_by(set,order) %>% count()  %>% arrange(-order,n) %>% ungroup() %>% mutate(order2=row_number())

# Age/weight: combine sequenced-dog subsets with all-dogs cohort and track source via alpha
pd <- indata %>% left_join(labels)
pd <- pd %>% select(dog_id,set,weight_kg,age_years)  %>% mutate(source="seqDAP") 
pd <- alldogs %>% select(dog_id,set,weight_kg,age_years)  %>% mutate(source="allDAP") %>% bind_rows(pd)
pdAge <- pd %>% select(dog_id,set,age_years,source) %>% filter(!is.na(age_years))
pdWeight <- pd %>% select(dog_id,set,weight_kg,source) %>% filter(!is.na(weight_kg))

# Sex/sterilization: compute within-set fractions for each category and source
pd <- indata  %>% left_join(labels) 
pd <- pd %>% select(dog_id,set,sex,sterilization_status) %>% mutate(source="seqDAP")
pd <- alldogs %>% select(dog_id,set,sex,sterilization_status) %>% mutate(source="allDAP") %>% bind_rows(pd) 
pd <- pd %>% pivot_longer(c(-dog_id,-set,-source)) %>% filter(!is.na(value)) %>% group_by(set,name,value,source) %>% count() 
pd <- pd %>% group_by(set,name,source) %>% summarize(ntot=sum(n)) %>% full_join(pd) %>% mutate(frac=n/ntot)

pdSex <- pd %>% ungroup() %>% filter(name=="sex") %>% select(source,set,frac,value) %>% rename(sex=value)
pdSter <- pd %>% ungroup() %>% filter(name=="sterilization_status") %>% select(source,set,frac,value) %>% rename(sterilization_status=value) 

# Create y-axis labels and breaks for plots (ordered by `counts$order2`)
counts <- counts %>% mutate(ylabel=paste(set," (N=",comma(n, accuracy = 1),")",sep="")) 
ylabels <- counts %>% arrange(order2) %>% pull(ylabel)
ybreaks <- counts %>% arrange(order2) %>% pull(set)

# Comparisons are defined by y-axis category labels (used by stat_compare_means)
my_comparisons <- list(c(ybreaks[3],ybreaks[5]),c(ybreaks[3],ybreaks[1]),c(ybreaks[3],ybreaks[2]))

# Lock factor levels to the computed y-axis ordering
pdAge$set <- factor(pdAge$set,levels=ybreaks)
pdWeight$set <- factor(pdWeight$set,levels=ybreaks)
pdSex$set <- factor(pdSex$set,levels=ybreaks)
pdSter$set <- factor(pdSter$set,levels=ybreaks)

# ------------------------------------------------------------------------------
# Panel A: Age distribution by dataset set (seqDAP vs allDAP via alpha)
# ------------------------------------------------------------------------------

p <- ggplot(pdAge,aes(x=age_years,y=set)) 
p <- p + geom_violin(aes(alpha=source),color=NA,fill="#969696")
p <- p + geom_boxplot(aes(alpha=source),fill=NA,outlier.size=0.5,outlier.alpha=0.25,outlier.shape=16,linewidth=0.25,width=0.5) 
p <- p + scale_alpha_manual(values=c(0.9,0.2))
p <- p + stat_compare_means(comparisons = my_comparisons,size=1.75,method="wilcox.test")
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels,expand = expansion(mult = c(0.01, 0.1)))
p <- p + scale_x_continuous("age (years)")
p <- p + theme_histo()
pbox_Age <- p 

# ------------------------------------------------------------------------------
# Panel A (right): Weight distribution by dataset set
# ------------------------------------------------------------------------------

p <- ggplot(pdWeight,aes(x=weight_kg,y=set)) 
p <- p + geom_violin(aes(alpha=source),color=NA,fill="#969696")
p <- p + geom_boxplot(aes(alpha=source),fill=NA,outlier.size=0.5,outlier.alpha=0.25,outlier.shape=16,linewidth=0.25,width=0.5) 
p <- p + stat_compare_means(comparisons = my_comparisons,size=1.5)
p <- p + scale_alpha_manual(values=c(0.9,0.2))
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels,expand = expansion(mult = c(0.01, 0.1)))
p <- p + scale_x_continuous("weight (kg)")
p <- p + theme_histo()
p <- p + theme(axis.title.y=element_text(margin=margin(l=0)),
               axis.text.y=element_blank())
pbox_Weight <- p 

# ------------------------------------------------------------------------------
# Panel B: Sex composition by set (fractions)
# ------------------------------------------------------------------------------

p <- ggplot(pdSex,aes(x=frac,y=set))
p <- p + geom_bar(aes(fill=factor(sex),alpha=source),width=0.6,stat = "identity",color="#252525",linewidth=0.25)
p <- p + scale_fill_manual(values=c("#969696","#FFFFFF"))
p <- p + scale_alpha_manual(values=c(1,0.5))
p <- p + geom_text(aes(label=paste(sex," (",round(frac*100,1),"%)",sep=""),x=0.01),lineheight=0.9,hjust=0,vjust = 0.5,size=2,color="#252525",data=pdSex %>% filter(sex=="male"))
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("fraction of dogs",breaks=c(0,0.5,1))
p <- p + theme_histo()
pbar_Sex <- p

# ------------------------------------------------------------------------------
# Panel B (right): Sterilization composition by set (fractions)
# ------------------------------------------------------------------------------

p <- ggplot(pdSter,aes(x=frac,y=set))
p <- p + geom_bar(aes(fill=factor(sterilization_status),alpha=source),width=0.6,stat = "identity",color="#252525",linewidth=0.25)
p <- p + scale_fill_manual(values=c("#969696","#FFFFFF"))
p <- p + scale_alpha_manual(values=c(1,0.5))
p <- p + geom_text(aes(label=paste(sterilization_status," (",round(frac*100,1),"%)",sep=""),x=0.01),lineheight=0.9,hjust=0,vjust = 0.5,size=2,color="#252525",data=pdSter %>% filter(sterilization_status=="sterilized"))
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("fraction of dogs",breaks=c(0,0.5,1))
p <- p + theme_histo()
p <- p + theme(axis.text.y=element_blank(),
               axis.title.y=element_text(margin=margin(l=0)))
pbar_Ster <- p

# ------------------------------------------------------------------------------
# Compare mixed vs single-breed within sequenced dogs
# - Age/weight: t-tests
# - Sex/sterilization: chi-squared tests; sex test also stratified by sterilization
# ------------------------------------------------------------------------------

pd <- indata %>% filter(!is.na(breed_type)) 
pd <- pd %>% select(dog_id,breed_type,weight_kg,age_years,sex,sterilization_status)
pdAge_mix <- pd %>% select(dog_id,breed_type,age_years) %>% filter(!is.na(age_years))
pdWeight_mix <- pd %>% select(dog_id,breed_type,weight_kg) %>% filter(!is.na(weight_kg))
pdSex_mix <- pd %>% select(dog_id,breed_type,sex) %>% filter(!is.na(sex))
pdSter_mix <- pd %>% select(dog_id,breed_type,sterilization_status) %>% filter(!is.na(sterilization_status))
pdSex_intact <- pdSter_mix %>% filter(sterilization_status!="sterilized") %>% select(dog_id) %>% inner_join(pdSex_mix)
pdSex_ster <- pdSter_mix %>% filter(sterilization_status=="sterilized") %>% select(dog_id) %>% inner_join(pdSex_mix)

# Build y-axis labels with counts for breed_type comparisons (plus a blank row for p-value text)
counts <- pd %>% group_by(breed_type) %>% count() %>% mutate(ylabel=paste0(breed_type," (N=",n,")"))
counts <- counts %>% arrange(breed_type)
ylabels <- c(counts %>% pull(ylabel),"")
ybreaks <- c(counts %>% pull(breed_type),"")

# ------------------------------------------------------------------------------
# Panel C: Age by breed_type (t-test)
# ------------------------------------------------------------------------------

pdAge_mix_ttest <- t.test(age_years ~ breed_type, data = pdAge_mix)
pdAge_mix.pval <- pdAge_mix_ttest$p.value
pdAge_mix.pval <- paste0("t-test p = ", signif(pdAge_mix.pval, 2))
pdAge_mix$breed_type <- factor(pdAge_mix$breed_type,levels=ybreaks)
p <- ggplot(pdAge_mix,aes(x=age_years,y=breed_type)) 
p <- p + geom_boxplot(aes(fill=breed_type),outlier.size=0.5,outlier.alpha=0.25,outlier.shape=16,linewidth=0.25,width=0.75) 
p <- p + geom_text(aes(x = 0, y = breed_type, label = pdAge_mix.pval),data = tibble(breed_type = ""),hjust = 0, size = 2)
p <- p + scale_fill_manual(values=c("single-breed"="grey100","mixed breed"="#969696"))
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("age (years)")
p <- p + theme_histo()
p <- p + theme(legend.position = "none")
pbox_Age_mix <- p 

# ------------------------------------------------------------------------------
# Panel D: Weight by breed_type (t-test)
# ------------------------------------------------------------------------------

pdWeight_mix_ttest <- t.test(weight_kg ~ breed_type, data = pdWeight_mix)
pdWeight_mix.pval <- pdWeight_mix_ttest$p.value
pdWeight_mix.pval <- paste0("t-test p = ", signif(pdWeight_mix.pval, 2))
pdWeight_mix$breed_type <- factor(pdWeight_mix$breed_type,levels=ybreaks)

p <- ggplot(pdWeight_mix,aes(x=weight_kg,y=breed_type))
p <- p + geom_boxplot(aes(fill=breed_type),outlier.size=0.5,outlier.alpha=0.25,outlier.shape=16,linewidth=0.25,width=0.75) 
p <- p + geom_text(aes(x = 0, y = breed_type, label = pdWeight_mix.pval),data = tibble(breed_type = ""),hjust = 0, size = 2)
p <- p + scale_fill_manual(values=c("single-breed"="grey100","mixed breed"="#969696"))
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("weight (kg)")
p <- p + theme_histo()
p <- p + theme(axis.text.y=element_blank())
pbox_Weight_mix <- p 

# ------------------------------------------------------------------------------
# Panel E: Fraction male by breed_type (chi-squared; also computed within sterilized dogs)
# ------------------------------------------------------------------------------

pdSex_mix_chisq <- chisq.test(table(pdSex_mix$breed_type, pdSex_mix$sex))
pdSex_mix.pval <- pdSex_mix_chisq$p.value 
pdSex_ster_chisq <- chisq.test(table(pdSex_ster$breed_type, pdSex_ster$sex))
pdSex_mix.pval.ster <- pdSex_ster_chisq$p.value
pdSex_mix.pval <- paste0("chisq p = ", signif(pdSex_mix.pval, 2)," (p = ", signif(pdSex_mix.pval.ster, 2)," in sterilized dogs)")

pdTot <- pdSex_mix %>% count(breed_type, name = "n_tot")
pdSex_mix <- pdSex_mix %>% filter(sex == "male") %>% count(breed_type, name = "n_male")
pdSex_mix <- pdSex_mix %>% left_join(pdTot) %>% mutate(frac=n_male/n_tot)
pdSex_mix$breed_type <- factor(pdSex_mix$breed_type,levels=ybreaks)

# Overlay intact-dog fraction males (semi-transparent) on top of overall bars
pdSex_intact <- pdSex_intact %>% filter(sex == "male") %>% count(breed_type, name = "n_male")
pdSex_intact <- pdSex_intact %>% left_join(pdTot) %>% mutate(frac=n_male/n_tot)
pdSex_intact$breed_type <- factor(pdSex_intact$breed_type,levels=ybreaks)

p <- ggplot(pdSex_mix,aes(x=frac,y=breed_type))
p <- p + geom_col(aes(fill = breed_type),color="#252525",width=0.75)
p <- p + geom_col(fill = "#B2182B",color=NA,width=0.75,alpha=0.25,data=pdSex_intact)
p <- p + scale_fill_manual(values=c("single-breed"="grey100","mixed breed"="#969696"))
p <- p + geom_text(aes(x = 0, y = breed_type, label = pdSex_mix.pval),data = tibble(breed_type = ""),hjust = 0, size = 2)
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("fraction males",limits=c(0,0.54))
p <- p + theme_histo()
pbar_Sex_mix <- p 

# ------------------------------------------------------------------------------
# Panel F: Fraction intact by breed_type (chi-squared)
# ------------------------------------------------------------------------------

pdSter_mix_chisq <- chisq.test(table(pdSter_mix$breed_type, pdSter_mix$sterilization_status))
pdSter_mix.pval <- pdSter_mix_chisq$p.value
pdSter_mix.pval <- paste0("chisq p = ", signif(pdSter_mix.pval, 2))
pdTot <- pdSter_mix %>% count(breed_type, name = "n_tot")
pdSter_mix <- pdSter_mix %>% filter(sterilization_status == "intact") %>% count(breed_type, name = "n_intact")
pdSter_mix <- pdSter_mix %>% left_join(pdTot) %>% mutate(frac=n_intact/n_tot)
pdSter_mix$breed_type <- factor(pdSter_mix$breed_type,levels=ybreaks)

p <- ggplot(pdSter_mix,aes(x=frac,y=breed_type))
p <- p + geom_col(aes(fill = breed_type),color="#252525",width=0.75)
p <- p + scale_fill_manual(values=c("single-breed"="grey100","mixed breed"="#969696"))
p <- p + geom_text(aes(x = 0, y = breed_type, label = pdSter_mix.pval),data = tibble(breed_type = ""),hjust = 0, size = 2)
p <- p + scale_y_discrete("",breaks=ybreaks,labels=ylabels)
p <- p + scale_x_continuous("fraction intact")
p <- p + theme_histo()
p <- p + theme(axis.text.y=element_blank())
pbar_Ster_mix <- p 

# ------------------------------------------------------------------------------
# Assemble multi-panel figure and save
# ------------------------------------------------------------------------------

row1 <- plot_grid(pbox_Age,pbox_Weight,ncol=2,rel_widths=c(1,0.5),labels=c("A"),label_size = 10)
row2 <- plot_grid(pbar_Sex,pbar_Ster,ncol=2,rel_widths=c(1,0.5),labels=c("B"),label_size = 10)
row3 <- plot_grid(pbox_Age_mix,pbox_Weight_mix,ncol=2,rel_widths=c(1,0.75),labels=c("C","D"),label_size = 10)
row4 <- plot_grid(pbar_Sex_mix,pbar_Ster_mix,ncol=2,rel_widths=c(1,0.75),labels=c("E","F"),label_size = 10)

grid <- plot_grid(row1,row2,row3,row4,ncol=1,rel_heights = c(1.75,1.75,1,1))
ggsave(plot=grid,filename=outpdf,width=6,height=5.5)

# ------------------------------------------------------------------------------
# Print formatted test results (for logs / reproducibility)
# ------------------------------------------------------------------------------

format_test_result(pdAge_mix_ttest)
format_test_result(pdWeight_mix_ttest)
format_test_result(pdSex_mix_chisq)
format_test_result(pdSex_ster_chisq)
format_test_result(pdSter_mix_chisq)
