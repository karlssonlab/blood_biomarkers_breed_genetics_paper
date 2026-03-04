library(tidyverse)
library(cowplot)

# Output file for supplemental figure PDF
outpdf_supp <- "../../../figures_pdf/fig_S4_TOPHER.pdf"

# =====================================
# Load precompiled supplemental data (creates raw_S2_*, raw_S3_* objects, etc.)
# =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # populate the global environment with objects stored in the RDS list
}

# =====================================
# Build phenotype labels and category annotations used for plotting
# =====================================
all_pheno <-  raw_S2_ALL_PHENOS %>% mutate(phenotype=tolower(phenotype)) %>% select(phenotype,paper_phenotype_category,plot_label) %>% distinct()

# Clean up a few long/awkward plot labels for display
all_pheno <- all_pheno %>% mutate(plot_label=str_replace(plot_label,"mean corpuscular","corpuscular"))
all_pheno <- all_pheno %>% mutate(plot_label=str_replace(plot_label," aminotransferase"," aminotransf."))

# Standardize category labels and create multi-line facet labels
all_pheno <- all_pheno %>% mutate(paper_phenotype_category=str_replace_all(paper_phenotype_category," \\([A-Z][A-Z][A-Z]\\)",""))
all_pheno <- all_pheno %>% mutate(survey_vs_blood=if_else(paper_phenotype_category=="Survey question",paper_phenotype_category,"Intermediate phenotype"))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=str_replace_all(paper_phenotype_category,"\ ","\n"))

# Add spacing for a specific facet label to improve readability
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=str_replace_all(pheno_cat_xlabel_her,"Plasma\nmetabolites","Plasma\nmetabolites\n"))

# Annotate special phenotype subsets in facet labels (owner/environment)
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=if_else(str_detect(phenotype,"mdors"),paste0(pheno_cat_xlabel_her,"\n(owner)"),pheno_cat_xlabel_her))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=if_else(str_detect(phenotype,"^de\\_"),paste0(pheno_cat_xlabel_her,"\n(environ)"),pheno_cat_xlabel_her))

# Keep only the fields needed downstream and ensure uniqueness
all_pheno <- all_pheno %>% select(phenotype,plot_label,paper_phenotype_category,pheno_cat_xlabel_her,survey_vs_blood) %>% distinct()

# Create a second version of facet labels for GWAS figures (remove owner/environment tags, special-case one label)
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=pheno_cat_xlabel_her)
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=if_else(phenotype=="dd_weight_lbs","Body\nweight\nquestion",str_remove(pheno_cat_xlabel_gwas,"\\([a-z]+\\)")))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_remove(pheno_cat_xlabel_gwas,"\n\\(owner\\)"))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_remove(pheno_cat_xlabel_gwas," \n\\(environment\\)"))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_trim(pheno_cat_xlabel_gwas))

# =====================================
# Extract top SNP per region/phenotype and construct a region table with unique IDs
# =====================================

# For each region-phenotype combination, keep the SNP with the smallest P value
top_snps <- raw_S3_GWAS_REGIONS %>% mutate(phenotype = tolower(phenotype)) %>%
  group_by(region, phenotype) %>% slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(region, phenotype, SNP, P)

# Create a distinct region table and assign an ID for stable referencing/order
unique_regions <- raw_S3_GWAS_REGIONS %>% mutate(phenotype = tolower(phenotype)) %>%
  select(CHR, region, region_P, phenotype, region_start, region_end) %>%
  distinct() %>% arrange(region_P, phenotype) %>% mutate(uniqueID = row_number())

# Attach the top SNP info to each region/phenotype row
unique_regions <- unique_regions %>% left_join(top_snps, by = c("region", "phenotype"))

# =====================================
# Plot theme for supplemental figures
# =====================================
theme_supp_fig <- function(){ 
  theme_cowplot(8) %+replace%
    theme(     plot.title = element_text(hjust=0,size=10,face="bold",margin=margin(0,0,4,0)),
               plot.subtitle = element_text(hjust=0,size=8,margin=margin(0,0,4,0),lineheight = 1),
               axis.title.x = element_text(size=9),
               axis.title.y = element_blank(),
               axis.line = element_line(linewidth = 0.25),
               axis.ticks = element_line(linewidth = 0.25),
               
               axis.text.x = element_text(size=8,vjust=1,margin=margin(2,0,0,4),hjust=0.5),
               axis.text.y = element_text(size=8,margin=margin(2,0,0,4),hjust=1),
               panel.grid.major.x=element_line(linewidth=0.2,color="grey70"),
               panel.grid.major.y=element_line(linewidth=0.2,color="grey90"),
               
               legend.position = c(0.02,0.98),
               legend.text = element_text(size=7),legend.title = element_blank(),
               legend.key.size = unit(0.25, 'cm'),
               legend.justification.inside = c(-0.05, 1),
               legend.box.just = "right",
               strip.text.y = element_text(lineheight=0.9,size=8,angle=90,hjust=0.5,face="bold",color="grey10",margin=margin(4,4,4,4))
    )
}

# =====================================
# Prepare heritability table with labels and standard errors
# =====================================
her <- raw_S2_HERITABILITY %>% mutate(phenotype=tolower(phenotype)) %>% left_join(all_pheno)
her <- her %>% select(phenotype,paper_phenotype_category,pheno_cat_xlabel_her,survey_vs_blood,plot_label,GREML_LDS_constrained_heritability,GREML_LDS_constrained_stderr)
her <- her  %>% filter(!is.na(GREML_LDS_constrained_heritability)) %>% mutate(her=GREML_LDS_constrained_heritability,her_err=GREML_LDS_constrained_stderr)

# =====================================
# Build supplemental figure: top heritability traits (80th percentile within group)
# =====================================
pd <- her  %>% filter(!is.na(her)) %>% select(pheno_cat_xlabel_her,plot_label,phenotype,her,her_err,paper_phenotype_category)
pd <- pd %>% rename(her=her,err=her_err)
pd <- pd %>% arrange(-her) %>% mutate(plotorder=row_number())
pd <- pd %>% arrange(-her) %>% mutate(rank=row_number())

# Compute error bar endpoints, clamped to [0, 1]
pd <- pd %>% mutate(xstart=if_else(her-err<0,0,her-err),xend=if_else(her+err>1,1,her+err))
dim(pd)

# Exclude owner/environment phenotypes from this figure
pd <- pd %>% filter(!str_detect(phenotype,"mdors"))
pd <- pd %>% filter(!str_detect(phenotype,"^de\\_"))
dim(pd)

# Harmonize facet labels for display
pd <- pd %>% mutate(survey_vs_blood=if_else(paper_phenotype_category=="Survey question","Survey questions","Intermediate phenotypes"))

# Keep only top 20% heritability phenotypes within each facet group
cutoffs <- pd %>% group_by(survey_vs_blood) %>% summarize(cutoff=quantile(her,0.80))
pd <- pd %>% inner_join(cutoffs) %>% filter(her>=cutoff)

# Determine x-axis lower bound aligned to 0.05 steps (used for consistent plotting if needed)
minx <- floor(min(pd$xstart)/0.05)*0.05

# Order y-axis (phenotype labels) by heritability rank
phenos <- pd %>% arrange(-plotorder) %>% pull(plot_label)
pd$plot_label <- factor(pd$plot_label,levels=phenos)

# =====================================
# Plot: point estimates with error bars, faceted by phenotype group
# =====================================
p <- ggplot(pd,aes(y=plot_label,x=her))
p <- p + geom_point(aes(shape=paper_phenotype_category,fill=paper_phenotype_category,color=paper_phenotype_category))
p <- p + geom_segment(aes(x=xstart,xend=xend,yend=plot_label,color=paper_phenotype_category))

p <- p + scale_color_manual(values = c(
  "Clinical analytes"    = "#f26c2a",
  "Plasma metabolites"   = "#66c2a5",
  "Survey question" = "#e78ac3"
)) 
p <- p + scale_fill_manual(values = c(
  "Clinical analytes"    = "#f26c2a",
  "Plasma metabolites"   = "#66c2a5",
  "Survey question" = "#e78ac3"
  
)) 

p <- p + scale_shape_manual(values = c(
  "Clinical analytes"    = 21,
  "Plasma metabolites"   = 22,
  "Survey question" = 24
)) 
p <- p + facet_grid(survey_vs_blood~.,scales="free_y",space="free_y")
p <- p + scale_x_continuous("heritability",limits=c(0,1))
p <- p + theme_supp_fig() 

# Save figure to PDF
ptop <- p

ggsave(plot=ptop,filename=outpdf_supp,width=6.5,height=6.5)
