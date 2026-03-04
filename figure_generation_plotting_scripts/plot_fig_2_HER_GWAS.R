library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(magick)
library(lme4)
library(emmeans)

# Project helper functions (e.g., pdf_to_plot and any cleaning utilities)
source("util.cleaning.R")

# Output figure path
outpdf <- "../../../figures_pdf/fig_2_HER_GWAS.pdf"

## =====================================
## Load all data (from precompiled RDS bundle)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # optional: creates raw_pheno, raw_anc, etc.
}

## =====================================
## Build phenotype label table for plotting
## - standardize naming/casing
## - create category labels for heritability and GWAS panels
## - add sample-size ranges per category for GWAS x-axis labels
## =====================================
all_pheno <-  raw_S2_ALL_PHENOS %>% mutate(phenotype=tolower(phenotype)) %>% select(phenotype,paper_phenotype_category,plot_label,sample_size) %>% distinct()

# Adjust phenotype plot labels for compact figure text
all_pheno <- all_pheno %>% mutate(plot_label=str_replace(plot_label,"mean corpuscular","corpuscular"))
all_pheno <- all_pheno %>% mutate(plot_label=str_replace(plot_label," aminotransferase"," aminotransf."))

# Clean category names and define survey vs intermediate (blood) grouping
all_pheno <- all_pheno %>% mutate(paper_phenotype_category=str_replace_all(paper_phenotype_category," \\([A-Z][A-Z][A-Z]\\)",""))
all_pheno <- all_pheno %>% mutate(survey_vs_blood=if_else(paper_phenotype_category=="Survey question",paper_phenotype_category,"Intermediate phenotype"))

# Line-break category labels for plotting
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=str_replace_all(paper_phenotype_category,"\ ","\n"))

# Add small qualifiers for specific survey subtypes to match manuscript labeling
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=if_else(str_detect(phenotype,"mdors"),paste0(pheno_cat_xlabel_her,"\n(owner)"),pheno_cat_xlabel_her))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=if_else(str_detect(phenotype,"^de\\_"),paste0(pheno_cat_xlabel_her,"\n(environment)"),pheno_cat_xlabel_her))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_her=if_else(pheno_cat_xlabel_her=="Survey\nquestion",paste0(pheno_cat_xlabel_her,"\n(dog)"),pheno_cat_xlabel_her))

# Final phenotype label table used across plots
all_pheno <- all_pheno %>% select(phenotype,plot_label,paper_phenotype_category,pheno_cat_xlabel_her,sample_size,survey_vs_blood) %>% distinct()

# GWAS x-axis labels differ slightly (collapse/rename survey labels and weight question)
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=pheno_cat_xlabel_her)
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=if_else(phenotype=="dd_weight_lbs","Body\nweight\nquestion",str_remove(pheno_cat_xlabel_gwas,"\\([a-z]+\\)")))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_remove(pheno_cat_xlabel_gwas,"\n\\(owner\\)"))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_remove(pheno_cat_xlabel_gwas," \n\\(environment\\)"))
all_pheno <- all_pheno %>% mutate(pheno_cat_xlabel_gwas=str_trim(pheno_cat_xlabel_gwas))

# Compute sample-size range per GWAS label and attach as multi-line x-axis text
all_pheno <- all_pheno %>% group_by(pheno_cat_xlabel_gwas) %>% 
  summarize(min_sample=min(sample_size),max_sample=max(sample_size)) %>% 
  mutate(size_str=if_else(min_sample==max_sample,paste0("N=",comma(min_sample)),
                                            paste0("N=",comma(min_sample),"-",comma(max_sample)))) %>% 
  mutate(pheno_cat_xlabel_gwas_size=paste(pheno_cat_xlabel_gwas,size_str,sep="\n")) %>% 
  select(pheno_cat_xlabel_gwas,pheno_cat_xlabel_gwas_size) %>% 
  mutate(pheno_cat_xlabel_gwas_size=str_replace_all(pheno_cat_xlabel_gwas_size,"Survey\nquestion","Other\nsurvey\nquestions")) %>% 
  right_join(all_pheno)

## =====================================
## Summarize GWAS regions:
## - pick top SNP per (region, phenotype)
## - create distinct region-phenotype entries with unique IDs
## =====================================
top_snps <- raw_S3_GWAS_REGIONS %>% mutate(phenotype = tolower(phenotype)) %>%
  group_by(region, phenotype) %>% slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(region, phenotype, SNP, P)

unique_regions <- raw_S3_GWAS_REGIONS %>% mutate(phenotype = tolower(phenotype)) %>%
  select(CHR, region, region_P, phenotype, region_start, region_end) %>%
  distinct() %>% arrange(region_P, phenotype) %>% mutate(uniqueID = row_number())

unique_regions <- unique_regions %>% left_join(top_snps, by = c("region", "phenotype"))

## =====================================
## Figure themes (consistent styling across panels)
## =====================================
theme_boxes <- function(){ 
  theme_cowplot(8) %+replace%
    theme(     plot.title=element_text(size=6.5,hjust=0,face="bold",margin=margin(b=5)),
               axis.title.y = element_text(angle=90,size=6,face="bold"),
               axis.title.x = element_blank(),
               axis.line.y = element_line(linewidth = 0.25),
               axis.ticks.y = element_line(linewidth = 0.25),
               axis.line.x = element_blank(),
               axis.text.x = element_text(size=6,vjust=1,margin=margin(t=5),hjust=0.5),
               axis.text.y = element_text(size=6,margin=margin(2,0,0,4),hjust=1),
               panel.grid.major.y=element_line(linewidth=0.25,color="grey80"),
               axis.ticks.x=element_blank(),
               legend.position = "none",
               plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}

theme_line <- function() { 
  theme_boxes() %+replace%
    theme(
      axis.text.x = element_text(size = 5.5, vjust = 1, margin = margin(t=3,b=5), hjust = 0.5,lineheight=1),
      axis.line = element_line(linewidth = 0.25),
      legend.position = "top",
      legend.justification = c("left", "top"),
      legend.direction = "horizontal",
      legend.box.just = "left",
      legend.key.spacing.x = unit(0.1, "cm"),
      legend.key.spacing.y = unit(0.02, "cm"),
      legend.margin = margin(1, 1, 1, 1),
      legend.spacing = unit(0.15, "cm"),
      legend.text = element_text(size = 5, vjust = 0.5),
      legend.key.size = unit(0.15, "cm"),
      legend.title = element_blank(),
      legend.box.background = element_rect(color = "grey40", fill = NA, linewidth = 0.25),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
      plot.margin = margin(t = 3, r = 5, b = 0, l = 0)
    )
}

theme_dots <- function() { 
  theme_boxes() %+replace%
    theme(
      # Axis formatting
      axis.title.x = element_text(
        size = 6,
        face = "bold",
        margin = margin(t = 0, b = 10)
      ),
      axis.text.x = element_text(
        size = 5.5,
        vjust = 1,
        margin = margin(t = 2, b = 2),
        hjust = 0.5,
        lineheight = 1
      ),
      axis.line = element_line(linewidth = 0.25),
      # Plot title
      plot.title = element_text(size = 5.5,hjust = 0,face = "plain",margin = margin(t = 4, b = 4)),
      # Legend inside lower-right
      legend.position = c(0.98, 0.02),
      legend.justification = c("right", "bottom"),
      legend.direction = "vertical",
      legend.box.just = "right",
      legend.key.spacing.x = unit(0.05, "cm"),
      legend.key.spacing.y = unit(0.05, "cm"),
      legend.margin = margin(2, 2, 2, 2),
      legend.spacing = unit(0.05, "cm"),
      legend.text = element_text(size = 5, vjust = 0.5),
      legend.key.size = unit(0.15, "cm"),
      legend.title = element_blank(),
      legend.box.background = element_rect(
        color = "grey40", fill = NA, linewidth = 0.25
      ),
      # Grid + margins
      panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
      plot.margin = margin(t = 3, r = 8, b = 0, l = 0)
    )
}

## =====================================
## Heritability data prep + ordering for boxplots
## =====================================
her <- raw_S2_HERITABILITY %>% mutate(phenotype=tolower(phenotype)) %>% left_join(all_pheno)
her <- her %>% select(phenotype,paper_phenotype_category,pheno_cat_xlabel_her,survey_vs_blood,plot_label,GREML_LDS_constrained_heritability,GREML_LDS_constrained_stderr,min_gwas_pvalue)
her <- her  %>% filter(!is.na(GREML_LDS_constrained_heritability)) %>% mutate(her=GREML_LDS_constrained_heritability,her_err=GREML_LDS_constrained_stderr)

# Order x-axis categories by median heritability (with stable tie-breaking)
ordered_xlabels <- her %>% group_by(paper_phenotype_category,pheno_cat_xlabel_her) %>% summarize(median=median(her)) %>% arrange(-median) %>% ungroup() %>% mutate(order=row_number())
ordered_xlabels <- ordered_xlabels %>% full_join(all_pheno %>% select(paper_phenotype_category,pheno_cat_xlabel_her,pheno_cat_xlabel_gwas) %>% distinct())
ordered_xlabels <- ordered_xlabels %>% ungroup() %>% arrange(order,pheno_cat_xlabel_gwas) %>% mutate(order=row_number())
ordered <- ordered_xlabels %>% arrange(order) %>% select(pheno_cat_xlabel_her) %>% distinct() %>% pull(pheno_cat_xlabel_her)

## =====================================
## Permutation test utility: two-group difference in median/mean
## =====================================
perm_p_two_group <- function(y, g, ref, B = 100000, seed = 1, stat = c("median","mean")) {
  stat <- match.arg(stat)
  set.seed(seed)
  
  g <- factor(g)
  y <- y[!is.na(y) & !is.na(g)]
  g <- g[!is.na(y) & !is.na(g)]
  
  y_ref <- y[g == ref]
  y_x   <- y[g != ref]
  
  stopifnot(length(y_ref) > 0, length(y_x) > 0)
  
  f_stat <- if (stat == "median") median else mean
  obs <- f_stat(y_ref) - f_stat(y_x)
  
  # Permute group labels to get a two-sided permutation p-value
  g_perm <- g
  perm <- replicate(B, {
    g_perm <- sample(g, replace = FALSE)
    f_stat(y[g_perm == ref]) - f_stat(y[g_perm != ref])
  })
  
  p_perm <- (sum(abs(perm) >= abs(obs)) + 1) / (B + 1)
  
  list(p_perm = p_perm, delta = obs)
}

## =====================================
## Panel A: Heritability distribution for survey questions
## - compute permutation p-values vs reference group
## - annotate brackets + p-values
## =====================================
pd_her_survey <- her %>% filter(survey_vs_blood=="Survey question")
ordered_here <- ordered_xlabels %>% arrange(order) %>% select(pheno_cat_xlabel_her) %>% filter(pheno_cat_xlabel_her %in% pd_her_survey$pheno_cat_xlabel_her) %>% distinct() %>% pull(pheno_cat_xlabel_her)
pd_her_survey$pheno_cat_xlabel_her <- factor(pd_her_survey$pheno_cat_xlabel_her,levels=ordered_here)

# Label a couple of highest-heritability survey traits
labelled_survey <- pd_her_survey %>% arrange(-her) %>% slice_head(n=2) %>% mutate(plot_label=str_replace(plot_label," person","\nperson"))

# Reference group and permutation settings for within-survey comparisons
ref <- "Survey\nquestion\n(dog)"
B <- 100000

levels_x <- levels(factor(pd_her_survey$pheno_cat_xlabel_her))

# Precompute permutation p-values (cached in perm_tbl if already present)
if(!exists("perm_tbl")){
perm_tbl <- tibble(pheno_cat_xlabel_her = levels_x) %>%
  filter(pheno_cat_xlabel_her != ref) %>%
  rowwise() %>%
  mutate(
    res = list(
      perm_p_two_group(
        y   = pd_her_survey$her[pd_her_survey$pheno_cat_xlabel_her %in% c(ref, pheno_cat_xlabel_her)],
        g   = pd_her_survey$pheno_cat_xlabel_her[pd_her_survey$pheno_cat_xlabel_her %in% c(ref, pheno_cat_xlabel_her)],
        ref = ref,
        B   = B,
        seed = 1,
        stat = "median"
      )
    ),
    p_perm = res$p_perm,
    delta  = res$delta
  ) %>%
  ungroup() %>%
  mutate(
    p_lab = case_when(
      p_perm <= 1/B ~ paste0("p≤",formatC(p_perm, format="fg", digits=2)),
      p_perm < 0.00000001 ~ paste0("p=", formatC(p_perm, format="fg", digits=2)),
      TRUE ~ paste0("p=", signif(p_perm, 2))
    ),
    group1 = ref,
    group2 = pheno_cat_xlabel_her
  )
}

# x-axis labels show per-category counts
xlabels <- pd_her_survey %>% group_by(pheno_cat_xlabel_her) %>% count()
xlabels <- xlabels %>% mutate(label_n=paste0(pheno_cat_xlabel_her,"\nN=",n)) %>% 
  mutate(label_n=str_replace(label_n,"Survey\nquestion\n","About\n")) %>% 
  mutate(label_n=str_replace_all(label_n,"[\\(\\)]",""))

# Map factor levels to numeric x positions for bracket drawing
x_map <- tibble(pheno_cat_xlabel_her = levels_x, x = seq_along(levels_x))

ann <- perm_tbl %>%
  left_join(x_map %>% rename(group1 = pheno_cat_xlabel_her, x1 = x), by = "group1") %>%
  left_join(x_map %>% rename(group2 = pheno_cat_xlabel_her, x2 = x), by = "group2") %>%
  mutate(
    x_min = pmin(x1, x2),
    x_max = pmax(x1, x2)
  )

# Evenly spaced y-positions for brackets at top of panel
ys <- seq(0.84, 0.90, length.out = nrow(ann))
ann <- ann %>% mutate(y = ys, x_mid = (x_min + x_max)/2)

# Build survey heritability boxplot + bracket annotations
p <- ggplot(pd_her_survey, aes(x = pheno_cat_xlabel_her, y = her))
p <- p + geom_boxplot(fill = "#e78ac3", color = "#e78ac3", alpha = 0.25, width = 0.4, outlier.shape = NA)
p <- p + geom_point(color = "grey30", shape = 20, size = 0.75, alpha = 0.4)
p <- p + geom_segment(data = ann,aes(x = x_min, xend = x_max, y = y, yend = y),
    inherit.aes = FALSE,linewidth = 0.2,color = "grey30")
p <- p + geom_text(data = ann,aes(x = x_mid, y = y + 0.01, label = p_lab),
    inherit.aes = FALSE,size = 1.75,color = "grey20",vjust=0)
p <- p + geom_text_repel(
  aes(label = plot_label),
  color = "#e78ac3", lineheight = 0.9, segment.color = "grey40",
  data = labelled_survey, hjust = 0, vjust = 1, nudge_x = 0.15, direction = "y",
  force = 2, segment.size = 0.25, min.segment.length = 0, max.overlaps = 50, size = 1.75
)
p <- p + ggtitle("Survey questions")
p <- p + scale_x_discrete("",labels=xlabels$label_n,breaks=xlabels$pheno_cat_xlabel_her,expand = expansion(mult = c(0.2, 0.2)))
p <- p + scale_y_continuous(expression("heritability (" * h[SNP]^2 * ")"), limits = c(0, 1), expand = expansion(mult = c(0, 0)))
p <- p + theme_boxes() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10)
  )

pbox_her_survey <- p
p

## =====================================
## Panel A companion: Heritability boxplots for blood traits
## =====================================
pd_her_blood <- her %>% filter(paper_phenotype_category!="Survey question") %>% distinct() 

# Label a subset of highest-heritability blood traits
labelled_blood <- pd_her_blood %>% arrange(-her) %>% slice_head(n=10)
labelled_blood <- labelled_blood %>% mutate(plot_label=str_replace(plot_label," ","
"))
labelled_blood <- labelled_blood %>% mutate(plot_label=str_replace(plot_label,"droxy","droxy
"))

# Use global category order derived above
pd_her_blood$pheno_cat_xlabel_her <- factor(pd_her_blood$pheno_cat_xlabel_her,levels=ordered)

# x-axis labels show per-category counts
xlabels <- pd_her_blood %>% group_by(pheno_cat_xlabel_her) %>% count()
xlabels <- xlabels %>% mutate(label_n=paste0(pheno_cat_xlabel_her,"\nN=",n))

p <- ggplot(pd_her_blood,aes(x=pheno_cat_xlabel_her,y=her))
p <- p + geom_boxplot(aes(fill= paper_phenotype_category,color = paper_phenotype_category),alpha= 0.25,width= 0.3,outlier.shape = NA)
p <- p + geom_point(color="grey30",shape=20,size=0.75,alpha=0.4)
p <- p + ggtitle("Blood traits")
p <- p + geom_text_repel(aes(label = plot_label,color=paper_phenotype_category),lineheight=0.75,
                         segment.color = "grey40",data=labelled_blood,hjust = 1, vjust=0.5, nudge_x = -0.2,direction = "y",   
                         force = 5,    segment.size = 0.25,    min.segment.length = 0,max.overlaps =50 , size=1.75)
p <- p + scale_y_continuous(expression("heritability (" * h[SNP]^2 * ")"),limits=c(0,1),expand = expansion(mult = c(0, 0)))
p <- p + scale_color_manual(values = c("Clinical analytes" = "#f26c2a","Plasma metabolites" = "#2a8c7a")) 
p <- p + scale_fill_manual(values = c("Clinical analytes" = "#f26c2a","Plasma metabolites" = "#2a8c7a")) 
p <- p + scale_x_discrete("",labels=xlabels$label_n,breaks=xlabels$pheno_cat_xlabel_her,
                          expand = expansion(mult = c(0.75,0.1)))
p <- p + theme_boxes() + theme(legend.position = "none")
pbox_her_blood <- p 

## =====================================
## GWAS association strength by trait group (mixed model + EMMs)
## - response: region-level -log10(P)
## - predictors: trait group label and sample size
## - random intercept per phenotype (handles multiple regions per trait)
## =====================================
compare_types <- unique_regions %>% select(phenotype,region,region_P) %>% distinct()
compare_types <- all_pheno %>% select(phenotype,pheno_cat_xlabel_gwas_size,sample_size) %>% 
  left_join(compare_types) %>% 
  mutate(log10P = -log10(region_P))

fit_lmm <- lmer(log10P ~ pheno_cat_xlabel_gwas_size + log(sample_size) + (1 | phenotype), data = compare_types)
fit_lm <- lm(
  log10P ~ pheno_cat_xlabel_gwas_size + log(sample_size),
  data = compare_types
)
emm_lmm <- emmeans(fit_lmm, ~ pheno_cat_xlabel_gwas_size)
ppoint_gwas_stats <- pairs(emm_lmm, adjust = "BH") %>% as.data.frame() %>% 
  mutate(g1 = str_trim(str_split_fixed(contrast, " - ", 2)[, 1]),
         g2 = str_trim(str_split_fixed(contrast, " - ", 2)[, 2]),
         g1 = str_replace_all(g1, "^\\(|\\)$", ""),
         g2 = str_replace_all(g2, "^\\(|\\)$", ""),
         p_str = paste0("p=", signif(p.value, 2))
  ) 

## =====================================
## Panel C: GWAS significance scatter across trait categories
## - x: ordered trait-category groups
## - y: -log10(region-level P)
## - label selected points
## =====================================
ordered <- ordered_xlabels %>% arrange(order) %>% select(pheno_cat_xlabel_gwas) %>% distinct() %>% pull(pheno_cat_xlabel_gwas)
gwas_sug_sig <- unique_regions %>% select(uniqueID,region,region_P,phenotype) %>% distinct() %>% left_join(all_pheno %>% select(phenotype,plot_label,paper_phenotype_category,pheno_cat_xlabel_gwas) )

pd_gwas_sig <- gwas_sug_sig
pd_gwas_sig <- ordered_xlabels %>% select(pheno_cat_xlabel_gwas,paper_phenotype_category,order) %>% 
  mutate(order=if_else(paper_phenotype_category=="Survey question",order+0.5,order)) %>% 
  select(-paper_phenotype_category) %>% group_by(pheno_cat_xlabel_gwas) %>% summarize(order=min(order)) %>% 
  ungroup() %>% right_join(pd_gwas_sig)

# x-axis lookup table (numeric x positions) + the displayed label text
xlabels <- pd_gwas_sig %>% select(order,pheno_cat_xlabel_gwas) %>% distinct() 
xlabels <- xlabels %>% left_join(all_pheno %>% select(pheno_cat_xlabel_gwas,pheno_cat_xlabel_gwas_size) %>% distinct())
xlabels <- xlabels %>% mutate(label=pheno_cat_xlabel_gwas_size)

# Choose a small set of region hits to label, plus best survey hit (excluding body weight question)
labelled_gwas <- pd_gwas_sig %>% ungroup() %>% arrange(uniqueID) %>% filter(pheno_cat_xlabel_gwas!="Body\nweight\nquestion") %>% slice_head(n=10)
top_question <- pd_gwas_sig %>% filter(paper_phenotype_category=="Survey question"&phenotype!="dd_weight_lbs") %>% group_by(paper_phenotype_category) %>% summarize(region_P=min(region_P))  %>% inner_join(pd_gwas_sig)
top_question <- top_question %>% mutate(plot_label=str_replace_all(plot_label,"\n"," ")) %>% 
  mutate(plot_label=str_replace(plot_label,"familiar ","familar\n")) %>% 
  mutate(plot_label=str_replace_all(plot_label,"while ","while\n")) %>%
  mutate(plot_label=str_replace_all(plot_label," with a","\nwith a"))

labelled_gwas <- labelled_gwas %>% bind_rows(top_question)

# Manual line breaks for selected label strings
labelled_gwas <- labelled_gwas %>% mutate(plot_label=str_replace(plot_label," ","
"))
labelled_gwas <- labelled_gwas %>% mutate(plot_label=str_replace(plot_label,"droxy","droxy
"))
labelled_gwas <- labelled_gwas %>% mutate(plot_label=str_replace(plot_label,"Acetic","\nAcetic"))
labelled_gwas <- labelled_gwas %>% mutate(plot_label=str_replace(plot_label,"Lactate","\nLactate"))

# Map EMM contrasts to numeric x coordinates for optional bracket annotation
p_annot <- ppoint_gwas_stats %>% select(g1,g2,p.value,p_str) %>%
  left_join(xlabels %>% transmute(g1 = label, x1 = order), by = "g1") %>%
  left_join(xlabels %>% transmute(g2 = label, x2 = order), by = "g2") %>%
  filter(!is.na(x1) & !is.na(x2)) %>%
mutate(x_min = pmin(x1, x2),x_max = pmax(x1, x2)) 

# Mark which points will be labeled (so we can draw them with non-jittered points)
pd_gwas_sig <- pd_gwas_sig %>% 
  left_join(labelled_gwas %>% select(uniqueID,region,phenotype) %>% distinct() %>% mutate(labeled=TRUE)) %>% 
  replace_na(list(labeled=FALSE))

# Precompute y-positions for bracket overlay (disabled by default below)
y_annot_top <- -log10(min(pd_gwas_sig$region_P))
y_annot_step <- (diff(range(-log10(pd_gwas_sig$region_P))) * 0.04)
p_annot <- p_annot %>% arrange(-x_min, -x_max) %>%
  mutate(
    y = y_annot_top + row_number() * y_annot_step,
    x_mid = (x_min + x_max) / 2
  )

# Top of y-axis rounded to nearest 10 (keeps consistent tick layout)
y_top <- ceiling(max(-log10(pd_gwas_sig$region_P))/10)*10

p <- ggplot(pd_gwas_sig,aes(x=order,y=-log10(region_P))) 
p <- p + geom_jitter(aes(color=paper_phenotype_category),shape=16,data=pd_gwas_sig %>% filter(!labeled),size=1,alpha=0.5,width=0.15)
p <- p + geom_point(aes(color=paper_phenotype_category),shape=16,data=pd_gwas_sig %>% filter(labeled),size=1,alpha=0.5)
p <- p + geom_hline(yintercept=-log10(5e-8),color="grey50",linetype=2,linewidth=0.25)
p <- p + geom_text_repel(aes(label = plot_label,color=paper_phenotype_category),lineheight=0.75,
                         segment.color = "grey40",data=labelled_gwas,hjust = 0, vjust=0.5, nudge_x = 0.2,direction = "y",   
                         force = 5,    segment.size = 0.25,    min.segment.length = 0,max.overlaps =50 , size=1.75)

p <- p + scale_y_continuous(expression(-log[10]*p),
                            limits = c(-log10(1e-6),y_top),
                            expand = expansion(mult = c(0.01, 0.05)),
                            breaks = c(0:4)*20,
                            labels = c(0:4)*20) 

p <- p + scale_x_continuous("", breaks=xlabels$order,labels=xlabels$label,expand = expansion(mult = c(0.1,0.2)))

p <- p + scale_color_manual(values = c("Clinical analytes" = "#f26c2a","Plasma metabolites" = "#2a8c7a","Survey question"= "#e78ac3"))
p <- p + theme_boxes() 

ppoint_gwas <-p 

p

## =====================================
## Panel D: Fraction of dog GWAS regions overlapping human GWAS peaks
## - compute counts by p-value cutoff and match type
## - plot fraction of overlapped regions (two match stringencies)
## =====================================

# Define p-value cutoffs (used as discrete x-axis)
cut_vals <- c(1e-6, 5e-8, 1e-10, 1e-15, 1e-20)

cutoffs <- tibble(cut = cut_vals) %>%
  arrange(desc(cut)) %>%
  mutate(
    order = row_number(),
    # stable key for joins / factor levels
    label_key = paste0("<", format(cut, scientific = TRUE))
  )

# Assemble region list for intermediate phenotypes (kept for reproducibility; not plotted directly)
all_regions <- raw_S3_GWAS_REGIONS %>%
  mutate(phenotype = tolower(phenotype)) %>%
  select(phenotype, region, region_P) %>%
  inner_join(
    all_pheno %>%
      filter(survey_vs_blood == "Intermediate phenotype") %>%
      select(phenotype) %>%
      distinct(),
    by = "phenotype"
  ) %>%
  distinct()

# Overlap calls (human peak overlaps with optional "lifted over" filtering)
overlap <- raw_S4_GWAS_OVERLAP %>% filter(lifted_over)

# Count overlaps passing two match-score thresholds (stricter <=2 and looser <=3)
pd <- overlap %>%
  filter(overlap_human_peak, !is.na(match_score), match_score <= 2) %>%
  mutate(match_type = "same entity or directly related")

pd <- overlap %>%
  filter(overlap_human_peak, !is.na(match_score), match_score <= 3) %>%
  mutate(match_type = "functionally related") %>%
  bind_rows(pd)

pd <- pd %>% select(phenotype, region, region_P, match_type)

# For each cutoff, count regions overlapped by match type
pd <- crossing(pd, cutoffs) %>%
  filter(region_P <= cut) %>%
  group_by(match_type, label_key, cut) %>%
  count()

# Totals at each cutoff (denominator): regions and traits represented
tots <- crossing(overlap %>% select(phenotype, region, region_P) %>% distinct(), cutoffs) %>%
  filter(region_P <= cut) %>%
  group_by(label_key, cut) %>%
  count() %>%
  rename(nregions = n)

tots <- crossing(overlap %>% select(phenotype, region_P) %>% distinct(), cutoffs) %>%
  filter(region_P <= cut) %>%
  select(-region_P) %>%
  distinct() %>%
  group_by(label_key, cut) %>%
  count() %>%
  rename(ntraits = n) %>%
  full_join(tots, by = c("label_key", "cut"))

pd <- pd %>%
  left_join(tots, by = c("label_key", "cut")) %>%
  mutate(frac = n / nregions)

# Factor order for x-axis (largest p-value on left, as in manuscript)
label_levels <- cutoffs %>% arrange(order) %>% pull(label_key)
pd$label_key <- factor(pd$label_key, levels = label_levels)

# X-axis labels as plotmath expressions: p<10^{-6}, p<5 x 10^{-8}, etc.
xlabels <- cutoffs %>%
  mutate(
    xlabel = purrr::map(cut, function(x) {
      exp10  <- as.integer(floor(log10(x)))
      coeff <- x / 10^exp10
      
      if (abs(coeff - 1) < 1e-12) {
        # <10^{-6}  (no 1x)
        bquote("p<" * 10^.(exp10))
      } else {
        # <5x10^{-8}
        bquote("p<" * .(coeff) * x * 10^.(exp10))
      }
    })
  )

# In-panel annotation with denominators at each cutoff (regions and traits)
annot <- cutoffs %>%
  left_join(
    pd %>%
      ungroup() %>%
      select(label_key, nregions, ntraits) %>%
      distinct(),
    by = "label_key"
  ) %>%
  mutate(
    y = 0.02,
    ann_label = paste0("N=", nregions, "\n", ntraits, " traits")
  )

p <- ggplot(pd, aes(x = label_key, y = frac)) +
  geom_line(aes(linetype = match_type, group = match_type, color = match_type), linewidth = 0.3)
p <- p + geom_point(aes(shape = match_type, color = match_type), size = 1.25)

p <- p + geom_text(
    data = annot,
    aes(x = label_key, y = y, label = ann_label),
    inherit.aes = FALSE,
    size = 1.35,
    vjust = 0,hjust=0.5,
    lineheight = 0.9
  )
p <- p + scale_y_continuous(
    "fraction overlap. human GWAS",
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
    expand = expansion(mult = c(0, 0.02))
  )
p <- p + scale_x_discrete(
    "dog GWAS p-value cutoff",
    breaks = cutoffs$label_key,
    labels = xlabels$xlabel,
    expand = expansion(mult = c(0.1, 0.1))
  ) 
p <- p + scale_color_manual(
    values = c(
      "functionally related" = "grey40",
      "same entity or directly related" = "#cb181d"
    ),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) 
p <- p + scale_shape_manual(
    values = c(
      "functionally related" = 17,
      "same entity or directly related" = 16
    ),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) 
p <- p + scale_linetype_manual(
    values = c(
      "functionally related" = 1,
      "same entity or directly related" = 1
    ),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) 
p <- p +  theme_line()

p
p_overlap <- p
p

## =====================================
## Panel B: Heritability vs top GWAS association strength
## - compute Spearman rho + permutation p-value
## - scatter by trait group
## =====================================
pd <- her %>% select(phenotype,paper_phenotype_category,her,min_gwas_pvalue) %>% distinct() %>%
  mutate(group=if_else(phenotype=="dd_weight_lbs","Body weight question",paper_phenotype_category))
pd <- pd %>% mutate(group=if_else(group=="Survey question","Other survey questions",group))
pd$group <- factor(pd$group,levels=c("Clinical analytes","Plasma metabolites","Body weight question","Other survey questions"))
pd <- pd %>% select(phenotype,paper_phenotype_category,group,her,best_p=min_gwas_pvalue) %>% distinct() 

# Observed Spearman correlation between heritability and -log10(best p)
obs_rho <- cor(pd$her,
               -log10(pd$best_p),
               method = "spearman",
               use = "complete.obs")

B <- 100000

# Cache permutation null distribution vector
perm_stats_her_vs_topSNP <- NULL

if (!exists("perm_stats_her_vs_topSNP", inherits = FALSE) ||
    length(perm_stats_her_vs_topSNP) != B) {
  
  set.seed(1)
  
  x <- pd$her
  y <- -log10(pd$best_p)
  
  perm_stats_her_vs_topSNP <- replicate(B, {
    cor(x, sample(y), method = "spearman", use = "complete.obs")
  })
}

# Two-sided permutation p-value for rho
b <- sum(abs(perm_stats_her_vs_topSNP) >= abs(obs_rho))

p_perm <- (b + 1) / (B + 1)
p_min  <- 1 / (B + 1)

# Format p-value text (cap at minimum possible value when b==0)
p_txt <- if (b == 0) {
  paste0("≤", formatC(p_min, format = "fg", digits = 1))
} else {
  paste0("=", formatC(p_perm, format = "fg", digits = 1))
}

rho_txt <- sprintf("Spearman rho = %.2f (p%s)", obs_rho, p_txt)

# Optional regression (printed to console) adjusting for sample size
df <- pd %>% select(phenotype,her,best_p) %>% 
  inner_join(all_pheno %>% select(phenotype,sample_size) %>% distinct()) %>%
  mutate(logp = -log10(best_p))
model <- lm(logp ~ her + sample_size, data = df)
summary(model)

p <- ggplot(pd,aes(y=her,x=-log10(best_p)))
p <- p + geom_point(aes(color=group,shape=group),alpha=0.8,size=1)
p <- p + ggtitle(rho_txt)
p <- p + geom_vline(xintercept=-log10(5e-8),color="grey50",linetype=2,linewidth=0.25)
p <- p + scale_y_continuous(expression("heritability (" * h[SNP]^2 * ")"), limits = c(0, 1), 
                            breaks=c(0:4)/4,expand = expansion(mult = c(0, 0)))
p <- p + scale_x_continuous(expression(-log[10](p) * " (top SNP)"),
                            expand = expansion(mult = c(0.05, 0.05)),
                            breaks = c(0:4)*20)

p <- p + scale_color_manual(values = c(
  "Clinical analytes"    = "#f26c2a",
  "Plasma metabolites"   = "#66c2a5",
  "Body weight question" = "#e78ac3",
  "Other survey questions" = "#e78ac3"
)) 
p <- p + scale_shape_manual(values = c(
  "Clinical analytes"    = 21,
  "Plasma metabolites"   = 24,
  "Other survey questions" = 22,
  "Body weight question" =15
)) 
p <- p + theme_dots()
p_gwas_v_her <- p
p

## =====================================
## Panel E/F: Example Manhattan plots and multi-panel assembly
## =====================================
p_Cys <- pdf_to_plot("source_figures/Cystathionine.manhattan.pdf")
p_Thr <- pdf_to_plot("source_figures/Threonine.manhattan.pdf")
row3 <- plot_grid(
  NULL, p_Thr, p_Cys,
  nrow = 1,
  labels = c("E", "", ""),      # put label on spacer
  label_size = 8,
  label_x = 0.5,                # center of spacer
  label_y = 0.98,
  hjust = 0.5,
  vjust = 1,
  rel_widths = c(0.035, 1, 1.012),  # spacer wide enough for label
  align = "v",
  axis  = "b"
)
row3 <- ggdraw(row3) +
  draw_plot_label(
    label = "F",
    x = 0.465, y = 0.98, hjust = 0, vjust = 1,
    size = 8
  )

# Assemble final multi-panel figure grid and save to PDF
row1 <- plot_grid(pbox_her_blood,pbox_her_survey,p_gwas_v_her,nrow=1,label_size=8,labels=c("A","","B"),rel_widths = c(1,0.8,1.2))
row2 <- plot_grid(ppoint_gwas,p_overlap,nrow=1,ncol=2,label_size=8,labels=c("C","D"),rel_widths = c(1,0.7))
grid <- plot_grid(row1,row2,row3,ncol=1,rel_heights = c(1,0.9,0.6))

ggsave(plot=grid,filename=outpdf,width=5.5,height=6)
