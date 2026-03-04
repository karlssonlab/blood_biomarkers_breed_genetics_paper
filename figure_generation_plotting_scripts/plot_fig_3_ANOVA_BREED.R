# ============================================================
# Figure 3: ANOVA effect sizes (GES), heritability, and dog–human comparisons
# ============================================================

# Load required libraries (only those used in this script)
library(tidyverse)      # Data manipulation (dplyr/tidyr), plotting (ggplot2), and string helpers
library(cowplot)        # Multi-panel plot assembly
library(ggrepel)        # Non-overlapping text labels in scatter plots
library(ggpubr)         # stat_pvalue_manual for significance brackets
library(lme4)           # Mixed-effects models
library(emmeans)        # Estimated marginal means and contrasts
source("util.cleaning.R")  # Project utilities (e.g., replace_middle_break)

# Output path for the final assembled figure
outpdf <- "../../../figures_pdf/fig_3_ANOVA_BREED.pdf"


## =====================================
## Load all data (from precomputed RDS)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # expose tables like raw_S5_ANOVA, raw_S2_HERITABILITY, etc.
}

# Subset ANOVA results for two models: with breed and without breed
anova_all_xbreed <- raw_S5_ANOVA  %>% mutate(Effect=str_remove(Effect," status")) %>% filter(model=="age + weight + sex status")
anova_all <- raw_S5_ANOVA  %>% mutate(Effect=str_remove(Effect," status")) %>% filter(model== "breed + age + weight + sex status")
 
# ============================================================
# Plot theming helpers (shared formatting across panels)
# ============================================================
theme_box <- function() {
  theme_cowplot(8) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text       = element_text(size = 6),
      axis.title.y    = element_text(size = 6),
      axis.title.x    = element_blank(),
      strip.text = element_text(size = 6, hjust = 0.5, vjust = 0,face = "bold",
        margin = margin(t = 2, b = 5),color = "black"),
      strip.background = element_rect(fill = "white", color = NA),
      
      panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text            = element_text(size = 5),
      legend.title           = element_blank(),
      legend.key.size        = unit(0.15, "cm"),
      legend.margin          = margin(1,1,1,1),
      legend.box.background  = element_rect(fill = "white", color = "grey50", linewidth = 0.15)
    )
}

theme_stacked <- function() {
  theme_box() %+replace%
    theme(
      plot.title       = element_text(size = 5.5, face = "bold", hjust = 0),
      axis.ticks.x     = element_line(linewidth = 0.25),
      axis.line.x      = element_line(linewidth = 0.25),
      axis.title.y     = element_blank(),
      axis.title.x     = element_text(size = 6),
      axis.text.y     = element_text(size = 5,hjust=1,vjust=0.5,margin=margin(l=2,r=3)),
      strip.text = element_text(size = 5.5,face="bold",hjust=0.5,vjust=0.5,margin=margin(l=2,r=3)),
      legend.position  = c(0.98, 0.02),
      legend.justification = c("right", "bottom"),
      legend.text      = element_text(size = 4.5),
      legend.key.spacing.y = unit(0.05,"cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.2),
      plot.margin      = margin(t = 5, r = 5, l = 5, b = 0),
    )
}

theme_stacked_comp <- function() {
  theme_stacked() %+replace%
    theme(
     legend.position  = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
     axis.ticks.x     = element_line(linewidth = 0.25),
     axis.ticks.y     = element_blank(),
     axis.line.x      = element_line(linewidth = 0.25),
     axis.title.y     = element_blank(),
     axis.title.x     = element_text(size = 6),
     
    )
}

theme_scatter <- function() {
  theme_box() %+replace%
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      axis.line.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.title.x = element_text(size = 6),
      panel.grid.major.x = element_line(color = "grey70", linewidth = 0.25),
      axis.line.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25),
      
    )
}

# ============================================================
# Prepare analysis table: merge ANOVA effect sizes with heritability
# ============================================================

# Extract heritability estimates and keep only complete entries
her <- raw_S2_HERITABILITY
her <- her %>% select(phenotype,GREML_LDS_constrained_heritability,GREML_LDS_constrained_stderr)
her <- her  %>% filter(!is.na(GREML_LDS_constrained_heritability)) %>% rename(her=GREML_LDS_constrained_heritability,her_err=GREML_LDS_constrained_stderr)

# Join phenotype metadata, ANOVA results, and heritability estimates
aov <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label) %>% 
  right_join(anova_all)
aov <- aov %>% select(phenotype,paper_phenotype_category,plot_label,Effect,ges,pFDR)
aov <- aov %>% left_join(her)

# Set plotting/analysis factors
pd <- aov
pd$paper_phenotype_category <- factor(pd$paper_phenotype_category,levels=c("Plasma metabolites","Clinical analytes"))
pd$Effect <- factor(pd$Effect,levels=c("breed","age", "sex","weight"))

# ============================================================
# Mixed model: compare GES across effects while accounting for phenotype
# ============================================================
fit_ges <- lmer(ges ~ Effect + (1 | phenotype), data = pd)
emm_ges <- emmeans(fit_ges, ~ Effect)
contr_ges <- contrast(emm_ges, method = "pairwise", adjust = "tukey")
contr_breed <- contrast(emm_ges, method = list(
  "breed vs age"    = c( 1, -1,  0,  0),
  "breed vs sex"    = c( 1,  0, -1,  0),
  "breed vs weight" = c( 1,  0,  0, -1)
))
cb_df <- as.data.frame(contr_breed)

# ============================================================
# Compare model variants (with vs without breed) summary stats
# ============================================================
df <- anova_all %>% select(-`p<.05`,-model) %>% 
  drop_na() %>% mutate(model="with_breed")
df <- anova_all_xbreed %>% select(all_of(intersect(names(anova_all_xbreed), names(df)))) %>% 
  mutate(model="without_breed") %>% bind_rows(df)
df_w_2models <- df
stats <- df %>% group_by(Effect,model) %>% 
  summarize(mean_ges = mean(ges, na.rm = TRUE),sd_ges=sd(ges, na.rm = TRUE),n= n())


raw_anova_xbreed_selected <- anova_all_xbreed

# Build annotation table for plotting breed-vs-others significance brackets
y_max <- max(pd$ges, na.rm = TRUE)

stat_df <- cb_df %>%
  mutate(
    group1 = "breed",
    group2 = case_when(
      grepl("age", contrast)    ~ "age",
      grepl("sex", contrast)    ~ "sex",
      grepl("weight", contrast) ~ "weight"
    ),
    p         = p.value,
    label     = sprintf("p = %.2g", p),
    y.position = y_max + c(0.1, 0.06, 0.02)
  ) %>%
  select(group1, group2, p, label, y.position)

# ============================================================
# Panel A: Boxplots of variance explained (GES) by effect
# ============================================================
y_max <- max(pd$ges, na.rm = TRUE)
p <- ggplot(pd,aes(x=Effect,y=ges)) 
p <- p + geom_boxplot(aes(fill=paper_phenotype_category,color=paper_phenotype_category),outlier.size=0.25,linewidth=0.25,alpha=0.25,outlier.alpha=0.75)
p <- p + stat_pvalue_manual(data=stat_df,label = "label",
                            bracket.size = 0.2,tip.length = 0.01, size = 1.5)
p <- p + scale_x_discrete(expand = expansion(mult = c(0.04,0.02)))
p <- p + scale_y_continuous("variance explained",limits=c(0,0.55))
p <- p + scale_color_manual(name="category",values = c(
  "Plasma metabolites"="#2D8070","Clinical analytes"="#E66100"))

p <- p + scale_fill_manual(name="category",values = c(
  "Plasma metabolites"="#2D8070","Clinical analytes"="#E66100"))
p <- p + theme_box()
p_box <- p
p

# ============================================================
# Identify phenotypes present in human aging models (for labeling/joins)
# ============================================================
in_aging_models <- raw_S4_VS_HUMAN %>% select(phenotype,`Models_AgingAI`|`Models_phenoAge`) %>% pivot_longer(-phenotype)
in_aging_models <- in_aging_models %>% filter(value) %>% select(-value) %>% mutate(name=str_remove(name,"Models_"))
in_aging_models <- in_aging_models %>% arrange(phenotype,name)
in_aging_models <- in_aging_models %>% group_by(phenotype) %>% summarize(models=paste(name,collapse="; ")) %>% mutate(models=str_remove(models,"; others$"))

# ============================================================
# Mixed model: test whether GES–heritability slope differs by effect
# ============================================================
fit_ges_her <- lmer(ges ~ her * Effect + (1 | phenotype), data = pd)
trends <- emtrends(fit_ges_her, ~ Effect, var = "her")
trends <- as.data.frame(summary(trends, infer = TRUE))
stat_text <- trends %>% filter(Effect=="breed") %>% mutate(str=paste0("slope=",round(her.trend,3),"\np=",formatC(p.value, format = "e", digits = 1))) %>% pull(str)
print(trends)

# ============================================================
# Panel B: Breed GES vs heritability scatter (with FDR significance)
# ============================================================
pd <- aov %>% filter(Effect=="breed") %>% drop_na() %>% mutate(sig=if_else(pFDR<=0.05,TRUE,FALSE))
pd <- pd %>% mutate(rank=min_rank(desc(ges))) 
labeled <- pd %>% inner_join(in_aging_models) %>% mutate(plot_label=map_chr(plot_label, replace_middle_break))
maxy <- ceiling(max(pd$her)*100)/100
maxx <- ceiling(max(pd$ges)*100)/100

p <- ggplot(pd,aes(y=her,x=ges))
p <-p + geom_point(data=pd %>% filter(!sig),color="grey30",shape=16,size=0.75,alpha=0.5)

p <-p + geom_point(data=pd %>% filter(sig),color="#3969AC",shape=16,size=1,alpha=0.75)
p <- p + annotate("text",label=stat_text,x=min(pd$ges),y=max(pd$her),size=1.75,vjust=1,hjust=0)
p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "grey30") 
p <- p + scale_y_continuous( "heritability",limits=c(min(pd$her),max(pd$her)),
                            breaks = seq(0, 1, by = 0.25),
                            expand = expansion(mult = c(0.05, 0.05))    ) 

p <- p + scale_x_continuous( "variance explained by breed",
                             breaks = seq(0, 0.4, by = 0.05),
                             expand = expansion(mult = c(0.05, 0.01))    ) 
p <- p + theme_scatter() + theme(plot.margin = margin(l=10,r=10,t=5,b=5))
p_her <- p
p

# ============================================================
# Panel D: Breed vs other effects (per phenotype) scatter
# ============================================================
df_breed <- aov %>% filter(Effect=="breed") %>% select(phenotype,ges.breed=ges,pFDR.breed=pFDR)
pd <- aov %>% filter(Effect!="breed") %>% inner_join(df_breed) %>% select(-her,-her_err) %>% distinct()
pd <- pd %>% mutate(sig=if_else(pFDR<=0.05,TRUE,FALSE))
maxy <- ceiling(max(pd$ges)*100)/100
maxx <- ceiling(max(pd$ges.breed)*100)/100
miny <- floor(min(pd$ges)*100)/100
minx <- floor(min(pd$ges.breed)*100)/100
pd <- pd %>% mutate(rank=min_rank(desc(ges))) 
labeled <- pd %>% filter(rank<=10) %>% 
  mutate(plot_label=map_chr(plot_label, replace_middle_break))

p <- ggplot(pd,aes(x=ges.breed,y=ges))
p <-p + geom_point(data=pd %>% filter(!sig),shape=16,color="grey40",size=0.5,alpha=0.25)
p <-p + geom_point(data=pd %>% filter(sig),aes(color=Effect,shape=Effect),size=1,alpha=0.75)
p <- p + geom_text_repel(
  aes(label = plot_label,color=Effect),
  lineheight = 0.9, segment.color = "grey40",
  data = labeled, hjust = 0, vjust = 1, nudge_x = 0.005, direction = "y",
  force = 3, segment.size = 0.2, min.segment.length = 0, max.overlaps = 50, size = 1.25)
p <- p + scale_y_continuous( "variance explained by other effects",
                             limits=c(miny,maxy*1.1),
                             breaks = seq(0, 0.4, by = 0.05),
                             expand = expansion(mult = c(0.01, 0.01))    ) 
p <- p + scale_x_continuous( "variance explained by breed",
                             limits=c(minx,maxx),
                             breaks = seq(0, 0.4, by = 0.05),
                             expand = expansion(mult = c(0.05, 0.05))    ) 

p <- p + scale_color_manual(values = c(
  "age" = "#A50026",
  "weight" = "#E66100",
  "sex" = "#CC79A7"
)) 
p <- p + theme_scatter() + theme(axis.line = element_line(color = "black", linewidth = 0.25))

p_scatter <- p

# ============================================================
# Panel C: Effect attenuation example (weight): with vs without breed model
# ============================================================
df <- df_w_2models
df <- df %>% filter(Effect=="weight") %>% 
  select(phenotype, model, ges, pFDR) %>%
  pivot_wider(names_from  = model,values_from = c(ges, pFDR),names_sep   = "__")
df <- df %>%
  mutate(ratio_ges = ges__with_breed / ges__without_breed,
         attenuation = 1 - ratio_ges)

df <- df %>% mutate(flagged=if_else(pFDR__without_breed < 0.05,"Sig. (-breed)","Not sig. (-breed)"))
df <- df %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label))
df <- df %>% mutate(sig=if_else(pFDR__with_breed<=0.05,"Sig. (+breed)","Not sig. (+breed)"))
xlabel = "variance exp. by weight (-breed)"
ylabel = "var. exp. by weight (+breed)"
ct <- cor.test(
  df$ges__without_breed,
  df$ges__with_breed,
  method = "spearman",
  use = "pairwise.complete.obs"
)
p <- ggplot(df,aes(x = ges__without_breed,y = ges__with_breed))
p <- p + geom_abline(slope = 1, intercept = 0,linetype = 2, color = "grey60")
p <-p + geom_point(data=df,aes(color=sig,shape=sig,size=flagged),alpha=0.65)
p <-p + annotate("text",x = 0.01, y = 0.1,
                 label = sprintf("Spearman rho = %.2f\np = %.2g",unname(ct$estimate),ct$p.value),
                 hjust = 0, vjust = 1.1,size = 1.75,face="bold")

p <- p + geom_text_repel( aes(label = plot_label),color="#E66100",lineheight = 0.9, segment.color = "grey40",
                          data = df %>% filter(sig=="Sig. (+breed)"), hjust = 0, vjust = 1, nudge_x = 0.005, direction = "y",
                          force = 4, segment.size = 0.2, min.segment.length = 0, max.overlaps = 50, size = 1.5)
p <- p + scale_color_manual(
  values = c("Not sig. (+breed)" = "grey40", "Sig. (+breed)" = "#E66100"),
  name = "Sig. (breed in model)")
p <- p + scale_shape_manual(
  values = c("Not sig. (+breed)" = 16, "Sig. (+breed)" = 15),
  name = "Sig. (breed in model)")
p <- p + scale_size_manual(
  values = c("Not sig. (-breed)" = 0.75, "Sig. (-breed)" = 1.25),
  name = "Sig. (w/o breed in model)")
p <- p + scale_x_continuous(xlabel,limits=c(0,0.15),breaks=c(0,0.05,0.10,0.15),expand = expansion(mult = c(0, 0)))
p <- p + scale_y_continuous(ylabel,limits=c(0,0.1),breaks=c(0,0.05,0.10,0.15),expand = expansion(mult = c(0, 0.02)))
p <- p + theme_scatter() + theme(legend.title=element_blank(),
                                 legend.position = c(0.98, 0.98),
                                 legend.justification = c("right", "top"),
                                 legend.box.spacing = unit(2, "pt"),
                                 legend.spacing.y = unit(2, "pt"),
                                 legend.spacing.x = unit(2, "pt"),
                                  plot.margin = margin(l=10,r=10,t=5,b=5))
                                 
                                 
p_wo_breed <- p 
p

# ============================================================
# Panel E: Top-N phenotypes by breed GES with stacked contributions
# ============================================================

# Keep only significant phenotype/effect pairs; then shorten some long labels
df <- aov %>% select(-her,-her_err) %>% filter(pFDR <= 0.05)

df <- df %>% mutate(
      plot_label = str_replace(plot_label, "corpuscular", "corp."),
      plot_label = str_replace(plot_label, "distribution", "dist."),
      plot_label = str_replace(plot_label, "relative", "rel."),
      plot_label = str_replace(plot_label, "absolute", "abs."),
      plot_label = str_replace(plot_label, "transferase", "trans."))
empty <- crossing(raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label),tibble(Effect=c("age","breed","weight","sex")))
                                             
df <- df %>% mutate(Effect = factor(Effect, levels = rev(c("breed", "age", "sex", "weight"))))
  
# Compute total explained variance per phenotype for ordering/annotation
sums <- df %>% group_by(phenotype) %>%
    summarise(total_ges = sum(ges, na.rm = TRUE), .groups = "drop")
df <- df %>% inner_join(sums)

# Rank effects within each Effect group, then keep phenotypes with highest ranks across effects
df <- df %>% group_by(Effect) %>% mutate(rank=min_rank(desc(ges))) 
df <- df %>% group_by(phenotype) %>% summarize(min_rank=min(rank)) %>% inner_join(df)

# Human association table (hazard ratio -> signed log2 scale)
pd_human <- raw_S4_VS_HUMAN %>% select(phenotype,PMID,human_hazard_ratio) %>% distinct() %>% drop_na()
pd_human <- pd_human %>% mutate(direc=if_else(human_hazard_ratio<1,"protect","risk"))
pd_human <- pd_human %>% mutate(logHR=if_else(human_hazard_ratio<1,log2(1/human_hazard_ratio),log2(human_hazard_ratio)))
pd_human <- pd_human %>% filter(phenotype %in% df$phenotype) %>% arrange(desc(logHR)) %>% mutate(rank=row_number()) 

plot_in <- df %>% left_join(in_aging_models) 

# Restrict to top 20 by breed explained variance, then keep all effects for those phenotypes
plot_in <- plot_in  %>% ungroup() %>% filter(Effect=="breed") %>% select(phenotype,ges) %>% 
  distinct() %>% mutate(breed_ges_rank=rank(-ges,ties.method="min")) %>% 
  filter(breed_ges_rank<=20) %>% select(phenotype,breed_ges_rank) %>% 
  inner_join(plot_in)

# Reorder y-axis labels by descending breed rank (for the horizontal bar plot)
order_levels <- plot_in %>% ungroup() %>% select(plot_label,breed_ges_rank) %>% distinct() %>% arrange(-breed_ges_rank) %>%
    pull(plot_label)
  
plot_in <- plot_in %>% mutate(plot_label = factor(plot_label, levels = order_levels))
  
# Compute max total (not used directly in scale limits here, but kept for parity)
maxy <- plot_in %>%
  group_by(phenotype) %>%
  summarise(tot = sum(ges), .groups = "drop") %>%
  pull(tot) %>%
  max()
maxy <- maxy+0.05

plot_in$Effect <- factor(plot_in$Effect, levels = rev(c("breed", "age", "sex", "weight")))
plot_in <- plot_in %>% mutate(paper_phenotype_category=str_replace(paper_phenotype_category," ","\n"))

p <- ggplot(plot_in, aes(x = plot_label, y = ges))
p <- p + geom_bar(aes(fill = Effect),color = "grey30",linewidth = 0.1,stat = "identity",width = 0.5)
p <- p +   coord_flip() +
  scale_y_continuous(
    "variance explained",
    breaks = c(0, 0.2, 0.4, 0.6),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 0.55)
  ) +
  scale_fill_manual(values = c(
    "age"    = "#A50026",
    "breed"  = "#7AA0D4",
    "weight" = "#E66100",
    "sex"    = "#CC79A7"
  )) 
p <- p + facet_grid(
  paper_phenotype_category ~ .,
  scales = "free_y",
  space  = "free_y",
  switch = "y"
)
p <- p + theme_stacked() +
  theme(
    strip.placement = "outside"
  ) 
p_topN  <- p 
p

# ============================================================
# Panels F: Human hazard ratio vs dog explained variance (paired layout)
# ============================================================

# Phenotypes included in human aging models (set used for dog+human comparison)
inc <- in_aging_models %>% select(phenotype) %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype)) %>% distinct()

# Build long-form table with: human (risk/protect logHR) and dog (ANOVA GES by effect)
pd <- inc %>% inner_join(pd_human %>% filter(!is.na(human_hazard_ratio)&human_hazard_ratio > 0) %>% select(phenotype,direc,logHR) %>% 
                         rename(Effect=direc,logHR_human=logHR) %>% 
                         mutate(species="human",logHR_human=-logHR_human) %>% distinct())

pd <- inc %>% inner_join(anova_all %>% filter(pFDR<=0.05) %>% select(phenotype,Effect,ges) %>%  distinct()) %>% mutate(species="dog") %>% bind_rows(pd) 
pd <- pd %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype,plot_label) %>% distinct())
pd <- pd %>% mutate(Effect = factor(Effect, levels = rev(c("breed", "age", "sex", "weight","risk","protect"))))

# Add ranks to annotate human panel with dog breed-rank for the same traits
dog_rank_for_printing <- pd %>% filter(species=="dog"&Effect=="breed") %>% select(phenotype,ges) %>% mutate(dog_rank_n=rank(-ges,ties.method = "first")) %>%
   mutate(dog_rank=if_else(dog_rank_n==2,paste0("breed effect rank in dog = ",dog_rank_n),as.character(dog_rank_n))) %>% select(-ges)
human_rank_for_printing <- pd %>% filter(species=="human") %>% select(phenotype,plot_label,logHR_human) %>% mutate(human_rank=rank(logHR_human,ties.method = "first")) %>% select(-logHR_human)

pd <- pd %>% left_join(dog_rank_for_printing) %>% inner_join(human_rank_for_printing)

# Enforce a shared y-axis order between the two paired plots (based on human rank)
ordered <- human_rank_for_printing %>% 
  arrange(-human_rank) %>% pull(plot_label)

pd$plot_label <- factor(pd$plot_label,levels=ordered)

# Alternating background guide lines (vertical after coord_flip)
hlines<- seq(from = 0.5, to = length(unique(pd$plot_label)), by = 2)

# Left: human hazard ratios (signed log scale) with dog rank annotation
p <- ggplot(pd %>% filter(species=="human") %>% select(plot_label,logHR_human,dog_rank) %>% distinct(), aes(x = plot_label, y = logHR_human)) 
p <- p + geom_vline(xintercept=hlines,color="grey80",linewidth=0.2)
p <- p + geom_bar(aes(fill = Effect), data=pd %>% filter(species=="human"),color="grey10",linewidth=0.1,stat = "identity", width = 0.5) 
p <- p + geom_text(aes(y=logHR_human-0.05,label=dog_rank),size=1,color="#4A74B5",hjust=1,vjust=0.5)
p <- p + coord_flip() 
p <- p + scale_y_continuous( "human hazard ratio",
                             breaks = -log2(c(1,1.4,2,2.8,4)),
                             labels = c(1,1.4,2,2.8,4),
                             expand = expansion(mult = c(0.05,0))    ) 
p <- p + scale_fill_manual(
  values = c("age" = "#A50026","breed" = "#7AA0D4","weight" = "#E66100","sex" = "#CC79A7","risk" = "black","protect" = "grey50"),
  labels = c(age = "age",breed  = "breed",weight = "weight",sex = "sex", risk = "risk", protect = "protect")
)
p <- p + theme_stacked_comp() + theme(
  legend.position = "none",
  plot.margin      = margin(t = 5, r = 5, l = 10, b = 0),
)
p_left <- p
p

# Right: dog variance explained for the same phenotypes/effects
pd <- pd %>% filter(species!="human") 
pd$plot_label <- factor(pd$plot_label,levels=ordered)

p <- ggplot(pd, aes(x = plot_label, y = ges) )
p <- p + geom_vline(xintercept=hlines,color="grey80",linewidth=0.2)

p <- p + geom_bar(aes(fill = Effect), color="grey10",linewidth=0.1,stat = "identity", width = 0.5) 
p <- p + coord_flip() 
p <- p + scale_y_continuous( "variance explained",
                             breaks = c(0,0.2,0.4,0.6),limits=c(0,0.6),
                             expand = expansion(mult = c(0,0))    ) 
p <- p + scale_fill_manual(
  values = c("age" = "#A50026","breed" = "#7AA0D4","weight" = "#E66100","sex" = "#CC79A7","risk" = "black","protect" = "grey50"),
  labels = c(age = "age",breed  = "breed",weight = "weight",sex = "sex", risk = "risk", protect = "protect")
)
p <- p + theme_stacked_comp() + theme(axis.text.y=element_blank(),
                                      legend.position  = c(0.98, 0.02),
                                      legend.justification = c("right", "bottom"))
                                      
p_right <- p

# ============================================================
# Assemble multi-panel figure and save
# ============================================================
grid_paired <- plot_grid(p_left,p_right,nrow=1,label_size=8,rel_widths=c(1,0.7))

col1 <- plot_grid(p_her,p_wo_breed,ncol=1,rel_heights = c(1,1),label_size=8,labels=c("","C"))

row1 <- plot_grid(p_box,col1,p_scatter,nrow=1,ncol=3,label_size=8,labels=c("A","B","D"),rel_widths=c(0.6,1,1))
row2 <- plot_grid(p_topN,grid_paired,nrow=1,labels=c("E","F"),rel_widths=c(0.8,1),label_size = 8,label_x = 0.01,label_y = 1,hjust = 0)
grid <- plot_grid(row1,row2,nrow=2,label_size=8,rel_heights = c(1.1,1))
ggsave(plot=grid,filename=outpdf,width=6,height=4.5)


# ============================================================
# Post-figure analyses: rank-based dog–human concordance and enrichment tests
# ============================================================

# Build stats table: phenotype metadata + inclusion in human models + human logHR + dog ranks
stats <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category)
stats <- stats %>% left_join(in_aging_models) %>% left_join(pd_human %>% select(phenotype,logHR) %>% distinct())
dog_ranks <- anova_all %>% select(phenotype,Effect,ges) %>% group_by(Effect) %>% mutate(dog_rank=rank(-ges,ties.method="min"))
stats <- stats %>% left_join(dog_ranks)
stats <- stats %>% mutate(in_model=if_else(is.na(models),FALSE,TRUE))

# ============================================================
# 1) Dog–human concordance controlling for category (per Effect)
#    Rank-based: correlate rank(ges) with rank(residualized logHR)
# ============================================================

res_spearman_cov <- stats %>%
  filter(!is.na(logHR), !is.na(ges), !is.na(paper_phenotype_category)) %>%
  group_by(Effect) %>%
  filter(n() >= 5, n_distinct(paper_phenotype_category) >= 1) %>%
  group_modify(~{
    dat <- .x
    
    # Residualize human logHR on phenotype category to remove category-level shifts
    fit_h <- lm(logHR ~ paper_phenotype_category, data = dat)
    r_logHR <- resid(fit_h)
    
    # Spearman correlation between dog effect size and residualized human signal
    ct <- cor.test(rank(dat$ges), rank(r_logHR),
                   method = "spearman",
                   alternative = "less")  # directionality note: adjust if sign convention changes
    
    broom::tidy(ct)
  }) %>%
  ungroup() %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "BH"))

res_spearman_cov


# ============================================================
# 2) Enrichment: human-supported traits have larger ges,
#    controlling for category (per Effect)
#    (lm is the natural covariate-adjusted version)
# ============================================================

res_ges_enrich_cov <- stats %>%
  filter(!is.na(ges), !is.na(paper_phenotype_category)) %>%
  group_by(Effect) %>%
  filter(n_distinct(in_model) == 2) %>%  # need TRUE and FALSE
  group_modify(~{
    fit <- lm(ges ~ in_model + paper_phenotype_category, data = .x)
    broom::tidy(fit)
  }) %>%
  ungroup() %>%
  filter(term == "in_modelTRUE") %>%      # adjusted effect of being in human model(s)
  mutate(p_adj_fdr = p.adjust(p.value, method = "BH"))

res_ges_enrich_cov


# Optional: one-sided p-value for positive enrichment (in_modelTRUE > 0)
res_ges_enrich_cov <- res_ges_enrich_cov %>%
  mutate(
    p_one_sided = if_else(estimate > 0, p.value/2, 1 - p.value/2),
    p_one_sided_fdr = p.adjust(p_one_sided, method = "BH")
  )

res_spearman_cov
res_ges_enrich_cov
