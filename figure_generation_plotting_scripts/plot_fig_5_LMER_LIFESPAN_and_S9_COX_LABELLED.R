# ============================================================
# Libraries and project helpers
# ============================================================
library(tidyverse)      # Data manipulation + ggplot2
library(ggrepel)        # Non-overlapping text labels for plots

# Helper functions used throughout (e.g., name standardization, label formatting)
source("util.cleaning.R")

# ============================================================
# Output paths
# ============================================================
outpdf <- "../../../figures_pdf/fig_5_LMER_LIFESPAN.pdf"
outpdf_labelled <- "../../../figures_pdf/fig_S9_COX_LABELLED.pdf"


# ============================================================
# Load analysis objects (creates raw_S* tables in the global env)
# ============================================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_anc, etc.
}

# Naming map for effects used in joins, ordering, and plot legends
effect_names <- tibble::tibble(
  Effect      = c("breed.weight.kg", "REML.t.val", "lifespan"),
  Effect_name = c("breed weight",  "Breed Ancestry Score", "breed lifespan"),
  order       = c(3, 4, 1)
)

# ============================================================
# Plot themes (shared styling across panels)
# ============================================================
theme_point <- function() {
  theme_cowplot(6) +
    theme(
      axis.line = element_line(linewidth = 0.25),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      
      axis.ticks = element_line(linewidth = 0.5),
      
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.25),
      panel.grid.major = element_blank(),
      
      strip.text = element_text(
        size = 6, hjust = 0, vjust = 0,
        face = "bold",
        margin = margin(t = 2, b = 2),
        color = "black"
      ),
      strip.background = element_rect(fill = "white", color = NA),
      
      legend.position = "top",
      legend.box = "vertical",
      legend.text = element_text(size = 4),
      legend.title = element_text(size = 4, face = "bold"),
      legend.key.size = unit(0.15, "cm"),
      legend.key.spacing.y = unit(0.05, "cm"),
      legend.box.background = element_rect(
        fill = "white", color = "grey50", linewidth = 0.15
      ),
      legend.margin = margin(1, 1, 1, 1)
    )
}

# Theme for compact bar/estimate-style panels
theme_bar_estimate <- function(){
  theme_cowplot(6) %+replace%
    theme(
      plot.title = element_text(size = 5,face="bold",hjust=0),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(linewidth=0.25),
      axis.line.y=element_blank(),
      axis.line.x = element_line(linewidth = 0.25),    
      axis.text = element_text(size = 5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 6),
      legend.position = "none",  
     panel.grid.major.x=element_blank(),
      strip.text = element_text(size = 6, hjust = 0.5, vjust = 0, face = "bold",margin = margin(t = 0, b = 3)),
      strip.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}


# ============================================================
# Breed covariates used in models/joins (standardize + recode)
# ============================================================
breeds <- raw_S12_BREED_PHENOS %>% select(breed,breed.height.cm,breed.weight.kg,curl,furnishings,HR=`HR.(McMillan.2024)`,lifespan=`median_age_death.(McMillan.2024)`) %>% distinct()
breeds <- standardize_breed_names(breeds,breed)
breeds <- convert_yesno_to_logical(breeds)


# ============================================================
# Assemble LMER/ANOVA summary table and attach metadata
# ============================================================
model_to_use <- "REML.t.val ~ breed.weight.kg + lifespan"
pdAll <- raw_S9_ANOVA_ON_LMER %>% filter(model_used==model_to_use) %>% 
  select(phenotype,Effect,ges=ges.anova,p.anova,p_adj_fdr.anova,term.lm,beta_std.lm,p.value.lm, conf.low.lm, conf.high.lm) %>% 
  drop_na()

# Add phenotype labels/categories and map Effect -> Effect_name/order
pdAll <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label) %>% 
  distinct() %>% right_join(pdAll) %>% inner_join(effect_names) 

# Add combined-statistics table (drops model_used since model already filtered)
pdAll <- pdAll %>% left_join(raw_S11_COMBINED %>% select(-model_used))

# Add Cox hazard results and a log(HR) column
pdAll <- pdAll %>% left_join(
  raw_S10_COX_HAZ %>% select(phenotype,exp_coef.cox = `exp(coef)`,p.cox = `Pr(>|z|)`,
                     p_adj_fdr.cox = p_adj_fdr,in_ben = metabolite_significant_in_Harrison_et_al)) %>% 
  mutate(logHR = log(exp_coef.cox))

# Summarize which human aging models include each phenotype (for annotation/flagging)
human <- raw_S4_VS_HUMAN %>% select(phenotype, AgingAI = `Models_AgingAI`, phenoAge = `Models_phenoAge`) %>%
  pivot_longer(-phenotype, names_to = "name", values_to = "value") %>%
  filter(value) %>% arrange(name) %>% group_by(phenotype) %>% summarize(models = paste(name, collapse = "\n"), .groups = "drop")
pdAll <- pdAll %>% left_join(human,by = "phenotype") %>% distinct() %>% 
  filter(Effect !="odd_fur")
dim(pdAll)


# ============================================================
# Plot palettes (subset to only used levels to avoid scale warnings)
# ============================================================
pal <- c(
  "breed weight"           = "#3070ad",
  "blood trait (LMER)"     = "#971d2b",
  "breed mortality hazard" = "#971d2b",
  "breed lifespan"         = "#971d2b"
)
pal.labels <- c(
  "breed weight"           = "weight",
  "blood trait (LMER)"     = "blood trait",
  "breed lifespan"         = "lifespan"
)
pal_use <- pal[names(pal) %in% unique(pdAll$Effect_name)]
pal.labels_use <- pal.labels[names(pal.labels) %in% unique(pdAll$Effect_name)]

# Significant subsets for ordering/filtering
pd_nom_sig <- pdAll %>% ungroup() %>% filter(p.anova <= 0.05) 
pd_fdr_sig <- pdAll %>% ungroup() %>% filter(p_adj_fdr.anova <= 0.05) 

# ============================================================
# Panel A: stacked variance explained (GES) across significant phenotypes
# ============================================================

# Ensure consistent legend order for effects (only those present)
pd_nom_sig$Effect_name <- factor(pd_nom_sig$Effect_name, levels = effect_names %>% filter(Effect_name %in% pd_nom_sig$Effect_name) %>% arrange(desc(order)) %>% pull(Effect_name))

# Compute total variance explained per phenotype and keep per-effect rows
pd_nom_sig <- pd_nom_sig %>% group_by(phenotype) %>% summarize(ges_tot=sum(ges)) %>% full_join(pd_nom_sig)

# Restrict to phenotypes that pass FDR overall
pd <- pd_nom_sig %>% inner_join(pd_fdr_sig %>% select(phenotype) %>% distinct())

# Order phenotypes by total variance explained for plotting
order <- pd %>% select(phenotype, plot_label, ges_tot) %>% distinct() %>% arrange(ges_tot) %>%
  mutate(rank = row_number())

order_from_row1 <-order

pd$plot_label  <- factor(pd$plot_label, levels = order %>% arrange(rank) %>% pull(plot_label))

# Alternating background guide lines (in flipped coordinates this becomes horizontal separators)
maxy <- length(unique(pd$phenotype))
hlines<- seq(from = 0.5, to = maxy, by = 2)

p <- ggplot(pd, aes(y = ges, x = plot_label))
p <- p + geom_vline(xintercept=hlines,color="grey80",linewidth=0.2)
p <- p + geom_bar(aes(fill = Effect_name), stat = "identity", width = 0.5)
p <- p + coord_flip()
p <- p + scale_y_continuous("variance explained",expand = expansion(mult = c(0,0.02)),breaks=c(0,0.2,0.4))
p <- p + scale_x_discrete("",breaks=pd$plot_label,labels=pd$plot_label,expand = expansion(mult = c(0.03,0.03))) #expand = expansion(mult = c(0.08,0.05)))
  
p <- p + scale_fill_manual(values = pal_use,labels=pal.labels_use)
p <- p + theme_bar_estimate() 
p_bar <- p
# Wrap in ggdraw to control panel scaling/alignment in final multi-panel figure
p_bar <- ggdraw() + draw_plot(p_bar, x = 0, y = 0, width = 1, height = 0.943)
p_bar

# ============================================================
# Panel B: standardized effect sizes (LM) with confidence intervals by effect
# ============================================================
pd <- pd %>% select(phenotype,plot_label) %>% distinct() %>% inner_join(pdAll)

# Axis range helpers based on CI limits
x_min <- floor(min(pd$conf.low.lm,  na.rm = TRUE)*10)/10
x_max <- ceiling(max(pd$conf.high.lm, na.rm = TRUE)*10)/10
pd$plot_label  <- factor(pd$plot_label, levels = order %>% arrange(rank) %>% pull(plot_label))

# Facet order for effects
pd$Effect_name <- factor(pd$Effect_name, levels = effect_names %>% arrange(order) %>% pull(Effect_name))

# Alpha encodes nominal significance (de-emphasize non-significant points)
pd <- pd %>% mutate(sig.anova=if_else(p.anova<=0.05,TRUE,FALSE))
maxy <- length(unique(pd$phenotype))

hlines<- seq(from = 0.5, to = maxy, by = 2)
unique(pd$plot_label)

p <- ggplot(pd, aes(y = plot_label, x = beta_std.lm))
p <- p + geom_vline(xintercept=0)
p <- p + geom_hline(yintercept=hlines,color="grey80",linewidth=0.2)
p <- p + geom_errorbar(aes(xmin = conf.low.lm, xmax = conf.high.lm,color = Effect_name,alpha=sig.anova), width=0.5,linewidth = 0.4)
p <- p + geom_point(aes(color = Effect_name,alpha=sig.anova),shape=16,size = 1.25)

# Mark FDR-significant effects with a small star placed beyond the upper CI
p <- p + geom_point(aes(x=conf.high.lm+0.05),shape=8,size=0.6,color="grey30",data=pd %>% filter(p_adj_fdr.anova<=0.05))
p <- p + scale_color_manual(values = pal_use,labels=pal.labels_use)
p <- p + scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.5))

# Expand limits to ensure star markers beyond CI remain visible
p <- p + scale_x_continuous("effect of breed trait on Breed Ancestry Score",breaks=c(-1,0,1),limits=c(min(c(-1.05,pd$conf.low.lm),na.rm=TRUE),max(c(0.95,pd$conf.high.lm),na.rm=TRUE)+0.05))
p <- p + facet_grid(.~Effect_name)
p <- p + theme_bar_estimate() + theme(legend.position = "none",axis.text.y=element_blank())
p_ci <- p
p 


# ============================================================
# Panel C: Cox model log(HR) bars aligned to row order (missing called out)
# ============================================================
pd <- pd %>% select(phenotype, plot_label, logHR, exp_coef.cox, p_adj_fdr.cox) %>% distinct() %>% 
  left_join(order, by = "phenotype") %>%
  left_join(human, by = "phenotype")

# Capture phenotypes without Cox result so they can be annotated as missing
missing <- pd %>% filter(is.na(exp_coef.cox)) %>% mutate(logHR = 0)

pd <- pd %>%
  filter(!is.na(exp_coef.cox)) %>%
  mutate(
    direction = if_else(exp_coef.cox < 1, "protect", "risk")
  )

# Precompute positions for text/markers based on bar direction and significance
annot <- pd %>%
  mutate(
    label_y = case_when(
      direction == "protect" ~ logHR - 0.03,
      direction == "risk"    ~ logHR + 0.1,
      TRUE ~ NA_real_
    ),
    hjust_lab = case_when(
      direction == "protect" ~ 1,
      direction == "risk"    ~ 0,
      TRUE ~ 0.5
    ),
    point_y = case_when(
      p_adj_fdr.cox <= 0.05 & direction == "protect" ~ logHR - 0.02,
      p_adj_fdr.cox <= 0.05 & direction == "risk"    ~ logHR + 0.02,
      TRUE ~ NA_real_
    )
  )

# Avoid exactly-zero bars (helps visibility when rounded to axis)
pd <- pd %>% mutate(logHR = case_when(
  abs(logHR) < 0.005 & logHR < 0  ~ -0.005,
  abs(logHR) < 0.005 & logHR > 0  ~  0.005,
  TRUE ~ logHR))
  
p <- ggplot(pd, aes(y = logHR, x = rank))
p <- p + geom_vline(xintercept = hlines, color = "grey80", linewidth = 0.2)
p <- p + geom_hline(yintercept = 0,color = "grey60", linewidth = 0.2)
p <- p + geom_bar(color="grey20",stat = "identity", width = 0.5)

# Explicit missing annotation for rows with no Cox estimate
p <- p + geom_text(
  data = missing,
  label = "missing", y = 0.02,
  color = "black", hjust = 0, vjust = 0.5, size = 1.25
)

# Significance stars for FDR-significant Cox effects
p <- p + geom_point(
  data = annot %>% filter(!is.na(point_y)),
  aes(y = point_y),
  shape = 8, color = "black", size = 0.75
)

# Use log(HR) axis but label with HR-fold-change semantics
p <- p + scale_y_continuous(
  "mortality hazard ratio",
  breaks = c(-1, -0.5, -0.24, 0, 0.25, 0.5, 1),
  labels = c("0.5×", "0.7×", "0.85×", "1×", "1.2×", "1.4×", "2×"),
  expand = expansion(mult = c(0.05, 0.05))
)

p <- p + scale_x_discrete("", expand = expansion(mult = c(0.0, 0.02)))
p <- p + coord_flip()

p <- p + theme_bar_estimate() + theme(
   axis.text = element_text(size = 4), #.y = element_blank(),
  legend.position = "top",
  legend.justification = c("right", "top"),
  legend.box = "vertical",
  legend.text = element_text(size = 5),
  legend.title = element_text(size = 5, face = "bold"),
  legend.key.size = unit(0.15, "cm"),
  legend.key.spacing.y = unit(0.05, "cm"),
  legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.15),
  legend.margin = margin(1, 1, 1, 1)
)
# Wrap in ggdraw for consistent sizing with other panels
p <- ggdraw() + draw_plot(p, x = 0, y = 0, width = 1, height = 0.943)

p
p_bar_cox <- p



# ============================================================
# Panel D: Cox volcano-like plot (HR vs -log10 FDR), labeled subset
# ============================================================
pd <- pdAll %>% filter(!is.na(exp_coef.cox)) %>% filter(Effect=="lifespan") %>% 
  select(phenotype,plot_label,paper_phenotype_category,p.anova,exp_coef.cox,p_adj_fdr.cox,in_ben,p_comb_fdr) %>% 
    distinct() %>% mutate(print_label=plot_label)
pd <- break_middle_label(pd,print_label)

# "Top group" used for bold highlighting and labeling
pd <- pd %>% mutate(top_group=if_else(p_comb_fdr<=0.05,TRUE,FALSE),
                       top_group_name=if_else(top_group,"pFDR (combined)<0.05","Not top 20%"))
pd$paper_phenotype_category = factor(pd$paper_phenotype_category,levels=c("Plasma metabolites","Clinical analytes")) 
legend_order <- c("Cox and LMER", "Cox only", "LMER only", "other")

# Fill missing Harrison significance flag with FALSE before grouping
pd <- pd %>% replace_na(list(in_ben=FALSE))

# Group controls point aesthetics and legend categories
pd <- pd %>% mutate(group=if_else(in_ben,"Signif. in Harrison 2025",if_else(paper_phenotype_category=="Plasma metabolites","Not signif.","Not included")))
pd$group <- factor(pd$group,levels=c("Signif. in Harrison 2025","Not signif.","Not included","pFDR (combined)<0.05","Not top 20%"))
pd$top_group_name <- factor(pd$top_group_name,levels=c("Signif. in Harrison 2025","Not signif.","Not included","pFDR (combined)<0.05","Not top 20%"))

p <- ggplot(pd, aes(x = exp_coef.cox , y = -log10(p_adj_fdr.cox)))
p <- p + geom_vline(xintercept =1, color = "gray40", linewidth = 0.25, linetype = 1)
p <- p + geom_hline(yintercept = -log10(0.05),color="grey40",linetype=2,linewidth = 0.25)
p <- p + geom_point(aes(color=group,shape=group,size=group), alpha = 0.75)

p <- p + scale_color_manual(name="Description",values=c("Not included"="#f26c2a","Signif. in Harrison 2025"="#66c2a5","Not signif."="#66c2a5","pFDR (combined)<0.05"="grey20","Not top 20%"=NA))
p <- p + scale_size_manual(name="Description",values=c("Not included"=1,"Signif. in Harrison 2025"=1,"Not signif."=0.5,"pFDR (combined)<0.05"=2.25))
p <- p + scale_shape_manual(name="Description",values=c("Not signif."=24,"Signif. in Harrison 2025"=17,"Not included"=16,"pFDR (combined)<0.05"=21))
p <- p + facet_wrap(~paper_phenotype_category,ncol=1)
p <- p + scale_y_continuous(expression(-log[10]*p[FDR]))
p <- p + scale_x_continuous("mortality hazard ratio")
p <- p + theme_point()
p <- p + guides(
  color = guide_legend(title = NULL,title.position = "top", direction = "horizontal",nrow = 2),
  shape = guide_legend(title = NULL,title.position = "top", direction = "horizontal",nrow = 2),
  size  = guide_legend(title = NULL,title.position = "top", direction = "horizontal",nrow = 2)
  )
p_cox <- p

# Overlay highlighted "top group" points with stronger alpha
p_cox <- p_cox + geom_point(aes(color=top_group_name,shape=top_group_name,size=top_group_name),alpha=0.5,data=pd %>% filter(top_group),alpha=1)

# Label the highlighted points
p_cox <- p_cox + geom_text_repel(aes(label = print_label),
                         data = pd %>% filter(top_group),
                         segment.color = "black",
                         force = 1.5,
                         force_pull = 0,
                         box.padding = 0.25,
                         segment.size = 0.25,
                         point.padding = 0.15,
                         size=1.5,
                         lineheight = 0.8,
                         min.segment.length = 0,
                         max.overlaps = 20)

# Supplement: label all Cox FDR-significant points and save
p_supp_cox_labeled <- p + geom_text_repel(aes(label = print_label),
                                          data = pd %>% filter(p_adj_fdr.cox<=0.05),
                                          segment.color = "black",
                                          force = 1.5,
                                          force_pull = 0,
                                          box.padding = 0.25,
                                          segment.size = 0.25,
                                          point.padding = 0.15,
                                          size=1.5,
                                          lineheight = 0.8,
                                          min.segment.length = 0,
                                          max.overlaps = 20)
ggsave(plot=p_supp_cox_labeled,filename=outpdf_labelled,width=5,height=5)

# ============================================================
# Panel E: Cox vs LMER combined Z-scores (colored by combined significance)
# ============================================================

# Aesthetic mappings for Cox-vs-LMER scatter
cox_lmer_shapes <- c("Plasma metabolites"=17,"Clinical analytes"=16,"pFDR (combined)<0.05"=21)
cox_lmer_sizes <- c( "other" = 0.75, "protective" = 1, "risk" = 1,"pFDR (combined)<0.05"=2)
cox_lmer_colors <- c( "other" = "gray40", "protective" = "#1F9BCF", "risk" = "#6A3D9A","pFDR (combined)<0.05"="grey20")
cox_lmer_labels <- c( "other" = "neither", "protective"= "protective", "risk"= "risk","pFDR (combined)<0.05"="pFDR (combined)<0.05")

# Restrict to lifespan effect with available Cox estimate
pd <- pdAll %>% filter(Effect=="lifespan") %>% filter(!is.na(exp_coef.cox)) %>% mutate(print_label=plot_label)
pd <- break_middle_label(pd,print_label)

# Combined-significance grouping uses concordant direction across LMER and Cox Z-scores
pd <- pd %>%
  mutate(
    group = case_when(
      p_comb_fdr <= 0.05 & z_lm < 0 & -z_cox < 0 ~ "risk",
      p_comb_fdr <= 0.05 & z_lm > 0 & -z_cox > 0 ~ "protective",
      TRUE ~ "other"
    )
  )

# Attach phenotype category for shape stratification
pd <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category) %>% right_join(pd)

# Category-adjusted correlation (Spearman) shown on plot
r_cox <- resid(lm(z_cox ~ paper_phenotype_category, data = pd))
r_lm  <- resid(lm(z_lm  ~ paper_phenotype_category, data = pd))
cor_result <- cor.test(r_cox, r_lm, method = "spearman")
r_val <- round(cor_result$estimate, 2)
p_val <- signif(cor_result$p.value, 2)
z_cor_spearman_str <- paste0("Spearman rho = ", r_val, ", p = ", p_val)

# Percentile ranks (used downstream for potential selection/diagnostics)
pd <- pd %>% arrange(exp_coef.cox) %>% 
  mutate(rank_indiv = rank(exp_coef.cox, ties.method = "first") / n())
pd <- pd %>% arrange(beta_std.lm) %>% 
mutate(rank_breed = rank(beta_std.lm, ties.method = "first") / n())

# Choose which points to label (union of row1 phenotypes + any non-"other" + any significant)
top_pd <- pd %>% filter(phenotype %in% order_from_row1$phenotype|group!="other"|p_comb_fdr<=0.05|p.anova<=0.05)
top_pd <- top_pd %>% mutate(orig_plot_label=plot_label)
top_pd <- break_middle_label(top_pd,plot_label)

# Set factor order for legends and facets
pd$group <- factor(pd$group,levels=c("risk","protective","other"))
pd$paper_phenotype_category <- factor(pd$paper_phenotype_category,levels=c("Plasma metabolites","Clinical analytes"))

pd %>% filter(p_comb_fdr<=0.05) %>% select(plot_label,exp_coef.cox,beta_std.lm,support_score,p_comb_fdr)

p <- ggplot(pd, aes(x = z_cox, y = z_lm))
p <- p + geom_smooth(method = "lm", se = TRUE, color = "#08306b",fill="#c6dbef",linewidth=0.25)
p <- p + geom_point(aes(color = group, shape = paper_phenotype_category, size= group), alpha = 0.75)

# Emphasize combined-FDR hits with an outlined point
p <- p + geom_point(aes(color = group),shape=21,size=3,data=pd %>% filter(p_comb_fdr<=0.05),alpha=1)

# Correlation annotation at lower-left corner of the panel
p <- p + annotate("text",y = min(pd$z_lm), x = min(pd$z_cox),color="black",hjust = 0, vjust = 0,label = z_cor_spearman_str,size = 2)

# Repelled labels for selected points
p <- p + geom_text_repel(aes(label = print_label,color=group),
                         data = top_pd,
                         segment.color = "black",
                         force = 1.5,
                         force_pull = 0,
                         box.padding = 0.25,
                         segment.size = 0.25,
                         point.padding = 0.15,
                         size=1.5,
                         lineheight = 0.8,
                         min.segment.length = 0,
                         max.overlaps = 20)
p <- p + scale_color_manual(values=cox_lmer_colors,labels = cox_lmer_labels)
p <- p + scale_shape_manual(values=cox_lmer_shapes,labels = cox_lmer_labels)
p <- p + scale_size_manual(values=cox_lmer_sizes,labels = cox_lmer_labels)
p <- p + scale_y_continuous("Z score (effect of breed lifespan on Breed Ancestry Score)")
p <- p + scale_x_continuous("Z score (mortality hazard ratio)")
p <- p + theme_cowplot(6)
p <- p + theme(
  axis.line = element_line(linewidth = 0.25),
  plot.margin = margin(t = 5, r = 5, l = 5, b = 5),
  axis.ticks = element_line(linewidth = 0.5),
  panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
  strip.text = element_text(size = 7, hjust = 0, vjust = 0, face = "bold", margin = margin(t = 2, b = 2), color = "black"),
  strip.background = element_rect(fill = "white", color = NA),
  legend.position = c(0.98, 0.98),
  legend.justification = c("right", "top"),
  legend.text = element_text(size = 4),
  legend.title = element_text(size = 4,face="bold"),
  legend.key.size = unit(0.1, "cm"),
  legend.key.spacing.y = unit(0.02, "cm"),
  legend.margin = margin(1, 1, 1, 1),
  legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.15)
)
p <- p + guides(
  color = guide_legend(title = NULL),
  shape = guide_legend(title = NULL),
  size  = guide_legend(title = NULL)
)
p_cox_v_lmer <- p
p


# ============================================================
# Combine panels and write final figure
# ============================================================
row1 <- plot_grid(p_bar,p_ci,p_bar_cox,nrow=1,label_size=8,labels=LETTERS,rel_widths = c(0.65,1,0.5))
row2  <- plot_grid(p_cox,p_cox_v_lmer,nrow=1,label_size=8,labels=c("D","E"),rel_widths = c(0.5,1))
grid <- plot_grid(row1,row2,ncol=1,rel_heights = c(1,1.6))
ggsave(plot=grid,filename=outpdf,width=6,height=5.5)
