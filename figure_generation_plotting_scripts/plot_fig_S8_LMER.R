# ==============================================================================
# Setup: load only the packages used in this script and shared utilities
# ==============================================================================
library(tidyverse)  # dplyr/tibble/ggplot2 pipelines and plotting
library(cowplot)    # theme_cowplot() and plot_grid() for multi-panel figures
source("util.cleaning.R")  # project-specific helper functions (loaded for consistency)

# Output path for the final figure PDF
outpdf <- "../../../figures_pdf/fig_S8_LMER.pdf"

# ==============================================================================
# Load all data objects from the RDS (only if not already present in the session)
# ==============================================================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects like raw_S2_ALL_PHENOS, raw_S9_ANOVA_ON_LMER, etc.
}

# ==============================================================================
# Plot theme helpers
# ==============================================================================

# Theme for the stacked bar chart (variance explained)
theme_bar <- function(){
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
      legend.position = "bottom",  
      
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5,face="bold"),
      legend.margin = margin(2,2,2,2), 
      legend.box.background = element_rect(fill="white",color = "grey50", linewidth = 0.15) ,
      
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.grid.major.x=element_line(color="grey60",linewidth=0.2),
      panel.grid.major.y=element_line(linewidth = 0.2,color="grey80")   
      
    )
}

# Theme for the standardized effect-size + CI panel
theme_bar_estimate <- function(){
  theme_cowplot(6) %+replace%
    theme(
      plot.title = element_text(size = 5,face="bold",hjust=0),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(linewidth=0.25),
      axis.line.y=element_blank(),
      axis.line.x = element_line(linewidth = 0.25),    
      axis.text.x = element_text(size = 5),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 6),
      panel.grid.major.x=element_line(color="grey60",linewidth=0.2),
      panel.grid.major.y=element_line(linewidth = 0.2,color="grey80"),    
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      legend.position = "bottom",  
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5,face="bold"),
      legend.margin = margin(2,2,2,2), 
      legend.box.background = element_rect(fill="white",color = "grey50", linewidth = 0.15),
      legend.box.spacing = unit(7, "pt")
      
    )
}

# ==============================================================================
# Prepare phenotype metadata and effect label mapping
# ==============================================================================
pheno <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label) 

# Human-readable names and plotting order for the two predictor effects
effect_names <- tibble(
  Effect      = c("breed.weight.kg","lifespan"),
  Effect_name = c("breed weight","breed lifespan"),
  order       = c(2,1)
)

# ==============================================================================
# Filter and assemble ANOVA/LMER summary for the target model
# ==============================================================================
pd <- raw_S9_ANOVA_ON_LMER %>% filter(model_used=="REML.t.val ~ breed.weight.kg + lifespan") %>% 
  select(phenotype,Effect,ges=ges.anova,p.anova,p_adj_fdr.anova,beta_std.lm,std.error.lm,conf.low.lm, conf.high.lm) %>% 
  drop_na() 
pd <- pd %>% mutate(sig.anova=if_else(p.anova<=0.05,TRUE,FALSE))

# Attach phenotype category/labels and effect display names
pd <- pheno %>% right_join(pd) 
pd <- effect_names %>% right_join(pd) 

# Keep only phenotypes that are significant for at least one effect in this model
ok_pheno <- pd %>% filter(p.anova<=0.05) %>% select(phenotype) %>% distinct()
pd <- pd %>% inner_join(ok_pheno)

# Rank phenotypes by total variance explained across effects (for ordering in plots)
pd_ranks <- pd %>% group_by(plot_label) %>% 
  summarize(ges_total = sum(ges, na.rm = TRUE), .groups = "drop") %>%
  arrange(-ges_total) %>%
  mutate(rank = row_number())
  
# Y-axis limit for variance-explained bars
maxy <- max(pd_ranks$ges_total, na.rm = TRUE)

# Add ranking info back to the long table
pd <- pd %>% left_join(pd_ranks, by = "plot_label")
  
# Set factor levels to control plot ordering
pd$plot_label  <- factor(pd$plot_label, levels = pd_ranks %>% arrange(-rank) %>% pull(plot_label))
pd$Effect_name <- factor(pd$Effect_name, levels = effect_names %>% arrange(-order) %>% pull(Effect_name))

# ==============================================================================
# Panel A: stacked bars of generalized eta-squared (variance explained)
# ==============================================================================
p <- ggplot(pd  %>% filter(p.anova<=0.05) , aes(y = ges, x = plot_label))
p <- p + geom_bar(aes(fill = Effect_name), stat = "identity", width = 0.5) + coord_flip()
p <- p + facet_grid(paper_phenotype_category ~ ., space = "free_y", scales = "free_y")
p <- p + scale_y_continuous("variance explained",
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
      expand = expansion(mult = c(0, 0)),
      limits = c(0, maxy + 0.01)
    )
p <- p + scale_fill_manual("Effect",values = c(
      "breed weight"="#3070ad",
      "breed lifespan"="#971d2b"
    ))
p <- p + theme_bar() + theme(plot.margin = margin(t = 20, r = 5, l = 0, b = 0))

p_effects <- p

# ==============================================================================
# Panel B: standardized effect estimates with confidence intervals
# ==============================================================================
# Used only to determine line positions/plot dimensions (kept for parity with original script)
maxy <- length(unique(pd$phenotype))
hlines<- seq(from = 0.5, to = maxy, by = 2)

# Debug/inspection print of labels present (kept as executable line)
unique(pd$plot_label)

p <- ggplot(pd, aes(y = plot_label, x = beta_std.lm))
p <- p + geom_vline(xintercept=0)

# CI bars and points; alpha indicates nominal significance
p <- p + geom_errorbar(aes(xmin = conf.low.lm, xmax = conf.high.lm,color = Effect_name,alpha=sig.anova), width=0.5,linewidth = 0.4)
p <- p + geom_point(aes(color = Effect_name,alpha=sig.anova),shape=16,size = 1.25)

# Mark effects passing FDR threshold with an "x" at the upper CI bound (nudged right)
p <- p + geom_point(aes(x=conf.high.lm+0.05),shape=8,size=0.6,color="grey30",data=pd %>% filter(p_adj_fdr.anova<=0.05))

p <- p + scale_color_manual("Effect",values = c(
  "breed weight"="#3070ad",
  "breed lifespan"="#971d2b"
))

# Map TRUE/FALSE significance to alpha levels for visual emphasis
p <- p + scale_alpha_manual("p < 0.05",values=c("TRUE"=1,"FALSE"=0.5))

# X-axis range padded to ensure CI markers fit within bounds
p <- p + scale_x_continuous("Standardized effect of breed trait on Breed Ancestry Score",breaks=c(-1,0,1),limits=c(min(c(-1.05,pd$conf.low.lm),na.rm=TRUE),max(c(0.95,pd$conf.high.lm),na.rm=TRUE)+0.05))
p <- p + facet_grid(paper_phenotype_category ~ Effect_name, space = "free", scales = "free")

p <- p + theme_bar_estimate() + theme(plot.margin = margin(t = 12, r = 5, l = 0, b = 0))

p_ci <- p

# ==============================================================================
# Combine panels and save
# ==============================================================================
grid <- plot_grid(p_effects,p_ci,ncol=2,rel_widths = c(1,1),labels=LETTERS,label_size=8)
grid
ggsave(plot=grid,filename=outpdf,width=6,height=6)
