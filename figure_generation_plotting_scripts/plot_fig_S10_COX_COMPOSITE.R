# ---- Setup: packages, helpers, and output path ----
library(tidyverse)      # dplyr/stringr/ggplot2 pipelines and plotting
library(cowplot)        # theme_cowplot()
source("util.cleaning.R")

# ---- I/O: output figure path ----
outpdf <- "../../../figures_pdf/fig_S10_COX_COMPOSITE.pdf"

# ---- Data load: read RDS once and unpack into the global environment ----
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects like raw_S2_ALL_PHENOS, raw_S5_COX_HAZ_W_BREED_ANC, etc.
}

# ---- Assemble plotting data: join phenotype labels onto Cox model results ----
pd <- raw_S2_ALL_PHENOS %>% select(phenotype,plot_label) %>% right_join(raw_S5_COX_HAZ_W_BREED_ANC)

# ---- Create display labels and y-axis ordering for terms/phenotypes ----
names <- pd %>% select(term,phenotype,plot_label) %>% distinct() 
names <- names %>% mutate(label=if_else(str_detect(term,"risk_breeds_composite_ancestry"),"Composite ancestry score",term))
names <- names %>% mutate(phenotype=str_remove(phenotype,"_ln_transformed"))

# If label matches phenotype (i.e., a categorical phenotype), annotate as "<plot_label> levels"
names <- names %>% mutate(label=if_else(phenotype==label,paste0(plot_label," levels"),label))

# Human-readable replacements for common covariates
names <- names %>% mutate(label=str_replace(label,"bestWt","Body weight")) %>% 
  mutate(label=str_replace(label,"recruitmentAge","Age at enrollment")) %>% 
  mutate(label=str_replace(label,"sexMale","Sex (male)")) 

# Put key covariates below the main block (higher order number = lower on the plot after reversing)
orderY <- tibble(label=c("Composite ancestry score",
                        "Sex (male)","Age at enrollment","Body weight"),order=c(2:5))
names <- names %>% left_join(orderY) %>% replace_na(list(order=1)) # default order=1 for all other terms
orderY <- names %>% select(label,order) %>% distinct() %>% arrange(order)

# ---- Apply labels and y-axis factor levels to plotting data ----
pd <- pd %>% left_join(names %>% select(term,label) %>% distinct())
orderY <- orderY %>% inner_join(pd %>% select(label) %>% distinct())
pd$label <- factor(pd$label,levels=orderY %>% arrange(-order) %>% pull(label))

# ---- Add p-value strings and significance flag for styling ----
pd <- pd %>% mutate(corr_str=if_else(p_value<0.001,"p<0.001",paste0("p=",round(p_value,3))))
pd <- pd %>% mutate(sig=if_else(p_value<=0.05,TRUE,FALSE))
max_x = max(pd$conf_high)+0.25

# ---- Plot: hazard ratio points with confidence intervals and annotations ----
p <- ggplot(pd,aes(x=estimate,y=label)) + geom_point(aes(color=sig),shape=16,size=1)
p <- p + geom_errorbar(aes(xmin = conf_low, xmax = conf_high,color=sig), width = 0.25)
p <- p + geom_text(aes(label = round(estimate,2),color=sig), nudge_y = 0.25, vjust = 0, hjust = 0.5, size=1.5)
p <- p + geom_text(aes(label = corr_str, x = conf_high+0.12,color=sig), vjust = 0.5, hjust = 0,size=1.9)
p <- p + geom_vline(xintercept = 1, linetype = 2, color = "black", linewidth = 0.25)

# Facet by phenotype panel label
p <- p + facet_grid(plot_label ~ ., space = "free_y", scales = "free_y")
p <- p + scale_color_manual(values = c("TRUE"  = "black","FALSE" = "grey40"))
p <- p + scale_y_discrete("",expand = expansion(mult = c(0.15,0.15)))
p <- p + scale_x_continuous("hazard ratio", limits=c(0,max_x),breaks=c(0:ceiling(max_x)),expand = expansion(mult = c(0,0.1)))

# ---- Styling ----
p <- p + theme_cowplot(7)
p <- p + theme(
  axis.line = element_line(linewidth = 0.25),
  legend.position="none",
  plot.margin = margin(t = 5, r = 5, l = 5, b = 5),
  axis.ticks.x = element_line(linewidth = 0.5),
  axis.ticks.y = element_blank(),
  strip.text = element_blank(),                     # remove strip labels
  strip.background = element_blank(),               # remove strip boxes
  panel.spacing.y = unit(3, "mm"),                   # increase vertical space between facets
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
)

# ---- Output: render and save PDF ----
p
ggsave(plot=p,filename=outpdf,width=5,height=6)
