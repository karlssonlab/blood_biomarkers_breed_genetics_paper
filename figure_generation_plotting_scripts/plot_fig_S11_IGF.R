# Libraries used for data wrangling and plotting
library(tidyverse)      # Data manipulation and ggplot2
library(cowplot)        # Theme and plot styling
library(ggpubr)         # stat_compare_means for significance testing

# Helper functions (project-specific cleaning/utilities)
source("util.cleaning.R")

# Output path for the final figure
outpdf <- "../../../figures_pdf/fig_S11_IGF.pdf"


## =====================================
## Load all data objects from the RDS (once per session)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_* objects in the global environment
}

## =====================================
## Build analysis table: select model + annotate GH-linked metabolites
## =====================================
model_to_use <- "REML.t.val ~ breed.weight.kg + lifespan"

# Start from S11 combined results for the chosen model and keep complete cases
pdAll <- raw_S11_COMBINED %>% filter(model_used==model_to_use) %>% 
  select(-model_used,-in_cox_analysis) %>% 
  drop_na()

# Attach metabolite plot labels (Plasma metabolites only)
pdAll <- raw_S2_ALL_PHENOS %>% filter(paper_phenotype_category=="Plasma metabolites") %>% 
  select(phenotype,plot_label) %>% distinct() %>% 
  distinct() %>% inner_join(pdAll) 

# Join GH/IGF1 linkage annotation; treat missing linkage as not GH-related
pdAll <- pdAll %>% left_join(raw_S6_IGF1 %>% select(phenotype,GH=`linked_with_GH-IGF1_in_rodent_model`)) %>%
  replace_na(list(GH=FALSE)) %>% ungroup() 

# Convenience group labels used for plotting
pdAll <- pdAll %>% mutate(set=if_else(GH,"GH","not GH")) 

pdAll <- pdAll %>%
  mutate(GH_group = if_else(GH, "GH", "non_GH"))

## =====================================
## Summary table: mean absolute z-scores by effect and score type, split by GH
## =====================================
means <- pdAll %>%
  # reshape z columns long
  pivot_longer(
    cols = c(z_lm, z_cox, z_comb),
    names_to = "z_score",
    values_to = "z"
  ) %>%
  mutate(abs_z = abs(z)) %>%
  filter(!is.na(abs_z), !is.na(GH), !is.na(Effect)) %>%
  group_by(Effect, z_score) %>%
  # means + n for each GH group
  mutate(
    n_groups = n_distinct(GH),
    mean_abs_z_GHTRUE  = mean(abs_z[GH], na.rm = TRUE),
    mean_abs_z_GHFALSE = mean(abs_z[!GH], na.rm = TRUE),
    n_GHTRUE  = sum(GH, na.rm = TRUE),
    n_GHFALSE = sum(!GH, na.rm = TRUE)
  ) %>%
  # one-row-per group with the summary columns
  summarise(
    mean_abs_z_GHTRUE  = first(mean_abs_z_GHTRUE),
    mean_abs_z_GHFALSE = first(mean_abs_z_GHFALSE),
    n_GHTRUE           = first(n_GHTRUE),
    n_GHFALSE          = first(n_GHFALSE),
    n_groups           = first(n_groups),
    .groups = "drop"
  ) 

## =====================================
## Create plotting data for violin/boxplots across three score definitions
## =====================================
# Cox model absolute z-score (hazard ratio association)
pdBox <- pdAll %>% mutate(abs_z_cox = abs(z_cox)) %>% select(phenotype,plot_label,GH,value=abs_z_cox) %>% mutate(set="Cox mortality\nhazard ratio") %>% distinct()

# Linear model absolute z-score for weight and lifespan effects (with readable facet labels)
pdBox <- pdAll %>% select(phenotype,plot_label,GH,set=Effect,value=z_lm) %>% mutate(value=abs(value)) %>% mutate(set=if_else(set=="lifespan","Effect of breed lifespan\non Breed Ancestry Score",
                                                                                                                             "Effect of breed weight\non Breed Ancestry Score")) %>% bind_rows(pdBox)

# Combined score (only defined for lifespan effect here)
pdBox <- pdAll %>% filter(Effect=="lifespan") %>% mutate(abs_z_comb=abs(z_comb)) %>% select(phenotype,plot_label,GH,value=abs_z_comb) %>% mutate(set="Combined score\n(Cox + breed lifespan)") %>% bind_rows(pdBox)

## =====================================
## Plot: GH vs non-GH distribution comparison (Wilcoxon test)
## =====================================
ggplot(pdBox,aes(x=GH,y=value)) + geom_boxplot() + facet_wrap(~set,nrow=1)

# Define pairwise comparison between the two GH groups (logical values shown as strings)
my_comparisons <- list(c("FALSE","TRUE"))

# Violin + boxplot with Wilcoxon comparison per facet
p <- ggplot(pdBox ,aes(x=value,y=GH))
p <- p + geom_violin(aes(color=GH,fill=GH),alpha=0.25)
p <- p + geom_boxplot(fill=NA,outlier.size=0.5,outlier.alpha=0.25,outlier.shape=16,linewidth=0.25,width=0.25) 
p <- p + scale_alpha_manual(values=c(0.9,0.2))
p <- p + stat_compare_means(comparisons = my_comparisons,size=1.75,method="wilcox.test")
p <- p + scale_color_manual(values=c("TRUE"="#8da0cb","FALSE"="#66c2a5"))
p <- p + scale_fill_manual(values=c("TRUE"="#8da0cb","FALSE"="#66c2a5"))
p <- p + scale_y_discrete("",breaks=c("TRUE","FALSE"),labels=c("Growth\nhormone\nrelated","others"),expand = expansion(add = 0.6) )
p <- p + scale_x_continuous("absolute z-score")
p <- p + facet_wrap(~set,ncol=2)
p <- p +
  theme_cowplot(font_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 7),
    
    legend.position = "none",
    
    # lighter axes
    axis.line = element_line(linewidth = 0.3),
    
    # y ticks off, x ticks longer
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.3),
    axis.ticks.length.x = unit(0.1, "cm"),  # ~2× default
    
    strip.text = element_text(
      size = 7,
      hjust = 0.5,
      face = "bold",
      lineheight = 1,
      margin = margin(t = 5, r = 3, b = 5, l = 3)
    ),
    strip.background = element_rect(fill = "white", color = NA),
    
    plot.margin = margin(r = 20)
  )
p_box <- p
p


# Save the assembled plot to PDF
ggsave(plot=p_box,filename=outpdf,width=4,height=3)
