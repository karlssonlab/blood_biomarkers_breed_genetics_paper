# =============================================================================
# Dependencies
# =============================================================================
library(tidyverse)  # Data manipulation + ggplot2 (also supplies stringr/purrr/tidyr)
library(cowplot)    # Publication-style themes and plot_grid composition
library(ggrepel)    # Non-overlapping text labels for ggplot
library(rstatix)    # kruskal_test and related helpers

# =============================================================================
# Load data objects and helper functions
# =============================================================================
# Load the supplemental data once and unpack it into the global environment
# (creates objects like raw_S2_ALL_PHENOS, raw_S8_LMER, raw_S12_BREED_PHENOS, etc.)
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # optional: creates raw_pheno, raw_anc, etc.
}

# Project utilities (data cleaning + breed image plotting)
source("util.cleaning.R")
source("util.make_breed_image.R")

# Output file for saving the composed figure objects
out_RData <- "fig_LMER.for_fig_4_FUR.RData"

# =============================================================================
# Prepare LM/REML results table and breed phenotype metadata
# =============================================================================
# Join phenotype labels onto the LM/REML results and drop incomplete rows
lmer <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label) %>% 
  distinct() %>% inner_join(raw_S8_LMER) %>% drop_na()

# Standardize breed naming across inputs so joins/labels are consistent
lmer <- standardize_breed_names(lmer,breed)
breed_pheno <- standardize_breed_names(raw_S12_BREED_PHENOS,breed)

# Convert yes/no phenotype encodings to logical where applicable
breed_pheno <- convert_yesno_to_logical(breed_pheno)

# Create compact plotting labels and derived coat-length indicators
breed_pheno <- breed_pheno %>% mutate(lmer_plot_breed = map_chr(breed, replace_middle_break))
breed_pheno <- breed_pheno %>% mutate(lmer_plot_breed=str_replace(lmer_plot_breed,"american","am."))
breed_pheno <- breed_pheno %>% mutate(length_med=if_else(length=="medium",TRUE,FALSE),
                                      length_long=if_else(length=="long",TRUE,FALSE))

# =============================================================================
# Plot styling helper
# =============================================================================
# Shared theme for small multi-panel figure elements
theme_lmer <- function(){ 
  theme_cowplot(6) %+replace%
    theme(axis.line=element_line(linewidth=0.25),
          plot.margin = margin(t=5,r=5,l=5,b = 5),
          axis.ticks=element_line(linewidth=0.25),
          legend.position = "NONE",
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.25),
          panel.grid.major.y=element_line(color="grey80",linewidth=0.2),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size=6,face = "bold")
          
    )}

# =============================================================================
# Function: scatter of breed ancestry score vs breed weight for one phenotype
# =============================================================================
plot_lmer_vs_weight <- function(phenotype_name,
                                lmer_df = lmer,
                                breed_pheno_df = breed_pheno) {
  
  # Subset LM/REML results for the phenotype and merge in breed weight + plot labels
  pd <- lmer_df %>%
    filter(phenotype == phenotype_name) %>%
    select(plot_label, breed, REML.t.val,
           ML.anova.p, ML.anova.p.adj_benjhoch_FDR) %>%
    inner_join(
      breed_pheno_df %>% select(breed, lmer_plot_breed, breed.weight.kg) %>% distinct(),
      by = "breed"
    ) %>%
    mutate(sig = ML.anova.p.adj_benjhoch_FDR <= 0.05)
  
  # Spearman correlation between ancestry score and breed weight
  corr <- cor.test(pd$REML.t.val, pd$breed.weight.kg, method = "spearman")
  corr.p <- if (corr$p.value < 0.001) format(corr$p.value, scientific = TRUE, digits = 2) else signif(corr$p.value, 2)
  # Format correlation p-value as plotmath (coefficient * 10^exponent)
  p_exp   <- floor(log10(corr$p.value))
  p_coeff <- corr$p.value / 10^p_exp
  
  corr_string <- paste0(
    "Spearman~r==", round(corr$estimate, 2),
    "*','~~p==", round(p_coeff, 1), "%*%10^", p_exp
  )
  print(length(unique(pd$breed)))
  print(corr)
  
  # Label the most extreme significant breeds (top 3 and bottom 3 by REML.t.val)
  top <- bind_rows(
    pd %>% filter(sig) %>% arrange(desc(REML.t.val)) %>% slice_head(n = 3),
    pd %>% filter(sig) %>% arrange(REML.t.val) %>% slice_head(n = 3)
  )
  
  # Axis labels & limits
  xlabel_raw <- unique(pd$plot_label)[1]
  xlabel <- paste0(toupper(substr(xlabel_raw,1,1)), substr(xlabel_raw,2,nchar(xlabel_raw)))
  minX <- floor(min(pd$REML.t.val,na.rm=T))
  maxX <- ceiling(max(3, max(pd$REML.t.val) * 1.05,na.rm=T))
  maxY <- max(pd$breed.weight.kg,na.rm=T)+5
  print(minX)
  print(maxY)
  
  # Plot: points (significant highlighted), linear trend, and text labels
  p <- ggplot(pd, aes(x = REML.t.val, y = breed.weight.kg))
  p <- p + geom_point(aes(color = sig, shape = sig), alpha = 0.7)
  p <- p + geom_smooth(method = 'lm', fill = NA, color = "#a50f15", alpha = 0.5, linewidth = 0.5)
  p <- p + geom_vline(xintercept = 0, color = "gray40", linewidth = 0.25, linetype = 2)
  p <- p + geom_text_repel(data = top, aes(label = lmer_plot_breed),
                                    force = 5, color = "#a50f15",
                                    segment.color = "#a50f15", box.padding = 0.5,
                                    segment.size = 0.25, point.padding = 0.3,
                                    size = 1.5, lineheight = 0.8,
                                    min.segment.length = 0, max.overlaps = 20)
  p <- p + annotate("text", label = corr_string, parse=TRUE,x = minX, y = maxY+2,
                    color = "#a50f15", size = 1.75, hjust = 0, vjust = 1)
  p <- p + scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#a50f15"))
  p <- p + scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 16))
  p <- p + scale_y_continuous("breed average weight (kg)",limits=c(-20,maxY+4))
  p <- p + scale_x_continuous(paste(xlabel, "Breed Ancestry Score"),
                              limits = c(minX, maxX),
                              expand = expansion(mult = c(0.01, 0.1)))
  p <- p + theme_lmer()
  
  return(p)
}

# =============================================================================
# Panel: creatinine (ln transformed) ancestry score vs breed weight
# =============================================================================
p_size_lmer <- plot_lmer_vs_weight("krt_cp_creatinine_value_ln_transformed")

# Identify smallest/largest breeds by ancestry score for the breed image strips
pd <- lmer %>% filter(phenotype == "krt_cp_creatinine_value_ln_transformed") 
  
if (!exists("small3")){
small3 <- pd %>% arrange(REML.t.val) %>% head(3) %>% pull(breed)
}
if (!exists("big3")){
big3 <- pd %>% arrange(REML.t.val) %>% tail(3) %>% pull(breed)
}
breed_size_plot <- make_breed_figure_custom_lists("breed.weight.kg", low_breeds=small3,high_breeds=big3,max_px=400)

# =============================================================================
# Panel: Cystathionine ancestry score vs coat type (hair length/furnishings)
# =============================================================================
# Mapping from phenotype indicator combinations to plot group labels and y positions
phenostrings <- tibble(
  length_long = c(TRUE, FALSE, TRUE),
  furnishings = c(FALSE,TRUE,  TRUE),
  phenostr    = c("long_no_furnishings","furnishings","furnishings"),
  ypos = c(2,3,3))

# Combine LM/REML results with breed phenotype fields needed to classify coat type
pd <- lmer %>% filter(phenotype=="Cystathionine") %>% select(plot_label,breed,REML.t.val,ML.anova.p,ML.anova.p.adj_benjhoch_FDR)
pd <- breed_pheno %>% select(breed,lmer_plot_breed,length_med,length_long,furnishings) %>% inner_join(pd)
pd <- pd %>% mutate(length_short=if_else(!length_med&!length_long, TRUE,FALSE))

# Assign each breed to a coat-type string; anything unmatched becomes "other breeds"
pd <- pd %>% left_join(phenostrings) %>% replace_na(list(phenostr="other\nbreeds",ypos=1))

# Kruskal-Wallis tests comparing REML.t.val between breeds with vs without each coat feature
# (builds a wide TRUE/FALSE table per breed, then tests each feature)
corr <- pd %>%
  select(breed, REML.t.val, phenostr) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = phenostr) %>%
  select(-`other\nbreeds`) %>%
  replace_na(list(
    long_no_furnishings  = FALSE,
    furnishings = FALSE
  )) %>%
  pivot_longer(c(long_no_furnishings, furnishings)) %>%
  group_by(name) %>%
  kruskal_test(REML.t.val ~ value) %>%
  ungroup() %>%
  mutate(
    # Plain-text p-value label (kept for debugging/inspection)
    p_label_plain = if_else(
      p < 0.001,
      format(p, scientific = TRUE, digits = 2),
      format(p, scientific = FALSE, nsmall = 3, drop0trailing = TRUE)
    ),
    
    # Plotmath label for ggplot (scientific format only when very small p)
    corr_string = if_else(
      p < 0.001,
      {
        exp10  <- floor(log10(p))
        coeff <- p / 10^exp10
        paste0("p==", signif(coeff, 2), "%*%10^", exp10)
      },
      paste0("p==", format(p, scientific = FALSE, nsmall = 3, drop0trailing = TRUE))
    )
  ) %>%
  left_join(phenostrings %>% select(name = phenostr, ypos))
print(corr)
print(length(unique(pd$breed)))

# Mark significant breeds by adjusted p-value and label the significant extremes
pd <- pd %>% mutate(sig=if_else(ML.anova.p.adj_benjhoch_FDR<=0.05,TRUE,FALSE))
top <- pd %>% filter(sig) %>% arrange(REML.t.val)

# Axis labels and limits
xlabel <- pd %>% select(plot_label) %>% distinct() %>% pull(plot_label)
minX <- floor(min(pd$REML.t.val))
maxX <- ceiling(max(3,max(pd$REML.t.val)*1.05))

# Pretty y-axis labels for the coat-type groupings
phenolabels <- pd %>% select(phenostr,ypos) %>% distinct() %>% 
  mutate(ylabel=str_replace_all(phenostr,"\\_"," ")) %>% 
  mutate(ylabel=str_replace_all(ylabel," and long","\nand long fur")) %>% 
  mutate(ylabel=str_replace_all(ylabel," no furnishings"," fur (no\nfurnishings)")) %>% 
  mutate(ylabel=if_else(ylabel=="furnishings","furnishings\n(long, coarse\nfacial hair)",ylabel))

# Median tick marks per group for visual reference
medians <- pd %>%
  group_by(ypos) %>%
  summarize(med_REML = median(REML.t.val, na.rm = TRUE))

# Add y-span columns used by ggrepel for label placement
top <- top %>%
  mutate(
    ymin = ypos - 0.4,
    ymax = ypos + 0.4
  )

# Plot: jittered points per group, medians, and p-values from Kruskal tests
p <- ggplot(pd,aes(x=REML.t.val,y=ypos))
p <- p + geom_segment(data = medians,aes(x = med_REML,xend = med_REML,y = ypos - 0.2,yend = ypos + 0.2),
                      inherit.aes = FALSE,color = "gray20",linewidth = 0.4)

p <- p + geom_jitter(color="gray40",shape=21,data=pd %>% filter(!sig),width = 0,alpha = 0.7, height = 0.1)
p <- p + geom_point(color="#a50f15",shape=16,alpha=0.7,data=pd %>% filter(sig))
p <- p + geom_vline(xintercept = 0, color = "gray40", linewidth = 0.25, linetype = 2) 
p <- p + geom_text_repel(data = top,aes(label = lmer_plot_breed), force = 3, color = "#a50f15",
                         segment.color = "#a50f15", box.padding = 0.5, segment.size = 0.2, 
                         point.padding = 0.3, size = 1.5, lineheight = 0.8, min.segment.length = 0,
                         max.overlaps = 20)

p <- p + geom_text(aes(label=corr_string,y=ypos+0.5),x=minX,color = "#a50f15",parse=TRUE,size=1.75,hjust=0,vjust=1,data=corr)
p <- p + scale_y_continuous("",labels=phenolabels$ylabel,breaks=phenolabels$ypos,limits=c(0,4))
p <- p + scale_x_continuous(paste(xlabel,"Breed Ancestry Score"),limits=c(minX,maxX),expand = expansion(mult = c(0.01,0.1)))
p <- p + theme_lmer() + theme(axis.text.y = element_text(size=5,face = "bold"))

# Store cystathionine/fur panel object
p_fur_lmer <- p
p

# =============================================================================
# Breed image strips for fur panel (choose extremes by ancestry score)
# =============================================================================
other3 <- pd %>% arrange(REML.t.val) %>% head(3) %>% pull(breed)
oddfur3 <- pd %>% filter(breed != "poodle (toy)") %>% arrange(REML.t.val) %>% tail(3) %>% pull(breed)
breed_fur_plot <- make_breed_figure_custom_lists("odd_fur", low_breeds=other3,high_breeds=oddfur3,max_px=400)

# =============================================================================
# Compose final multi-part figures and save
# =============================================================================
# Tighten margins on the breed-image strips before combining
breed_size_plot$low <- breed_size_plot$low + theme(plot.margin = margin(l=0,r=5,t=-5,b=10))
breed_size_plot$high <- breed_size_plot$high + theme(plot.margin = margin(l=5,r=0,t=-5,b=10))
breed_size_plot_all <- plot_grid(breed_size_plot$low ,breed_size_plot$high,nrow=1,align = "v",axis = "b")

breed_fur_plot$low <- breed_fur_plot$low + theme(plot.margin = margin(l=0,r=5,t=-5,b=10))
breed_fur_plot$high <- breed_fur_plot$high + theme(plot.margin = margin(l=5,r=0,t=-5,b=10))
breed_fur_plot_all <- plot_grid(breed_fur_plot$low ,breed_fur_plot$high,nrow=1, align = "v",axis = "b")

# Add a spacer column to shift the fur breed strip rightward
breed_fur_plot_all <- plot_grid(NULL,breed_fur_plot_all,nrow=1,rel_widths = c(0.1,1))

# Stack the scatter panel over the corresponding breed-image strip
grid_lmer_fur <- plot_grid(p_fur_lmer,breed_fur_plot_all,ncol=1,rel_heights = c(1.2,0.4))
grid_lmer_size <- plot_grid(p_size_lmer,breed_size_plot_all,ncol=1,rel_heights = c(1.2,0.4))

# Save the composed plot objects for downstream figure assembly
save(grid_lmer_size, grid_lmer_fur, file = out_RData)
