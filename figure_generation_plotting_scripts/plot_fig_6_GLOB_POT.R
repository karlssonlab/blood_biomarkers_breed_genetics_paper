# =============================================================================
# Libraries
# =============================================================================
library(tidyverse)      # Data manipulation (dplyr/tidyr/purrr/stringr) and plotting (ggplot2)
library(cowplot)        # Plot composition (plot_grid) and theming (theme_cowplot)
library(ggrepel)        # Non-overlapping text labels (geom_text_repel)
library(survival)       # Survival objects and model fitting (Surv, survfit)
library(broom)          # Tidy survival fit objects into data frames for ggplot

# =============================================================================
# Project utilities
# =============================================================================
source("util.cleaning.R")
source("util.make_breed_image.R",echo = TRUE)

# Output path for the assembled figure PDF
outpdf <- "../../../figures_pdf/fig_6_GLOB_POT.pdf"

# =============================================================================
# Load all data objects (from a single RDS list) into the global environment
# =============================================================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects like raw_S2_ALL_PHENOS, raw_S8_LMER, etc.
}

# =============================================================================
# Theme for LMER scatter panels
# =============================================================================
theme_lmer <- function(){ 
  theme_cowplot(6) %+replace%
    theme(axis.line=element_line(linewidth=0.25),
          plot.margin = margin(t=5,r=5,l=5,b = 5),
          axis.ticks=element_line(linewidth=0.5),
          legend.position = "NONE",
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.25),
          panel.grid.major.y=element_line(color="grey80",linewidth=0.2),
          strip.text = element_text(size = 5,hjust = 0,vjust = 0,face = "bold",margin = margin(t = 2, b = 2),color = "black"),
          strip.background = element_rect(fill = "white",color = NA)
          
    )}

# =============================================================================
# Prepare LMER results joined with phenotype labels and breed lifespan
# =============================================================================
pheno <- raw_S2_ALL_PHENOS %>% select(phenotype,paper_phenotype_category,plot_label) 

# Keep only the lifespan effect from the ANOVA summary and store standardized effect sizes/p-values
lmer_anova <- pheno %>% right_join(raw_S9_ANOVA_ON_LMER %>% filter(Effect=="lifespan") %>% select(phenotype,Effect,ges=ges.anova,p.anova,beta_std.lm,p.lm=p.value.lm))

# LMER per-breed scores for each phenotype
lmer <- raw_S8_LMER %>% select(phenotype,breed,N,REML.t.val,ML.anova.p)
lmer <- lmer %>% filter(!is.na(REML.t.val))
lmer <- standardize_breed_names(lmer,breed)

# Attach phenotype metadata/statistics to each per-breed LMER record
lmer <- lmer_anova %>% inner_join(lmer)

# Breed-level lifespan metadata; filter to breeds present in LMER results
breeds <- raw_S12_BREED_PHENOS %>% select(breed,lifespan=`median_age_death.(McMillan.2024)`) %>% distinct() 
breeds <- standardize_breed_names(breeds,breed) %>% filter(breed %in% lmer$breed)
breeds <- convert_yesno_to_logical(breeds)

# Precompute plotting labels for breeds (line breaks/abbreviations for repelled labels)
breeds <- breeds %>% mutate(lmer_plot_breed = map_chr(breed, replace_middle_break))
breeds <- breeds %>% mutate(lmer_plot_breed=str_replace(lmer_plot_breed,"american","am."))
breeds <- breeds  %>% filter(!is.na(lifespan))
dim(lmer)

# Final joined dataset used for LMER scatterplots
lmer <- lmer %>% inner_join(breeds) %>% distinct()

# Quick check of strongest standardized effects across phenotypes
lmer_anova %>% arrange(-beta_std.lm)

# =============================================================================
# Function: per-phenotype scatter of breed ancestry score vs breed lifespan
# =============================================================================
plot_lmer_effect <- function(phenotype) {
  # Subset to one phenotype and mark statistically significant points
  pd <- lmer %>%
    filter(phenotype == !!phenotype) %>%
    mutate(sig = if_else(ML.anova.p <= 0.05, TRUE, FALSE))
  
  # Spearman correlation annotation (rank correlation between score and lifespan)
  cor_result <- cor.test(pd$REML.t.val, pd$lifespan, method = "spearman")
  r_val <- signif(cor_result$estimate, 2)
  p_val <- if (cor_result$p.value < 0.001) {
    format(cor_result$p.value, scientific = TRUE, digits = 2)
  } else {
    signif(cor_result$p.value, 2)
  }
  spearman_string <- paste0("rho = ", r_val, ", p = ", p_val)
  
  # Human-readable trait label for axis title
  in_plot_label <- pd %>%
    distinct(plot_label) %>%
    pull(plot_label) %>%
    first()
  
  # Label a limited set of top significant breeds (largest absolute t-values)
  top <- pd %>%
    filter(sig) %>%
    arrange(-abs(REML.t.val)) %>%
    slice_head(n = 10)
  
  # Space above the top point for the correlation text
  yrange <- max(pd$lifespan, na.rm = TRUE) - min(pd$lifespan, na.rm = TRUE)
  yshift <- yrange / 20
  
  # Scatter + linear trendline + repelled labels for key breeds
  p <- ggplot(pd, aes(x = REML.t.val, y = lifespan)) +
    geom_point(aes(color = sig, shape = sig), alpha = 0.75, size = 1) +
    geom_smooth(method = 'lm', fill = NA, color = "#a50f15", alpha = 0.5, linewidth = 0.5) +
    annotate("text", label = spearman_string,
             x = min( pd$REML.t.val, na.rm = TRUE),
             y = max(pd$lifespan, na.rm = TRUE) + yshift * 2, size = 2, vjust = 1, hjust = 0, color = "#a50f15") +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.25, linetype = 2) +
    geom_text_repel(
      data = top,
      aes(label = lmer_plot_breed),
      force = 5,
      color = "#a50f15",
      segment.color = "black",
      box.padding = 0.5,
      segment.size = 0.25,
      point.padding = 0.3,
      size = 1.5,
      lineheight = 0.8,
      min.segment.length = 0,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "#a50f15")) +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 16)) +
    scale_y_continuous("breed median lifespan (years)",
                       expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_continuous(paste0("Breed Ancestry Score on ",in_plot_label))+
    theme_lmer()
  
  return(p)
}

# =============================================================================
# Load longitudinal (Kaplan-Meier) datasets if not already in memory
# =============================================================================
if (!exists("raw_metabolite_ldat")  | !exists("raw_cbc_chem_ldat")){
  load("../../data/kaplan_meier/longitudinalData2")
  raw_metabolite_ldat <- as_tibble(ldat)
  rm(ldat)
  load("../../data/kaplan_meier/longitudinalData2.addingCBC_chem_data5")
  raw_cbc_chem_ldat <- as_tibble(ldat %>% group_by(dog_id) %>% filter(tstop == max(tstop)) %>% ungroup())
  rm(ldat)
}

# =============================================================================
# Theme for Kaplan-Meier survival panels
# =============================================================================
theme_km <- function(){
  theme_cowplot(6) %+replace%
    theme(
      plot.title = element_text(size = 6.5,face="bold",hjust=0,margin=margin(b=4)),
      axis.ticks = element_line(linewidth=0.25),
      axis.line = element_line(linewidth = 0.25),    
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.position = "none",  
      panel.grid.major= element_line(linewidth = 0.2,color="grey80"),
      strip.text = element_text(size = 6, hjust = 0.5, vjust = 0, face = "bold",margin = margin(t = 3, b = 3)),
      strip.background = element_rect(fill = "white", color = NA)
    )
}

# =============================================================================
# Reshape longitudinal data to long format and harmonize phenotype naming
# =============================================================================
cbc_chem_ldat <- raw_cbc_chem_ldat %>% pivot_longer(c(-dog_id,-recruitmentAge,-sex,-Weight_Class_10KGBin_at_HLES,-tstart,-tstop,-sterile,-death)) %>% 
  rename(km_phenotype=name)

# Drop character krt* columns before pivoting (avoids mixing character with numeric measures)
metabolite_ldat <- raw_metabolite_ldat %>%
  select(-any_of(names(.)[grepl("^krt", names(.)) & vapply(., is.character, logical(1))]))
metabolite_ldat <- metabolite_ldat  %>% pivot_longer(c(-dog_id,-recruitmentAge,-sex,-sterilization,-Weight_Class_10KGBin_at_HLES,-tstart,-tstop,-death,-genomeSex,-coi,-size,-bestWt)) %>% 
  rename(km_phenotype=name)

# Combine the two longitudinal sources into one standardized long dataset
ldat <- metabolite_ldat %>% mutate(set="metabolite") %>% bind_rows(cbc_chem_ldat %>% rename(sterilization=sterile) %>% mutate(set="analyte"))

# Clean phenotype strings to match join keys (replace punctuation/spaces; remove parentheses)
ldat <- ldat %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\/", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, " ", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\,", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\(", "")) %>% 
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\)", ""))

# Create a matching KM phenotype key from the master phenotype table
pheno <- raw_S2_ALL_PHENOS %>% mutate(km_phenotype=phenotype) %>% 
  mutate(km_phenotype=gsub("_ln_transformed", "", km_phenotype)) %>% 
  mutate(km_phenotype=gsub("_sqrt_transformed", "", km_phenotype)) 

# Direction of association with lifespan for each phenotype (used to define "risk" tail)
traits <- raw_S11_COMBINED %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype,plot_label,paper_phenotype_category)) %>% 
  filter(Effect=="lifespan") %>% 
  mutate(direc=if_else(z_comb>0,"risk","protective")) %>% 
  select(phenotype,direc) %>% distinct()

# Join in phenotype metadata and keep only records with measurements and known direction
ldat <- pheno %>% inner_join(ldat) %>% filter(!is.na(value)) %>% inner_join(traits)

# =============================================================================
# Define recruitment age category and restrict to the last observation per dog
# =============================================================================
ldat <- ldat %>% mutate(recruitAge_categ = case_when(recruitmentAge <= 7 ~ "not_old",
                                                     recruitmentAge > 7 ~ "old"))
ldat$recruitAge_categ <- factor(ldat$recruitAge_categ, levels = c("not_old", "old"))
table(ldat$recruitAge_categ)

# Keep only each dog's final timepoint (KM curves use tstop and event status)
ldat <- ldat %>% group_by(dog_id) %>% filter(tstop == max(tstop)) %>% ungroup()

# Recompute age category after filtering (ensures factor levels are consistent)
ldat <- ldat %>% mutate(recruitAge_categ = case_when(recruitmentAge <= 7 ~ "not_old",recruitmentAge > 7 ~ "old"))
ldat$recruitAge_categ <- factor(ldat$recruitAge_categ, levels = c("not_old", "old"))

# =============================================================================
# Function: Kaplan-Meier plot for one phenotype (optionally stratified by age)
# =============================================================================
plot_one_kaplan <- function(in_phenotype,with_age=FALSE){
  pd <- ldat %>% filter(phenotype==in_phenotype) 
  risk_direction <- pd %>% select(direc) %>% distinct() %>% pull(direc)
  title <- pd %>% select(plot_label) %>% distinct() %>% pull(plot_label)
  
  # Compute within-phenotype quartiles to define risk vs not-risk groups
  cutoffs <- pd %>% group_by(km_phenotype) %>% 
    summarize(q25 = quantile(value, 0.25, na.rm = TRUE),
              q75 = quantile(value, 0.75, na.rm = TRUE))
  
  pd <- pd %>% left_join(cutoffs)
  if (risk_direction == "risk") {
    pd <- pd %>% mutate(risk_categ = if_else(value > q75,"risk","not_risk"))
  } else {
    pd <- pd %>% mutate(risk_categ = if_else(value < q25,"risk","not_risk"))
  }
  pd$risk_categ <- factor(pd$risk_categ, levels = c("not_risk","risk"))
  
  # Combined stratification label when splitting by both risk group and recruitment age
  pd <- pd %>% mutate(risk_age_categ = interaction(risk_categ, recruitAge_categ, sep = "_"))
  
  # Fit KM curves by risk group
  risk_fits <- pd %>%
    nest() %>%
    mutate(
      fit = map(data, ~ survfit(
        Surv(tstop, death == "Dead") ~ risk_categ,
        data = .x
      ))
    )
  
  # Convert survival fit into a tidy data frame for ggplot
  surv_df <- risk_fits %>%
    select(fit) %>%
    mutate(curve = map(fit, ~ broom::tidy(.x))) %>%
    unnest(curve) %>%
    mutate(strata = str_remove(strata, "^risk_categ="))
  
  # Optional: refit and plot curves stratified by risk group x age group
  if (with_age) {
    
    risk_fits <- pd %>%
      nest() %>%
      mutate(
        fit = map(data, ~ survfit(
          Surv(tstop, death == "Dead") ~ risk_age_categ,
          data = .x
        ))
      )
    
    surv_df <- risk_fits %>%
      select(fit) %>%
      mutate(curve = map(fit, ~ broom::tidy(.x))) %>%
      unnest(curve) %>%
      mutate(strata = str_remove(strata, "^risk_age_categ="))
  }
  
  # Human-friendly curve labels depend on whether "risk" corresponds to high or low values
  labels <- pd %>% select(strata=risk_categ) %>% distinct()
  labels <- labels %>% mutate(label=if_else(strata=="risk"&risk_direction=="risk","High",
                                            if_else(strata=="risk"&risk_direction=="protective","Low",
                                                    if_else(strata=="not_risk"&risk_direction=="risk","Not high",
                                                            if_else(strata=="not_risk"&risk_direction=="protective","Not low","error")))))
  
  # Expand labels to match age-stratified strata naming when with_age=TRUE
  labels2 <- labels %>% filter(!is.na(label)) %>% mutate(strata=paste0(strata,"_not_old"),label=paste0(label," (not old)")) 
  labels2 <- labels %>% mutate(strata=paste0(strata,"_old"),label=paste0(label," (old)")) %>% bind_rows(labels2)
  labels <- labels %>% bind_rows(labels2) %>% rename(gglab=label) %>% mutate(gglab=str_replace(gglab," \\(","\n("))
  
  # Build right-side line-end labels: compute curve endpoints and spread labels vertically to avoid overlap
  line_end <- surv_df %>% group_by(strata) %>% summarize(time=max(time)) %>% 
    inner_join(surv_df,relationship = "many-to-many") %>% left_join(labels)
  x_lab <- max(surv_df$time, na.rm = TRUE) * 1.06
  
  line_end <- line_end %>% mutate(yrank=rank(estimate,ties.method = "first"),yspan=max(estimate)-min(estimate),ymidpt=mean(estimate)) %>% ungroup()
  line_end <- line_end %>%
    mutate(incr=max(yspan/3,0.12)) %>% 
    mutate(
      y_labelpos_raw = ymidpt + (yrank - (max(yrank) + 1) / 2) * incr,
      shift_down = pmax(0, max(y_labelpos_raw, na.rm = TRUE) - 0.96),
      y_labelpos = y_labelpos_raw - shift_down
    ) %>%
    ungroup() %>%
    select(-y_labelpos_raw, -shift_down)
  
  # KM plot with confidence ribbons and right-side labels
  p <- ggplot(surv_df,aes(x = time, y = estimate, color = strata, fill = strata))
  p <- p + geom_step(linewidth = 0.4)
  p <- p + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA)
  p <- p + geom_text(data=line_end,aes(y=y_labelpos,x=x_lab,label = gglab,color = strata),
                     hjust=0,vjust=0.5,lineheight=0.9,size=2)
  p <- p + geom_segment(data=line_end,aes(x=time+0.05,xend=x_lab-0.05,y=estimate,yend=y_labelpos,color = strata),size=0.25)
  p <- p + scale_x_continuous("year",breaks=c(0:4),limits=c(0,5),expand = expansion(mult = c(0, 0)))
  p <- p + scale_y_continuous("survival probability",breaks=c(0,0.25,0.50,0.75,1),limits=c(0,1.05),expand = expansion(mult = c(0, 0)))
  p <- p + ggtitle(title)
  p <- p + scale_color_manual(values = c("not_risk" = "#2166ac", "risk" = "#b2182b","not_risk_old" = "#2166ac", "risk_old" = "#b2182b",  risk_not_old = "#762a83",not_risk_not_old = "#000000"), drop = TRUE)
  p <- p + scale_fill_manual(values = c("not_risk" = "#2166ac", "risk" = "#b2182b","not_risk_old" = "#2166ac", "risk_old" = "#b2182b",  risk_not_old = "#762a83",not_risk_not_old = "#000000"), drop = TRUE)
  p <- p + theme_km() 
  p
}

# =============================================================================
# Build panels for the two highlighted traits (Globulins and Potassium)
# =============================================================================
p_pot_km <- plot_one_kaplan("krt_cp_potassium_value_ln_transformed",TRUE)
p_glob_km <- plot_one_kaplan("krt_cp_globulins_value_ln_transformed",TRUE)

dot_plot1 <- plot_lmer_effect(phenotype = "krt_cp_globulins_value_ln_transformed")

# Breed icon grids for low vs high ends of the breed ancestry score distribution
breed_plot1 <- make_breed_figure_custom_lists("krt_cp_globulins_value_ln_transformed", low_breeds=c("labrador_retriever", "australian_shepherd", "english_setter"),
                                                      high_breeds=c("great_dane","english_bulldog","rottweiler"),
                                                      max_px=400,font_size=4)
dot_plot2 <- plot_lmer_effect(phenotype = "krt_cp_potassium_value_ln_transformed")
breed_plot2 <- make_breed_figure_custom_lists("krt_cp_potassium_value_ln_transformed", low_breeds=c("wire_fox_terrier","border_collie","belgian_malinois"),
                                                      high_breeds=c("airedale_terrier","pug","english_bulldog"),
                                                      max_px=400,font_size=4)

# Tweak margins so the two breed grids sit tightly as a combined row
breed_plot1$low <- breed_plot1$low + theme(plot.margin = margin(l=0,r=5,t=-5,b=10))
breed_plot1$high <- breed_plot1$high + theme(plot.margin = margin(l=5,r=0,t=-5,b=10))
breed_plot1_dogs <- plot_grid(breed_plot1$low ,breed_plot1$high,nrow=1,align = "v",axis = "b")

breed_plot2$low <- breed_plot2$low + theme(plot.margin = margin(l=0,r=5,t=-5,b=5))
breed_plot2$high <- breed_plot2$high + theme(plot.margin = margin(l=5,r=0,t=-5,b=5))
breed_plot2_dogs <- plot_grid(breed_plot2$low ,breed_plot2$high,nrow=1,align = "v",axis = "b")

# =============================================================================
# Assemble final multi-row figure and save to PDF
# =============================================================================
row1 <- plot_grid(dot_plot1,dot_plot2,nrow=1,label_size=8,labels=c("A","B"))
row2 <- plot_grid(breed_plot1_dogs,breed_plot2_dogs,nrow=1,align = "h", axis = "b")
row3 <- plot_grid(p_glob_km,p_pot_km,nrow=1,label_size=8,labels=c("C","D"))
combined <- plot_grid(row1,row2,NULL,row3,ncol=1,rel_heights = c(1,0.4,0.05,1))
ggsave(plot=combined,filename=outpdf,width=5,height=4.25)
