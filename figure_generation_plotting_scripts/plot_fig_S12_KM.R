# -------------------------------------------------------------------------
# Kaplan–Meier survival curves for lifespan-associated phenotypes
# - Load phenotype + longitudinal data (metabolites and clinical analytes)
# - Define "risk" vs "not_risk" groups using phenotype-specific quartiles
# - Fit KM curves and assemble multi-panel figure
# -------------------------------------------------------------------------

library(tidyverse)
library(survival)
library(stringr)
library(cowplot)
library(broom)

# Output paths (and unused sheet URL retained as a record of provenance)
supp_sheet_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/edit?usp=sharing"#
outpdf <- "../../../figures_pdf/fig_S12_KM.pdf"
outRdata_glob_pot <- "./fig_KM.for_fig6.Rdata"

# -------------------------------------------------------------------------
# Load prepackaged supplemental objects into the global environment (if needed)
# -------------------------------------------------------------------------
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_S2_ALL_PHENOS, raw_S10_COX_HAZ, raw_S11_COMBINED, etc.
}

# -------------------------------------------------------------------------
# Load longitudinal datasets (if not already present in the session)
# - Metabolite longitudinal data
# - Clinical analyte CBC/chemistry data (restricted to each dog's last timepoint)
# -------------------------------------------------------------------------
if (!exists("raw_metabolite_ldat")  | !exists("raw_cbc_chem_ldat")){
  load("../data/kaplan_meier/longitudinalData2")
  raw_metabolite_ldat <- as_tibble(ldat)
  rm(ldat)
  load("../data/kaplan_meier/addingCBC_chem_forVista_data5")
  raw_cbc_chem_ldat <- as_tibble(ldat %>% group_by(dog_id) %>% filter(tstop == max(tstop)) %>% ungroup())
  rm(ldat)
}

# -------------------------------------------------------------------------
# Reshape both datasets to long format:
# - One row per dog per phenotype (km_phenotype) with a single `value`
# -------------------------------------------------------------------------
cbc_chem_ldat <- raw_cbc_chem_ldat %>% pivot_longer(c(-dog_id,-recruitmentAge,-sex,-Weight_Class_10KGBin_at_HLES,-tstart,-tstop,-sterile,-death)) %>% 
  rename(km_phenotype=name)

# Drop character krt* columns (these can interfere with pivoting / downstream joins)
metabolite_ldat <- raw_metabolite_ldat %>%
  select(-any_of(names(.)[grepl("^krt", names(.)) & vapply(., is.character, logical(1))]))

metabolite_ldat <- metabolite_ldat  %>% pivot_longer(c(-dog_id,-recruitmentAge,-sex,-sterilization,-Weight_Class_10KGBin_at_HLES,-tstart,-tstop,-death,-genomeSex,-coi,-size,-bestWt)) %>% 
  rename(km_phenotype=name)

# Combine into a single long dataset; align sterilization column naming
ldat <- metabolite_ldat %>% mutate(set="metabolite") %>% bind_rows(cbc_chem_ldat %>% rename(sterilization=sterile) %>% mutate(set="analyte"))

# Clean phenotype names so they are safe for matching/plotting
ldat <- ldat %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\/", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, " ", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\,", "-")) %>%
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\(", "")) %>% 
  mutate(km_phenotype = str_replace_all(km_phenotype, "\\)", ""))

# -------------------------------------------------------------------------
# Build phenotype mapping table used to join plot labels and standardize names
# -------------------------------------------------------------------------
pheno <- raw_S2_ALL_PHENOS %>% mutate(km_phenotype=phenotype) %>% 
  mutate(km_phenotype=gsub("_ln_transformed", "", km_phenotype)) %>% 
  mutate(km_phenotype=gsub("_sqrt_transformed", "", km_phenotype)) 

# Join phenotype metadata, drop missing values, and keep the last timepoint per dog
ldat <- pheno %>% inner_join(ldat) %>% filter(!is.na(value))
ldat <- ldat %>% group_by(dog_id) %>% filter(tstop == max(tstop)) %>% ungroup()

# -------------------------------------------------------------------------
# Define age category at recruitment (used optionally for stratification)
# -------------------------------------------------------------------------
ldat <- ldat %>% mutate(recruitAge_categ = case_when(recruitmentAge <= 7 ~ "not_old",recruitmentAge > 7 ~ "old"))
ldat$recruitAge_categ <- factor(ldat$recruitAge_categ, levels = c("not_old", "old"))

# -------------------------------------------------------------------------
# Identify traits to plot:
# - Significant in combined analysis and/or Cox model (FDR <= 0.05)
# - Assign direction (risk vs protective) using HR or combined z
# -------------------------------------------------------------------------
traits_comb <- raw_S11_COMBINED %>% 
  filter(Effect=="lifespan") %>% 
  select(phenotype,z_comb,p_comb_fdr)

traits_cox <- raw_S10_COX_HAZ %>% select(phenotype,cox_p_adj_fdr=p_adj_fdr,cox_HR=`exp(coef)`)

traits <-   raw_S2_ALL_PHENOS %>% select(phenotype,plot_label,paper_phenotype_category) %>% 
  inner_join(traits_comb %>% full_join(traits_cox))

traits <- traits %>%
  filter(!is.na(cox_HR)|!is.na(z_comb)) %>% 
  filter(p_comb_fdr<=0.05|cox_p_adj_fdr<=0.05) %>% 
  mutate(direc=if_else(!is.na(cox_HR),if_else(cox_HR>1,"risk","protective"),
                       if_else(z_comb>0,"risk","protective")))

# Dataset used for plotting: one row per dog per phenotype with direction and survival info
plotdata <- ldat %>% inner_join(traits)
plotdata <- plotdata %>% select(phenotype,km_phenotype,direc,dog_id,recruitAge_categ,tstop,death,value)

# -------------------------------------------------------------------------
# Plot a single phenotype:
# - Define risk group using phenotype-specific quartiles (Q75 for risk traits; Q25 for protective traits)
# - Fit KM curve by risk group (or risk-by-age group if with_age=TRUE)
# - Add end-of-line labels with non-overlapping vertical positioning
# -------------------------------------------------------------------------
plot_kaplan_one <- function(phenotype_name,with_age = FALSE) {
  
  # Filter to a single phenotype
  pd <- plotdata %>%
    filter(phenotype == phenotype_name)
  risk_direction <- pd %>% pull(direc) %>% first()  
  
  if (nrow(pd) == 0) {
    warning(paste0("No rows after filtering for phenotype='", phenotype_name),
            call. = FALSE)
    return(NULL)
  }
  
  # Compute per-phenotype quartile cutoffs
  cutoffs <- pd %>%
    group_by(km_phenotype) %>%
    summarize(
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  pd <- pd %>% left_join(cutoffs, by = "km_phenotype")
  
  # Risk group assignment depends on whether higher values are harmful (risk) or beneficial (protective)
  if (risk_direction == "risk") {
    pd <- pd %>% mutate(risk_categ = if_else(value > q75, "risk", "not_risk"))
  } else {
    pd <- pd %>% mutate(risk_categ = if_else(value < q25, "risk", "not_risk"))
  }
  
  pd$risk_categ <- factor(pd$risk_categ, levels = c("not_risk", "risk"))
  pd <- pd %>% mutate(risk_age_categ = interaction(risk_categ, recruitAge_categ, sep = "_"))
  
  # Fit KM curves, optionally stratified by age category as well
  if (!with_age) {
    fit <- survfit(Surv(tstop, death == "Dead") ~ risk_categ, data = pd)
    
    surv_df <- broom::tidy(fit) %>%
      mutate(strata = stringr::str_remove(strata, "^risk_categ=")) %>%
      mutate(km_phenotype = pd$km_phenotype[1])
  } else {
    fit <- survfit(Surv(tstop, death == "Dead") ~ risk_age_categ, data = pd)
    
    surv_df <- broom::tidy(fit) %>%
      mutate(strata = stringr::str_remove(strata, "^risk_age_categ=")) %>%
      mutate(km_phenotype = pd$km_phenotype[1])
  }
  
  # Build human-readable legend/label text for strata, customized by risk direction
  labels <- pd %>%
    select(strata = risk_categ) %>%
    distinct()
  
  labels <- labels %>%
    mutate(label = if_else(strata == "risk" & risk_direction == "risk", "High",
                           if_else(strata == "risk" & risk_direction == "protective", "Low",
                                   if_else(strata == "not_risk" & risk_direction == "risk", "Not high",
                                           if_else(strata == "not_risk" & risk_direction == "protective", "Not low", "error")))))
  
  # Precompute age-specific variants of labels (used when plotting with age strata)
  labels2_young <- labels %>%
    filter(!is.na(label)) %>%
    mutate(strata = paste0(strata, "_not_old"),
           label  = paste0(label, " (young)"))
  
  labels2_old <- labels %>%
    mutate(strata = paste0(strata, "_old"),
           label  = paste0(label, " (old)"))
  
  labels <- labels %>%
    bind_rows(labels2_old, labels2_young) %>%
    rename(gglab = label) %>%
    mutate(gglab = stringr::str_replace(gglab, " \\(", "\n("))
  
  # Join plot labels for title text (via km_phenotype key)
  surv_df <- pheno %>%
    select(km_phenotype, phenotype, plot_label) %>%
    right_join(surv_df, by = "km_phenotype")
  
  # Find last point per curve to place labels and connecting segments
  line_end <- surv_df %>%
    group_by(km_phenotype, strata) %>%
    summarize(time = max(time), .groups = "drop") %>%
    inner_join(surv_df, by = c("km_phenotype", "strata", "time")) %>%
    left_join(labels, by = "strata")
  
  # Label x-position placed slightly beyond the max follow-up time
  x_lab <- max(surv_df$time, na.rm = TRUE) * 1.06
  
  # Compute non-overlapping y positions for end labels, with bounds protection
  line_end <- line_end %>%
    group_by(km_phenotype) %>%
    mutate(
      yrank = rank(estimate, ties.method = "first"),
      yspan = max(estimate) - min(estimate),
      ymidpt = mean(estimate)
    ) %>%
    ungroup()
  
  max_spacing <- if (with_age) 0.25 else 0.2
  
  line_end <- line_end %>%
    group_by(km_phenotype) %>%
    mutate(
      incr = max(yspan / 3, max_spacing),
      y_labelpos_raw = ymidpt + (yrank - (max(yrank) + 1) / 2) * incr,
      shift_down = pmax(0, max(y_labelpos_raw, na.rm = TRUE) - 0.96),
      shift_up   = pmax(0, 0.1 - min(y_labelpos_raw, na.rm = TRUE)),
      shift      = shift_up - shift_down,
      y_labelpos = y_labelpos_raw + shift
    ) %>%
    ungroup() %>%
    select(-y_labelpos_raw, -shift_down, -shift_up, -shift)
  
  # Titles/subtitles for the panel
  title <- str_to_sentence(traits %>% filter(phenotype == phenotype_name) %>%
    pull(plot_label) %>% first() ) 
  subtitle <- traits %>% filter(phenotype == phenotype_name) %>%
    pull(subtitle) %>% first()  
  
  # Plot KM curve with confidence ribbon and end-of-line labels
  p <- ggplot(surv_df, aes(x = time, y = estimate, color = strata, fill = strata)) +
    geom_step(linewidth = 0.4) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
    geom_text(
      data = line_end,
      aes(y = y_labelpos, x = x_lab, label = gglab, color = strata),
      hjust = 0, vjust = 0.5, lineheight = 0.75, size = 1.5
    ) +
    geom_segment(
      data = line_end,
      aes(x = time + 0.05, xend = x_lab - 0.05, y = estimate, yend = y_labelpos, color = strata),
      linewidth = 0.15
    ) +
    scale_x_continuous("year", breaks = c(0:4), limits = c(0, 5), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous("survival probability",
                       breaks = c(0, 0.25, 0.50, 0.75, 1),
                       limits = c(0, 1.05),
                       expand = expansion(mult = c(0.05, 0.05))) +
    ggtitle(title, subtitle = subtitle) +
    theme_cowplot(6) +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 5, hjust = 0, margin = margin(b = 2)),
      axis.ticks = element_line(linewidth = 0.25),
      axis.line = element_line(linewidth = 0.25),
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 5),
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.2, color = "grey80"),
      strip.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(2, "mm")
    )
  
  # Manual palette for risk strata (kept identical to original)
  p <- p +
    scale_color_manual(
      values = c(
        "not_risk" = "#2166ac", "risk" = "#b2182b",
        "not_risk_old" = "#2166ac", "risk_old" = "#b2182b",
        "risk_not_old" = "#762a83", "not_risk_not_old" = "#000000"
      ),
      drop = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "not_risk" = "#2166ac", "risk" = "#b2182b",
        "not_risk_old" = "#2166ac", "risk_old" = "#b2182b",
        "risk_not_old" = "#762a83", "not_risk_not_old" = "#000000"
      ),
      drop = TRUE
    )
  
  p
}

# -------------------------------------------------------------------------
# Helper for formatting p-values in subtitles
# -------------------------------------------------------------------------
format_p <- function(p) {
  case_when(
    is.na(p)          ~ "n.s.",
    p >= 0.05         ~ "n.s.",
    p < 1e-3          ~ formatC(p, format = "e", digits = 1),
    TRUE              ~ as.character(signif(p, 2))
  )
}

# -------------------------------------------------------------------------
# Construct per-phenotype subtitle text combining direction, Cox, and combined results
# -------------------------------------------------------------------------
traits <- traits %>% mutate(subtitle=str_to_sentence(paste0(direc," ",paper_phenotype_category))) 
traits <- traits %>%
  mutate(
    combined_label = if_else(
      format_p(p_comb_fdr) == "n.s.",
      "Combined pFDR = n.s.",
      paste0("Combined z = ",round(z_comb,2)," (pFDR = ", format_p(p_comb_fdr),")")
    ),
    cox_p_label = if_else(
      format_p(cox_p_adj_fdr) == "n.s.",
      "Cox pFDR = n.s.",
      paste0("Cox HR = ", round(cox_HR,2)," (pFDR = ", format_p(cox_p_adj_fdr),")")
    )
  )

traits <- traits %>% mutate(subtitle=paste(subtitle,cox_p_label,combined_label,sep="\n")) 

# Phenotype ordering for two panels:
# A: significant in combined analysis; B: Cox-only significant among clinical analytes
phenos_vec_A <- traits %>% filter(p_comb_fdr<=0.05) %>% arrange(-cox_HR,p_comb_fdr) %>% pull(phenotype)
phenos_vec_B <- traits %>% filter(p_comb_fdr>0.05&cox_p_adj_fdr<=0.05&paper_phenotype_category=="Clinical analytes") %>% arrange(-cox_HR,p_comb_fdr) %>% pull(phenotype)

# Build one plot per phenotype (in order)
plots_A <- map(phenos_vec_A, ~ plot_kaplan_one(phenotype_name = .x, with_age = FALSE))
plots_B <- map(phenos_vec_B, ~ plot_kaplan_one(phenotype_name = .x, with_age = FALSE))

# Drop any NULLs (e.g., phenotypes with no data) while preserving order
plots_A <- compact(plots_A)
plots_B <- compact(plots_B)

if (length(plots_A) == 0) stop("plots_A is empty (all NULL).")
if (length(plots_B) == 0) stop("plots_B is empty (all NULL).")

# -------------------------------------------------------------------------
# Assemble grids and save figure
# -------------------------------------------------------------------------
in_ncol <- 3
in_nrow_A <- ceiling(length(phenos_vec_A)/in_ncol)
in_nrow_B <- ceiling(length(phenos_vec_B)/in_ncol)

p_grid_A <- plot_grid(
  plotlist = plots_A,
  ncol = in_ncol,
  nrow = in_nrow_A,
  align = "hv",
  axis = "tblr"
)

p_grid_B <- plot_grid(
  plotlist = plots_B,
  ncol = in_ncol,
  nrow = in_nrow_B,
  align = "hv",
  axis = "tblr"

)

grid <- plot_grid(p_grid_A,p_grid_B,ncol=1,rel_heights = c(in_nrow_A,in_nrow_B),  labels=LETTERS,label_size = 8)

ggsave(plot=grid,filename=outpdf,width=6,height=8)
