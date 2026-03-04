library(tidyverse)
library(cowplot)
library(googlesheets4)
source("util.cleaning.R")

# Output path for the assembled Figure 4 (FUR) PDF
outpdf <- "../../../figures_pdf/fig_4_FUR.pdf"

# =====================================
# Load all data objects (if not already in the environment)
# =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects such as raw_S3_FIFTY_FUR, raw_S15_BREED_ANC, pop, etc.
}

# =====================================
# Import GWAS locus zoom panels from PDF and arrange as a single A/B row
# =====================================
p_fur_1 <- pdf_to_plot("source_figures/Cystathionine.locus_zoom_chr13.pdf", scale = 1)
p_fur_2 <- pdf_to_plot("source_figures/Cystathionine.locus_zoom_chr32.pdf", scale = 1)

# Helper: add consistent plot margins so the PDF-to-ggplot images align cleanly
pad_plot <- function(p, margin = 6) {
  ggdraw(p) +
    theme(plot.margin = margin(margin, margin, margin, margin))
}

p_fur_1p <- pad_plot(p_fur_1, margin = 6)
p_fur_2p <- pad_plot(p_fur_2, margin = 6)

# Add a small bottom spacer to standardize heights between the two imported PDF panels
p_fur_1b <- plot_grid(p_fur_1p, NULL, ncol = 1, rel_heights = c(1, 0.02))
p_fur_2b <- plot_grid(p_fur_2p, NULL, ncol = 1, rel_heights = c(1, 0.02))

# Combine locus zoom panels side-by-side with labels A/B
p_fur_gwas <- plot_grid(
  p_fur_1b, p_fur_2b,
  nrow = 1,
  align = "hv",
  axis = "b",
  rel_widths = c(1, 1.01),
  labels = c("A", "B"),
  label_size = 8,
  label_x = c(0.02, 0.02),
  label_y = c(0.98, 0.98),
  hjust = 0,
  vjust = 1
)

# =====================================
# Breed-frequency (MAF) boxplots for selected SNPs by phenotype category
# =====================================
make_maf_boxplot_snps <- function(in_snps, in_trait, print_title = FALSE) {
  
  # Human-readable trait labels
  print_trait <- case_when(in_trait == "length" ~ "Fur length", 
                           in_trait == "furnishings" ~ "Furnishings", TRUE ~ in_trait)
  # Lazy-load breed phenotypes and SNP breed frequency tables from the supplement sheet
  if (!exists("raw_S12_BREED_PHENOS")) {
    raw_S12_BREED_PHENOS <- read_sheet(supp_sheet_url, range = "Data S12_BREED_PHENOS", na = c("","NA","#N/A"), skip = 0) %>% 
      as_tibble() %>% rename_with(~ gsub("[?'~ /]", "_", .x)) 
  }
  if (!exists("raw_S7_CYST_BREED_FREQ")) { 
    raw_S7_CYST_BREED_FREQ <- read_sheet(supp_sheet_url, range = "Data S7_CYST_BREED_FREQ", na = c("","NA","#N/A"), skip = 0) %>% 
      as_tibble() %>% rename_with(~ gsub("[?'~ /]", "_", .x))
  }
  
  # One row per breed with its trait categories (furnishings and length)
  maf_breed_phenos <- raw_S12_BREED_PHENOS %>% select(breed, furnishings, length) %>% distinct()
  
  # Merge SNP frequency table to breed phenotype categories and reshape long by trait
  maf_cyst_breed <- raw_S7_CYST_BREED_FREQ %>% mutate(breed = str_replace_all(CLST, "_", " ")) %>% 
    inner_join(maf_breed_phenos, by = "breed") %>% select(-CLST) %>% 
    pivot_longer(c(furnishings, length), names_to = "trait", values_to = "phenotype") %>% filter(phenotype!="variable")
  
  # Filter to requested SNPs and trait
  maf_pd <- maf_cyst_breed %>% filter(trait == in_trait & SNP %in% in_snps)
  
  # Extract position from "chr:pos:ref:alt" and create Mb label for x axis
  maf_pd <- maf_pd %>% mutate(pos = as.integer(str_extract(SNP, "(?<=:)[0-9]+"))) %>% 
    mutate(pos_str=as.character(round(pos/1000000,3)))
  
  # Flip allele frequency for specific SNPs so the same (effect) allele is plotted consistently
  maf_pd <- maf_pd %>% mutate(MAF = if_else(SNP %in% c("chr32:35494936:C:T","chr32:35496522:G:A"), 1 - MAF, MAF))
  
  # Standardize phenotype labels for display
  maf_pd <- maf_pd %>% mutate(phenotype=if_else(phenotype=="no","no furn.",if_else(phenotype=="yes","furnishings",phenotype)))
  
  # Build ordered facet labels including per-facet N
  order <- tibble(phenotype = c("short","no furn.","medium","long","furnishings"), order = 1:5) %>% inner_join(maf_pd %>% distinct(phenotype), by = "phenotype") %>% arrange(order)
  order <- maf_pd %>% select(trait,phenotype,breed) %>% distinct() %>% group_by(trait,phenotype) %>% count() %>% full_join(order)
  order <- order %>% mutate(pheno_name=paste0(phenotype,"\n(N=",n,")")) %>% ungroup() %>% arrange(order)
  maf_pd <- order %>% select(pheno_name,phenotype) %>% right_join(maf_pd)
  maf_pd$phenotype <- factor(maf_pd$phenotype, levels = order$phenotype)
  maf_pd$pheno_name <- factor(maf_pd$pheno_name, levels = order$pheno_name)
  
  # Per-SNP permutation-based independence test of MAF vs phenotype groups (approximate distribution)
  perm_p <- maf_pd %>% group_by(SNP) %>% group_modify(~{ tmp <- .x; if (n_distinct(tmp$phenotype) < 2 || nrow(tmp) < 2) return(tibble(p = NA_real_)); tt <- coin::independence_test(MAF ~ phenotype, data = tmp, distribution = "approximate"); tibble(p = as.numeric(coin::pvalue(tt))) }) %>% ungroup()
  maf_pd <- maf_pd %>% left_join(perm_p, by = "SNP") %>% mutate(formatted_p = case_when(is.na(p) ~ "p = NA", p < 0.0001 ~ "p < 0.0001", TRUE ~ paste0("p = ", signif(p, 2))))
  xlabels <- maf_pd %>% group_by(SNP, phenotype) %>% summarize(n = n(), mean = mean(MAF), sd = sd(MAF), .groups = "drop") %>% mutate(label = paste0(phenotype, "\n(", n, ")"))
  
  # Draw boxplots of breed frequencies for the selected SNPs, facetted by phenotype category
  p <- ggplot(maf_pd, aes(x = pos_str, y = MAF))
  p <- p + geom_boxplot(fill = "grey90", color = "grey30", outlier.shape = NA, width = 0.75, linewidth = 0.25, alpha = 0.25)
  p <- p + geom_jitter(fill = "grey90", color = "grey30", shape = 16, width = 0.1, height = 0, alpha = 0.5)
  p <- p + scale_x_discrete(name = "position (Mb)")
  p <- p + scale_y_continuous("breed freq.", breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.05)))
  p <- p + facet_wrap(~ pheno_name, scales = "fixed")
  p <- p + ggtitle(paste0("chr. ",unique(maf_pd$CHR)))
  p <- p + labs(subtitle = NULL)
  p <- p + theme_cowplot(6)
  p <- p + theme(plot.title = element_text(size = 6.5, face = "bold",
                                           margin = margin(b = 1,t=4)), 
                 plot.subtitle = element_text(size = 5, margin = margin(t = 0, b = 2)), 
                 axis.line.y = element_line(linewidth = 0.25), 
                 axis.line.x = element_blank(), 
                 axis.ticks.y = element_line(linewidth = 0.25),
                 axis.ticks.x = element_blank(), axis.text.y = element_text(size = 5), 
                 axis.text.x = element_text(size = 5,angle=90,vjust=1,hjust=0.5), 
                 axis.title.y = element_text(size = 6, face = "bold"), 
                 axis.title.x = element_blank(), 
                 strip.text = element_text(size = 5, hjust = 0.5, vjust = 0, face = "bold", margin = margin(t = 2, b = 5), color = "black"), strip.background = element_rect(fill = "white", color = NA), panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25), legend.position = "none")
  p
}

# =====================================
# Build per-dog top ancestry assignment (used later for fur/trait context)
# =====================================
anc<- standardize_breed_names(raw_S15_BREED_ANC, pop)

anc_hair <- anc %>%
  group_by(dog_id) %>%
  summarize(pct = max(pct), .groups = "drop") %>%
  inner_join(anc, by = c("dog_id", "pct")) %>%
  rename(top_pop = pop)

# =====================================
# Prepare fur phenotype + cystathionine data for plotting
# =====================================
d <- raw_S3_FIFTY_FUR %>%
  mutate(dog_id = as.numeric(dog_id)) %>%
  left_join(anc_hair, by = "dog_id")

# Reshape trait columns long for faceting across furnishings vs fur length
df <- d %>%
  select(dog_id, Furnishings, Fur_length, Cystathionine) %>%
  pivot_longer(c(-dog_id, -Cystathionine))

# Clean trait names for display
df <- df %>%
  mutate(name = str_replace_all(name, "_", " "))

# Common theme for the cystathionine-by-trait boxplots
theme_boxes <- function() { 
  theme_cowplot(8) %+replace%
    theme(
      plot.title   = element_text(size = 7, face = "bold", hjust = 0.5,
                                  margin = margin(b = 2)),
      axis.title.y = element_text(angle = 90, size = 7,,face="bold"),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 6, vjust = 1,
                                  margin = margin(b = 0), hjust = 0.5),
      axis.text.y  = element_text(size = 6,
                                  margin = margin(2, 0, 0, 4), hjust = 1),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey80"),
      axis.ticks.x       = element_blank(),
      legend.position    = "none",
      axis.line.x        = element_blank(),
      axis.line.y        = element_line(linewidth = 0.25),
      axis.ticks.y       = element_line(linewidth = 0.25),
      strip.background   = element_rect(fill = "white", color = NA),
      strip.text         = element_text(
        color  = "black", hjust = 0.5, vjust = 0,
        face   = "bold", size  = 6, margin = margin(b = 5)
      ),
      plot.margin = margin(t = 1, r = 5, l = 5, b = 0)
    )
}

# =====================================
# Recode trait values to human-readable labels and compute annotations
# =====================================
pd <- df %>%
  mutate(name = tolower(name)) %>%
  mutate(name=if_else(name == "furnishings",paste0(name, "\n(long, coarse\nfacial hair)"),name))
  

pd <- pd %>%
  mutate(
    # Convert numeric codes into display labels for each facet
    label = case_when(
      name == "fur length" & value == 1 ~ "short",
      name == "fur length" & value == 2 ~ "med.",
      name == "fur length" & value == 3 ~ "long",
      name != "fur length" & value == 1 ~ "no",
      name != "fur length" & value == 2 ~ "yes",
      TRUE ~ NA_character_
    )
  )

# Mean cystathionine per facet/category (annotated at bottom of each panel)
means <- pd %>%
  group_by(name, label) %>%
  summarise(
    mean = round(mean(Cystathionine, na.rm = TRUE), 3),
    .groups = "drop"
  )

# Per-facet association tests and p-value labels
results <- pd %>%
  group_by(name) %>%
  summarise(
    # x position for the p-value annotation within each facet
    xpos = median(sort(unique(value))),
    test_result = list({
      # Use the current facet's data for testing
      d <- pick(everything())
      if (name[1] == "fur length") {
        # Ordinal association for fur length
        x <- d$Cystathionine
        y <- d$value
        if (is.factor(y)) y <- as.integer(y)
        stats::cor.test(x, y, method = "kendall")
      } else {
        # Two-group comparison for furnishings
        stats::wilcox.test(Cystathionine ~ label, data = d, exact = FALSE)
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(
    pval      = map_dbl(test_result, ~ .x$p.value),
    pval_text = paste0("p = ", signif(pval, 2))
  )

# Plot limits (used to place annotations consistently across facets)
min_y <- min(df$Cystathionine)
max_y <- max(df$Cystathionine)

# Control facet ordering and label ordering within facets
facet_levels <- c(
  "furnishings\n(long, coarse\nfacial hair)",
  "fur length"
)

pd <- pd %>%
  mutate(
    name  = factor(name, levels = facet_levels),
    label = factor(label, levels = c("no","yes","short","med.","long"))
  )

# Ensure annotation data frames use the same facet factor levels as the plotting data
results <- results %>% mutate(name = factor(name, levels = facet_levels))
means   <- means   %>% mutate(name = factor(name, levels = facet_levels))

# =====================================
# Panel C: cystathionine levels vs fur traits (with p-values and mean annotations)
# =====================================
p <- ggplot(pd, aes(x = label, y = Cystathionine))
p <- p + geom_boxplot(outlier.size = 0.25, fill = "white", width = 0.5, linewidth = 0.25, 
                      outlier.shape = NA)
p <- p + geom_jitter(width = 0.15, height = 0, shape = 16, size = 0.5, alpha = 0.5)

# Add per-facet p-value at the top
p <- p + geom_text(  aes(label = pval_text,x=xpos),  size = 1.75,  y = max(pd$Cystathionine, na.rm = TRUE) + 0.1,  data = results )

# Add per-category mean at the bottom
p <- p + geom_text(  aes(label = mean),  color = "#a84a3c",  
                     size = 1.75,  y   = min(pd$Cystathionine, na.rm = TRUE) - 0.12,  
                     data = means)

p <- p + facet_grid(.~name,scales="free_x",space="free_x")
p <- p + scale_y_continuous(  "cystathionine levels",  limits = c(min_y - 0.08, max_y + 0.05),  breaks = c(-1, 0, 1) )
p <- p + theme_boxes() + theme(plot.margin = margin(b=15,t=5))

p_hair <- p
p

# =====================================
# Load or build the LMER-based panels used elsewhere in the figure
# =====================================
rdata_file <- "fig_LMER.for_fig_4_FUR.RData"
build_script <- "plot_fig_4_LMER_panel.R"

if (file.exists(rdata_file)) {
  load(rdata_file)
} else {
  source(build_script)
}

# =====================================
# Panels D: breed-frequency plots for lead SNPs at chr13 (furnishings) and chr32 (length)
# =====================================
p13_furn <- make_maf_boxplot_snps(in_snps = c("chr13:8839016:G:A","chr13:8865645:G:A","chr13:8867526:C:T"), in_trait = "furnishings")
p32_length <- make_maf_boxplot_snps(in_snps = c("chr32:35494497:C:A","chr32:35494936:C:T","chr32:35496522:G:A"),"length")

# =====================================
# Assemble final multi-panel layout and save
# =====================================
row_gwas <- plot_grid(p_fur_gwas,p_hair,nrow=1,label_size=8,labels=c("","","C"),rel_widths = c(2,0.9))
row_lmer <- plot_grid(grid_lmer_size,grid_lmer_fur,nrow=1,labels=c("E","F"),label_size=8,rel_widths = c(0.8,1))
row_freq <- plot_grid(p13_furn,p32_length,ncol=1,labels=c("D"),label_size=8)

row_freq_lmer <- plot_grid(row_freq,row_lmer,nrow=1,rel_widths = c(0.55,2))

grid <- plot_grid(row_gwas,row_freq_lmer,ncol=1,rel_heights = c(1,1))

ggsave(plot=grid,filename=outpdf,width=6.5,height=5)
