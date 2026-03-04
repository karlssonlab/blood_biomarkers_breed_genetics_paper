## =====================================
## Libraries and global setup
## =====================================
library(tidyverse)  # data manipulation + ggplot2
library(ggpubr)     # stat_compare_means()
library(cowplot)    # plot_grid() and theme_cowplot()
library(scales)     # comma() label formatting

# Project-specific helper functions (e.g., breed name cleanup, yes/no conversion)
source("util.cleaning.R")

# Output figure path
outpdf <- "../../../figures_pdf/fig_1_DEMOG.pdf"

# Maximum weight shown on weight-based panels (kg)
max_weight_plotted <- 70

## =====================================
## Load data and harmonize key fields
## =====================================
# Load pre-saved DAP supplement data if not already present in the workspace
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_S15_BREED_ANC, etc.
}

## Standardize breed naming across input tables to a consistent format
raw_S15_BREED_ANC        <- standardize_breed_names(raw_S15_BREED_ANC, pop)
raw_S12_BREED_PHENOS <- standardize_breed_names(raw_S12_BREED_PHENOS, breed)
raw_S14_SEQ_DOGS <- standardize_breed_names(raw_S14_SEQ_DOGS, standardized_breed)

# Start from sequencing cohort dogs table
indata <- raw_S14_SEQ_DOGS

## Convert yes/no columns to TRUE/FALSE (logical) for consistent filtering
indata <- convert_yesno_to_logical(raw_S14_SEQ_DOGS)

# Ancestry table: normalize breed labels for joining/plotting
anc <- raw_S15_BREED_ANC %>% rename(Breed = pop) %>%
  mutate(Breed = tolower(str_replace_all(Breed, "_", " ")))

# Breed size grouping table (one row per breed)
breedsizes <- raw_S12_BREED_PHENOS %>%
  rename(Breed = breed) %>%
  mutate(Breed = tolower(Breed)) %>%
  select(Breed, Size.Group) %>%
  distinct()

# Ensure every ancestry breed has a size group; default missing to "Variable"
breedsizes <- anc %>% select(Breed) %>% distinct() %>%
  full_join(breedsizes) %>% replace_na(list(Size.Group = "Variable")) %>%
  select(Breed, Size.Group)

# Keep analyzable dogs (must have ID and weight) and harmonize age column name
indata <- indata %>%
  filter(!is.na(dog_id) & !is.na(weight_kg)) %>%
  mutate(source = "seqDAP") %>%
  rename(age_years = Estimated_Age_Years_at_HLES)

## =====================================
## Define cohorts/sets used in Age/Weight comparisons
## =====================================
# Dogs with survey data
labels <- indata %>%
  filter(survey_data) %>%
  select(dog_id) %>%
  mutate(set = "Surveys")

# Dogs with any blood trait data (CBC/chemistry/metabolites)
labels <- labels %>%
  bind_rows(
    indata %>%
      filter(complete_blood_counts | serum_chemistry_panel | plasma_metabolites) %>%
      select(dog_id) %>%
      mutate(set = "Blood traits")
  )

## Reduce to plotting variables and define population-type categories
indata <- indata %>%
  select(
    dog_id, source, age_years, Breed=standardized_breed, single_breed, sex,
    sterilization_status, weight_kg, top_breed, max_percent_ancestry,
    num_breeds_greater_than_five_percent
  ) %>%
  mutate(
    breed_type = if_else(single_breed, "single-breed", "mixed breed"),
    poptype    = if_else(
      breed_type == "single-breed", "single-breed",
      if_else(!is.na(max_percent_ancestry) & max_percent_ancestry >= 0.4,
              "simple mix", "complex mix")
    )
  )

# Histogram bin width used later for weight binning
binwidth <- 3
indata <- indata %>%
  mutate(weight_rounded = floor(weight_kg / binwidth) * binwidth)

## =====================================
## Plot themes (small-format figure panels)
## =====================================
theme_categ_y <- function() {
  theme_cowplot(12) %+replace%
    theme(
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      legend.position = c(0.98, 0.98),
      legend.text     = element_text(size = 5.5),
      legend.title    = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.justification.inside = c(1.05, 1.25),
      legend.box.just = "right",
      axis.ticks.y    = element_blank(),
      axis.line.y     = element_blank(),
      axis.line.x     = element_line(linewidth = 0.25),
      axis.ticks.x    = element_line(linewidth = 0.25),
      strip.text.x    = element_text(size = 7, hjust = 0.5,
                                     face = "bold", lineheight = 1),
      plot.margin     = margin(5, 3, 5, 3)
    )
}

theme_quant <- function() {
  theme_cowplot(8) %+replace%
    theme(
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6),
      legend.position     = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.text         = element_text(size = 5),
      legend.title        = element_blank(),
      legend.key.size     = unit(0.2, "cm"),
      legend.box.just     = "right",
      axis.line           = element_line(linewidth = 0.25),
      axis.ticks          = element_line(linewidth = 0.25),
      panel.grid.major.y  = element_line(linewidth = 0.2, color = "grey80"),
      strip.text.x        = element_text(size = 7, hjust = 0.5,
                                         face = "bold", lineheight = 1),
      plot.margin         = margin(5, 3, 5, 3)
    )
}

## Quick check: weight differences between single-breed vs mixed-breed (not plotted)
wilcox.test(weight_kg ~ single_breed, data = indata)

## =====================================
## Prepare age/weight comparison data and y-axis labels
## =====================================
labels <- labels %>%
  as_tibble() %>%                          # ensure it's a tibble/data.frame
  mutate(set = as.character(set))         # force set to atomic character

# Counts per set used for y-axis labels with N
counts <- labels %>%
  group_by(set) %>%
  count() %>%
  ungroup()

# Helper: upper whisker used to set reasonable x-axis limits for boxplots
get_upper_whisker <- function(x) {
  x <- x[!is.na(x)]
  qs <- quantile(x, probs = c(0.25, 0.75))
  iqr <- qs[2] - qs[1]
  max(x[x <= qs[2] + 1.5 * iqr])
}

# Join set labels onto dog-level data and split into age/weight frames
pd <- indata %>% inner_join(labels, by = "dog_id")
pd <- pd %>% select(dog_id, set, weight_kg, age_years)
pdAge    <- pd %>% select(dog_id, set, age_years)   %>% filter(!is.na(age_years))
pdWeight <- pd %>% select(dog_id, set, weight_kg)   %>% filter(!is.na(weight_kg))

# Build y-axis labels that include sample size per set
counts <- counts %>%
  mutate(ylabel = paste(set, "\n(N=", comma(n, accuracy = 1),")", sep = ""))
ylabels <- counts %>% arrange(set) %>% pull(ylabel)
ybreaks <- counts %>% arrange(set) %>% pull(set)

# Ensure consistent ordering across plots
pdAge$set    <- factor(pdAge$set,    levels = ybreaks)
pdWeight$set <- factor(pdWeight$set, levels = ybreaks)

# Pairwise comparison: Surveys vs Blood traits
my_comparisons <- list(c(ybreaks[1], ybreaks[2]))

## =====================================
## Panel A: Age distributions by set
## =====================================
maxX <- ceiling(get_upper_whisker(pdAge$age_years)) * 1.1

p <- ggplot(pdAge, aes(x = age_years, y = set))
p <- p + geom_boxplot(fill = NA, outlier.size = 0.5, outlier.alpha = 0.25,
                      outlier.shape = NA, linewidth = 0.25, width = 0.5)
p <- p + stat_compare_means(comparisons = my_comparisons, size = 1.75, label.y = maxX * 0.9,
                            tip.length = 0.01)
p <- p + scale_y_discrete("", breaks = ybreaks, labels = ylabels)
p <- p + scale_x_continuous("age (years)", limits=c(0,maxX), expand = expansion(mult = c(0.01, 0.1)))
p <- p + theme_categ_y()
pbox_Age <- p

# Summary stats (console output)
pdAge %>%
  group_by(set) %>%
  summarize(mean = mean(age_years), sd = sd(age_years), n = n())

## =====================================
## Panel B: Weight distributions by set
## =====================================
maxX <- ceiling(get_upper_whisker(pdWeight$weight_kg)) * 1.1

p <- ggplot(pdWeight, aes(x = weight_kg, y = set))
p <- p + geom_boxplot(fill = NA, outlier.size = 0.5, outlier.alpha = 0.25,
                      outlier.shape = NA, linewidth = 0.25, width = 0.5)
p <- p + stat_compare_means(comparisons = my_comparisons, size = 1.75, label.y = maxX * 0.9,
                            tip.length = 0.01)
p <- p + scale_alpha_manual(values = c(0.9, 0.2))
p <- p + scale_y_discrete("", breaks = ybreaks, labels = ylabels)
p <- p + scale_x_continuous("dog's weight (kg)", limits=c(0,maxX), expand = expansion(mult = c(0.01, 0.1)))
p <- p + theme_categ_y()
pbox_Weight <- p

## =====================================
## Panel C: Cumulative fraction by number of contributing breeds (>5% ancestry)
## =====================================
pd <- indata %>%
  select(
    dog_id, source, single_breed, Breed, top_breed,
    max_percent_ancestry, poptype, num_breeds_greater_than_five_percent
  ) %>%
  mutate(num_breeds_greater_than_five_percent =
           if_else(single_breed, 1L, num_breeds_greater_than_five_percent))

# Total dogs per source (used downstream for joins/compatibility with prior code)
pd_tot <- pd %>%
  select(source, dog_id) %>%
  distinct() %>%
  group_by(source) %>%
  count() %>%
  rename(number_total = n)

# Total dogs used as denominator for fractions (note: uses all dogs, not per-source)
number_total <- length(pd$dog_id)

# Count dogs by poptype and number-of-breeds, compute cumulative fraction
pd <- pd %>%
  filter(!is.na(num_breeds_greater_than_five_percent)) %>%
  group_by(poptype, source, num_breeds_greater_than_five_percent) %>%
  count() %>%
  arrange(source, poptype, -num_breeds_greater_than_five_percent) %>%
  ungroup() %>%
  group_by(source, poptype) %>%
  mutate(cum = cumsum(n)) %>%
  filter(num_breeds_greater_than_five_percent <= 6 &
           num_breeds_greater_than_five_percent > 0) %>%
  left_join(pd_tot, by = "source") %>%
  mutate(frac = cum / number_total)

pd$poptype <- factor(pd$poptype, levels = c("single-breed", "simple mix", "complex mix"))

# Totals used for text labels above bars (for bins > 1)
pdsum <- pd %>%
  group_by(source, num_breeds_greater_than_five_percent) %>%
  summarize(total = sum(cum), .groups = "drop") %>%
  filter(num_breeds_greater_than_five_percent > 1)

p <- ggplot(pd, aes(x = num_breeds_greater_than_five_percent)) +
  geom_bar(aes(fill = poptype, y = frac), stat = "identity", width = 0.75,
           color = "#252525", linewidth = 0.25)
p <- p + geom_text(aes(label = total, y = total / number_total),
                   vjust = 0, hjust = 0.5, nudge_y = 0.02, size = 1.75,
                   data = pdsum)
p <- p + scale_fill_manual(values = c(
  "single-breed" = "#fafafa",
  "simple mix"   = "#969696",
  "complex mix"  = "#252525"
))
p <- p + scale_x_continuous("# breeds (>5%)",
                            breaks = c(1:6),
                            labels = paste(c(1:6), "+", sep = ""))
p <- p + scale_y_continuous(
  "% all dogs", breaks = c(0:4 / 4), limits = c(0, 1),
  labels = paste(c(0:4 / 4) * 100, "%", sep = ""),
  expand  = expansion(mult = c(0, 0))
)
p <- p + theme_quant() +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_text(margin = margin(t = 2))
  )
psum <- p

## =====================================
## Prepare dog-level data for weight-binned breed-size plots (age > 2)
## =====================================
pd <- indata %>%
  select(
    dog_id, age_years, Breed, top_breed, breed_type,
    weight_rounded, single_breed, max_percent_ancestry
  ) %>%
  mutate(
    Breed = if_else(
      single_breed,
      if_else(breed_type == "mixed breed", top_breed, Breed),
      "mixed"
    )
  ) %>%
  left_join(breedsizes %>% select(Breed, Size.Group) %>% distinct(), by = "Breed") %>%
  filter(age_years > 2)

# Recreate poptype within this filtered dataset
pd <- pd %>%
  mutate(poptype = if_else(
    single_breed, "single-breed",
    if_else(max_percent_ancestry >= 0.4, "simple mix", "complex mix")
  ))

pd$poptype <- factor(pd$poptype, levels = c("single-breed", "simple mix", "complex mix"))

## =====================================
## Panel D: Weight density by ancestry complexity + Hartigan’s dip test annotation
## =====================================
# Build density dataset including an "all mixed" aggregation (duplicate mixed-breed rows)
pd_density <- indata %>% select(dog_id, poptype, weight_kg) %>% drop_na()
pd_density <- pd_density %>% filter(poptype!="single-breed") %>% mutate(poptype="all mixed") %>% bind_rows(pd_density)
pd_density <- pd_density %>% filter(poptype != "simple mix") %>% mutate(poptype = factor(poptype, levels = c("single-breed", "all mixed", "complex mix")))

# Format dip-test results into per-facet text annotations
make_bimodality_annotations <- function(results) {

  results %>%
    mutate(
      # Format dip p-value
      dip_label = paste0(
        "Dip p = ",
        ifelse(is.na(dip_p), "NA", ifelse(dip_p==0,"<1e5",signif(dip_p, 2)))
      ),
      dip_statistic_str = paste0(
        "Dip = ",
        ifelse(is.na(dip_statistic), "NA", round(dip_statistic, 4))
      ),
      n_label = paste0("N = ", comma(n,accuracy=1)),
      # Combine into final annotation
      annotation = paste(n_label, dip_statistic_str, dip_label, sep = "\n")
    ) %>%
    select(poptype, annotation)
}

## Grouped Hartigan’s dip test (per poptype)
results <- pd_density %>%
  filter(!is.na(weight_kg)) %>%
  group_by(poptype) %>%
  group_modify(~ {
    # Force to plain numeric vector (defensive)
    x <- .x %>%
      dplyr::pull(weight_kg) %>%
      unlist(use.names = FALSE) %>%
      as.numeric()

    n <- length(x)

    if (n < 3 || length(unique(x)) < 2) {
      return(tibble(
        n             = n,
        dip_statistic = NA_real_,
        dip_p         = NA_real_,
        dip_reject    = NA
      ))
    }

    # Hartigan’s dip test with bootstrap p-value
    dip_res <- diptest::dip.test(x,B = 100000)

    tibble(
      n             = n,
      dip_statistic = unname(dip_res$statistic),
      dip_p         = dip_res$p.value,
      dip_reject    = dip_res$p.value < 0.05
    )
  }) %>%
  ungroup()

annot <- make_bimodality_annotations(results)

# Annotation placement within facets
dens_max <- 0.04
annot <- annot %>% mutate(weight_kg=3, y=dens_max+0.002)

# Facet labels
label_map <- c(
  "single-breed" = "Single-\nbreed",
  "simple mix"   = "Simple\nmixes",
  "complex mix"  = "Complex\nmixes",
  "all mixed"  = "Mixed\nbreed"
)

p <- ggplot(pd_density, aes(x = weight_kg))
p <- p + geom_density(aes(fill = poptype), color = "black", linewidth = 0.25, alpha=0.5)
p <- p + geom_text(aes(y=y, label=annotation), data=annot, hjust=0, vjust=0, size=1.25, lineheight=0.9)
p <- p + scale_fill_manual(values = c("single-breed"="#fafafa", "complex mix"="#252525","all mixed"= "#969696"))
p <- p + scale_y_continuous("density", limits=c(0,dens_max+0.007), breaks=c(0,0.01,0.02,0.03,0.04), expand = expansion(0,0))
p <- p + scale_x_continuous("dog's weight (kg)", expand = expansion(mult = c(0.01,0)))
p <- p + coord_cartesian(xlim = c(0, max_weight_plotted))
p <- p + facet_wrap(~poptype, nrow=1, labeller = labeller(poptype = label_map))
p <- p +  theme_quant()
p <- p + theme(legend.position="none",
               strip.text.x = element_text(size = 5, hjust = 0, vjust = 0, face = "bold", margin = margin(t = 2, b = 2), color = "black"),
               strip.background = element_rect(fill = "white", color = NA))
p_density <- p

## =====================================
## Panel E: Single-breed weight distribution by breed size group
## =====================================
pdCounts <- pd %>%
  filter(single_breed & !is.na(Size.Group)) %>%
  group_by(poptype, weight_rounded, Size.Group) %>%
  count() %>%
  rename(ndogs = n)

pdCounts$Size.Group <- factor(pdCounts$Size.Group,
                              levels = c("Toy and Small", "Medium", "Standard",
                                         "Large", "Giant", "Variable"))

p <- ggplot(
  pdCounts %>% filter(weight_rounded <= max_weight_plotted),
  aes(x = weight_rounded, y = ndogs, fill = Size.Group)
)
p <- p + geom_bar(
  stat = "identity",
  width = binwidth,
  alpha = 1,
  color = "#252525",
  linewidth = 0.25
)
p <- p + scale_fill_manual(
  values = c(
    "Toy and Small" = "#8dd3c7",
    "Medium"        = "#ffffb3",
    "Standard"      = "#bebada",
    "Large"         = "#fb8072",
    "Giant"         = "#80b1d3",
    "Variable"      = "#bdbdbd"
  )
)
p <- p + scale_y_continuous(
  "# dogs (single-breed)",
  breaks = c(-5:6) * 100,
  labels = abs(c(-5:6) * 100)
)
p <- p + scale_x_continuous(
  "dog's weight (kg)",
  limits = c(-binwidth / 2, max_weight_plotted + binwidth / 2),
  expand = expansion(mult = c(0.025, 0))
)
p <- p + labs(fill = "Breed size")
p <- p + theme_quant() +
  theme(legend.title = element_text(size = 5, face = "bold"))
pbreed <- p

## =====================================
## Panel F: Complex-mix ancestry fraction by breed size group across weight bins
## =====================================
pdSums <- indata %>%
  filter(!single_breed & max_percent_ancestry <= 0.4) %>%
  select(dog_id, weight_rounded) %>%
  inner_join(anc, by = "dog_id")

# Sum ancestry percentages per size group and weight bin (fractional dogs by ancestry)
pdSums <- pdSums %>%
  left_join(breedsizes, by = "Breed") %>%
  group_by(Size.Group, weight_rounded) %>%
  summarize(ndogs = sum(pct), .groups = "drop")

# Convert to fraction of total ancestry mass across all bins
pdSums <- pdSums %>%
  mutate(frac = ndogs / sum(pdSums$ndogs))

pdSums$Size.Group <- factor(
  pdSums$Size.Group,
  levels = c("Toy and Small", "Medium", "Standard", "Large", "Giant", "Variable")
)

p <- ggplot(
  pdSums %>% filter(weight_rounded <= max_weight_plotted),
  aes(x = weight_rounded, y = frac)
)
p <- p + geom_bar(
  aes(fill = Size.Group),
  stat = "identity",
  width = binwidth,
  alpha = 1,
  color = "#252525",
  linewidth = 0.25
)
p <- p + scale_fill_manual(
  values = c(
    "Toy and Small" = "#8dd3c7",
    "Medium"        = "#ffffb3",
    "Standard"      = "#bebada",
    "Large"         = "#fb8072",
    "Giant"         = "#80b1d3",
    "Variable"      = "#bdbdbd"
  )
)
p <- p + scale_y_continuous(
  "frac. ancestry (complex mix)",
  breaks = c(0, 0.05, 0.1)
)
p <- p + scale_x_continuous(
  "dog's weight (kg)",
  limits = c(-binwidth / 2, max_weight_plotted + binwidth / 2),
  expand = expansion(mult = c(0.025, 0))
)
p <- p + theme_quant() + theme(legend.position = "none")
pmixed <- p

## =====================================
## Combine panels and save figure
## =====================================
panel1 <- plot_grid(pbox_Age, pbox_Weight, ncol = 1,
                    labels = LETTERS[1:2], label_size = 8)

row1 <- plot_grid(panel1, psum, p_density, labels = c("", "C", "D"), nrow = 1, label_size = 8, rel_widths = c(0.8, 0.8, 1) )
row2 <- plot_grid(pbreed, pmixed, labels = c("E", "F"), nrow = 1, label_size = 8, rel_widths = c(0.9, 0.9, 1) )
grid <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

ggsave(plot = grid, filename = outpdf, width = 5, height = 4)
