## =====================================
## Libraries and setup
## =====================================
# Core data manipulation + ggplot2 + purrr helpers used throughout
library(tidyverse)

# Read Google Sheets data (read_sheet)
library(googlesheets4)

# Plot styling and layout utilities (theme_cowplot, plot_grid)
library(cowplot)

# Percent label formatting for y-axis
library(scales)

# Project-specific helper functions (not necessarily used directly in this script)
source("util.cleaning.R")

# Input Google Sheet and output figure path
supp_sheet_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/edit?usp=sharing"
outpdf <- "../../../figures_pdf/fig_S13_NUM_BREEDS.pdf"


## =====================================
## Load all data
## =====================================
# Load breed ancestry table from Google Sheets once per session; sanitize column names
if (!exists("raw_S15_BREED_ANC")  ) {
   raw_S15_BREED_ANC <- read_sheet(supp_sheet_url, range = "Data S15_BREED_ANC",
                        na = c("", "NA", "#N/A"), skip = 0) %>% as_tibble() %>%
    rename_with(~ gsub("[?'~ /]", ".", .x))

}

# Standardize column naming and drop missing ancestry fractions
anc_all <- raw_S15_BREED_ANC %>% rename(breed=pop) %>% filter(!is.na(pct))

# Total ancestry contribution per breed across all dogs (used to define "top X" breeds)
breed_totals <- anc_all %>%
  group_by(breed) %>%
  summarise(total_frac = sum(pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_frac))



## ----------------------------
## Coverage function: for top X breeds, how many dogs have all breeds >= cutpct covered?
## ----------------------------
# For a given ancestry cutoff (cutpct), check each dog: are ALL of its >=cutpct breeds within the top X breeds?
count_covered_dogs <- function(X, cutpct) {
  
  # Keep only dog/breed entries meeting the ancestry fraction cutoff
  anc_cut <- anc_all %>%
    filter(pct >= cutpct) %>%
    select(dog_id, breed)
  
  # Define the set of top X breeds by total contribution
  top_breeds <- breed_totals %>%
    slice_head(n = X) %>%
    transmute(breed, fnd = TRUE)
  
  # Mark which cutoff-passing breeds are found in the top-X set (missing => FALSE)
  fnd <- anc_cut %>%
    left_join(top_breeds, by = "breed") %>%
    mutate(fnd = replace_na(fnd, FALSE))
  
  # A dog is "covered" if all of its cutoff-passing breeds are present in the top-X list
  fnd %>%
    group_by(dog_id) %>%
    summarise(all_fnd = all(fnd), .groups = "drop") %>%
    summarise(n_dogs_covered = sum(all_fnd), .groups = "drop") %>%
    pull(n_dogs_covered)
}

## ----------------------------
## Build coverage curves for cutpct = 0.05 and 0.02
## ----------------------------
# X ranges from 1..(#breeds), representing the "top X breeds" included in the model
max_X <- nrow(breed_totals)
X_vals <- unique(c(seq(1, max_X, by = 1), max_X))

# Evaluate two ancestry inclusion thresholds
cutpcts <- c(0.02, 0.05)

# Denominator: number of dogs that have at least one breed at/above the cutoff
n_total_by_cutpct <- tibble(cutpct = cutpcts) %>%
  mutate(
    n_total_dogs = map_int(
      cutpct,
      ~ n_distinct(anc_all %>% filter(pct >= .x) %>% pull(dog_id))
    )
  )

# For each (X, cutoff), compute covered dogs and convert to fraction of eligible dogs
coverage_curve <- crossing(
  X = X_vals,
  cutpct = cutpcts
) %>%
  mutate(
    n_dogs_covered = map2_int(X, cutpct, count_covered_dogs)
  ) %>%
  left_join(n_total_by_cutpct, by = "cutpct") %>%
  mutate(
    frac_dogs_covered = n_dogs_covered / n_total_dogs
  )


## ----------------------------
## Plot coverage curves and save to PDF
## ----------------------------
p <- ggplot(coverage_curve, aes(x = X, y = frac_dogs_covered,group=factor(cutpct))) 
p <- p + geom_line(aes(color=factor(cutpct)),linewidth = 0.5) 
p <- p + scale_x_continuous("Number of breeds (top X by total contribution)",breaks=c(0:5)*25)
p <- p + scale_y_continuous(labels = scales::percent,name = "% dogs with all breed ancestry modelled")
p <- p + scale_color_manual(name   = "Breed\nancestry\nincluded:",
                            values = c("0.05" = "#0072B2", "0.02" = "#D55E00"),
                            labels = c("0.05" = ">5%", "0.02" = ">2%"))
p <- p + theme_cowplot(8) +  theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
                                   legend.position = c(0.98,0.02),
                                   legend.text = element_text(size=8),legend.title = element_text(size=8,face="bold"),
                                   legend.key.size = unit(0.2, 'cm'),
                                   legend.justification.inside = c(1, 0),
                                   legend.box.just = "right",
                                   legend.background = element_rect(color = "grey20",fill = "white",linewidth = 0.25),
                                   legend.margin = margin(2,2,2,2),
                                   axis.text = element_text(margin=margin(t=5,b=5)),
                                   axis.title = element_text(face="bold")
)
p_cov <- p
print(p_cov)
grid <- plot_grid(p_cov, ncol = 1)
ggsave(plot = grid, filename = outpdf, width = 3.5, height = 3)

## ----------------------------
## Quick lookup: X needed for a target dog-coverage fraction (example: 50%, 80%, 90%)
## ----------------------------
# For each cutoff and target coverage, find the smallest X achieving the target (or NA if unreachable)
targets <- c(0.5, 0.8, 0.9,0.95)
needed_for_targets <- crossing(
  cutpct = sort(unique(coverage_curve$cutpct)),
  target = targets
) %>%
  group_by(cutpct, target) %>%
  group_modify(~{
    df <- coverage_curve %>% filter(cutpct == .y$cutpct[1])
    
    out <- df %>%
      filter(frac_dogs_covered >= .y$target[1]) %>%
      slice_min(X, with_ties = FALSE)
    
    if (nrow(out) == 0) {
      tibble(
        X_needed = NA_integer_,
        frac_dogs_covered = max(df$frac_dogs_covered, na.rm = TRUE)
      )
    } else {
      tibble(
        X_needed = out$X,
        frac_dogs_covered = out$frac_dogs_covered
      )
    }
  }) %>%
  ungroup()
needed_for_targets %>% arrange(-target,cutpct) %>% mutate(target=paste0(target*100,"%"))
