library(tidyverse)
library(cowplot)
source("util.cleaning.R")

# Output figure path
outpdf <- "../../../figures_pdf/fig_S3_SURVEYS.pdf"

## =====================================
## Load all data (RDS contains a list; expose elements into .GlobalEnv)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_pheno, raw_anc, raw_S1_SURVEYS, etc.
}

## =====================================
## Plot styling helper (compact cowplot theme for small multi-panel figure)
## =====================================
theme_bar <- function(){
  theme_cowplot(8) %+replace%
    theme(
      plot.margin = margin(t=5,r=5,l=10,b = 5),
      plot.title = element_text(size = 5,face="bold",hjust=0),
      axis.ticks = element_line(linewidth=0.25),
      axis.line = element_line(linewidth = 0.25),    
      axis.text = element_text(size = 7),
      axis.title = element_text(size =8 ),
      
      legend.position = "none",
      panel.grid.major=element_line(color="grey80",linewidth=0.2),
      
      
    )
}

## =====================================
## Prepare survey data
## =====================================
# Ensure 'included' is logical so ggplot groups correctly
survey_data <- raw_S1_SURVEYS %>%
  mutate(included = as.logical(included))

# Build x-axis labels with group sample sizes
xlabels <- survey_data %>% group_by(included) %>% count() %>% 
  mutate(label=if_else(included,"Included","Excluded")) %>%
  mutate(label=paste0(label,"\nN=",n))
  

## =====================================
## Panel A: Missingness by inclusion status (boxplot)
## =====================================
p <- ggplot(survey_data, aes(x = included, y = proportion_NA, fill = included)) +
  geom_boxplot(width = 0.4, outlier.size = 1, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#a50f15", "FALSE" = "#2171b5")) +
  labs(
    y = "Proportion Missing Responses"
  ) +
  scale_x_discrete("",labels=xlabels$label,breaks=xlabels$included) + 
  theme_bar()
p_missing <- p

## =====================================
## Summarize reasons for exclusion (collapsed into broader categories)
## =====================================
excluded <- survey_data %>% 
  filter(included == "FALSE") %>%
  replace_na(list(exclusion_reason="other")) 

# Collapse several detailed/rare reasons into "other"
excluded <- excluded %>% mutate(exclusion_reason=if_else(exclusion_reason %in% c("excluded any pre-defined score in survey",
                                                                            "date",
                                                                            "time in minutes by itself not helpful",
                                                                            "unbalanced phenotype distribution (kurtosis > 7)"),"other",exclusion_reason))

# Group owner-related exclusions
excluded <- excluded %>% mutate(exclusion_reason=if_else(exclusion_reason %in% c("owner contact section excluded",
                                                                            "not baseline related",
                                                                            "owner demographic section excluded"),
                                                                            "owner information",exclusion_reason))

# Group non-standard response format exclusions
excluded <- excluded %>% mutate(exclusion_reason=if_else(exclusion_reason %in% c("free-text",
                                                                                 "nested",
                                                                                 "date",
                                                                                 "health status question too broad or too specific"),
                                                         "non-standard response types",exclusion_reason))

# Clean up wording for readability in the plot
excluded <- excluded %>% mutate(exclusion_reason=str_replace(exclusion_reason,"dog demographic section excluded except weight","dog demography"))
excluded <- excluded %>% mutate(exclusion_reason=str_remove(exclusion_reason,"phenotype "))
excluded <- excluded %>% mutate(exclusion_reason=str_remove(exclusion_reason," excluded"))

# Count excluded questions per reason (sorted for plotting)
excluded <- excluded %>% 
  count(exclusion_reason, sort = TRUE)


## =====================================
## Panel B: Exclusion reasons (horizontal bar chart)
## =====================================
p <- ggplot(excluded, aes(x = reorder(exclusion_reason, n), y = n)) +
  geom_col(fill = "#2171b5",color="grey30",linewidth=0.5,width=0.5) +
  coord_flip() +
  labs(
    x = "Exclusion Reason",
    y = "Number of Excluded Questions"
  ) +
  theme_bar()

p_reason <- p

## =====================================
## Combine panels and save
## =====================================
grid <- plot_grid(p_missing,p_reason,nrow=1,label_size = 8,labels=LETTERS,rel_widths = c(0.5,1))
ggsave(grid,filename=outpdf,width=5,height=3)
grid
