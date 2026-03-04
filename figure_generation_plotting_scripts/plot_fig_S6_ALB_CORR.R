library(tidyverse)
library(cowplot)

# Load helper functions (e.g., label formatting utilities used below)
source("util.cleaning.R")

# Output path for the final figure
outpdf <- "../../../figures_pdf/fig_S6_ALB_CORR.pdf"

# =====================================
# Load data (once) and expose tables to the global environment
# =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects like raw_S7_ALB, raw_S2_ALL_PHENOS, etc.
}

# =====================================
# Reshape albumin correlation table into long format for plotting
# =====================================
alb <- raw_S7_ALB %>% rename(blood_trait = ...1) %>% pivot_longer(-blood_trait) %>% rename(label1=blood_trait,label2=name,correlation=value)

# Quick check: traits in ALB table that are not present in the phenotype label table
alb %>% filter(!(label1 %in% raw_S2_ALL_PHENOS$plot_label))

# =====================================
# Define the variable order as displayed in the figure
# =====================================
order <- tibble(phenotype=c(
  "krt_cp_albumin_value_ln_transformed","krt_cp_globulins_value_ln_transformed",
  "krt_cp_alb_glob_ratio_value_ln_transformed", "Tryptophan",
  "L-Kynurenine", "Indole-3-Lactate", "Indole-3-Acetic-Acid", "Indole-3-Propionate",
  "Margaric-Acid"),order=rev(c(1:9))) %>% left_join(raw_S2_ALL_PHENOS %>% select(phenotype,plot_label))

# Create all pairwise combinations of variables and carry both axis orders
ordering <- order %>% rename(label1=plot_label,order1=order) %>% select(-phenotype)
ordering <- crossing(ordering,order %>% rename(label2=plot_label,order2=order) )

# =====================================
# Prepare plotting data: pretty labels and ordering for both axes
# =====================================
pd <- alb
names <- tibble(name=c(pd$label1,pd$label2)) %>% distinct()

# Insert line breaks for long labels (including handling of "/" labels)
names <- names  %>%  mutate(label_newline = map_chr(name, replace_middle_break)) %>%
 mutate(label_newline=str_replace(label_newline,"/","/\n"))
names <- names

# Join pretty-printed labels for both dimensions
pd <- pd %>% left_join(names %>% select(label1=name,label1_pr=label_newline)) %>%
  left_join(names %>% select(label2=name,label2_pr=label_newline))

# Keep only pairs in the predefined ordering and attach axis ranks
pd <- pd %>% inner_join(ordering)

# Recompute the final ordering based on pretty labels
rm(ordering)
ordering <- pd %>% select(name=label1_pr,order=order1) %>% distinct()

# Use only one triangle (avoid duplicate mirrored pairs)
pd <- pd %>% filter(order2 < order1)

# Symmetric color scale limits based on maximum absolute correlation
max_corr <- max(abs(pd$correlation))

# Apply factor levels to enforce axis ordering in the plot
pd$label1_pr <- factor(pd$label1_pr,levels=ordering %>% arrange(-order) %>% pull(name))
pd$label2_pr <- factor(pd$label2_pr,levels=ordering %>% arrange(order) %>% pull(name))

# =====================================
# Build correlation dot plot
# =====================================
p <- ggplot(pd,aes(x=label1_pr,y=label2_pr))
p <- p + geom_point(aes(fill = correlation,size=abs(correlation)),shape=21)
p <- p + geom_text(aes(label=round(correlation,2)),color="black",vjust=0.5,hjust=0.5,size=1.5)

# Color scale: diverging around zero, bounded by the observed max absolute correlation
p <- p + scale_fill_gradient2(
  low = "#4F9CE6",      # light blue
  mid = "#E6E6E6",      # light grey (not white)
  high = "#F28E2B",     # light orange
  midpoint = 0,
  na.value = "grey90",
  name = "Correlation",
  limits = c(-max_corr, max_corr)
)

# Size scale: encode absolute correlation magnitude
p <- p + scale_size_continuous(
  name = "Correlation\n(absolute)",
  range = c(5, 10),
  limits = c(0, 1),
  breaks = c(0, 0.25, 0.5, 1),
  labels = c("0", "0.25", "0.5", "1")
)      

# Use the same label set on both axes
p <- p + scale_x_discrete("",breaks=ordering$name,labels=ordering$name)
p <- p + scale_y_discrete("",breaks=ordering$name,labels=ordering$name)

# Styling and legend placement
p <- p + theme_cowplot(8) 
p <- p + theme(axis.line=element_blank(),
               plot.margin = margin(t=5,r=5,l=5,b = 5),
               axis.ticks=element_blank(),
               panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
               panel.grid.major=element_line(color="grey80",linewidth=0.2),
                 legend.position =  "top", 
               legend.box = "horizontal",
               legend.box.just = "center" ,  
               legend.text = element_text(size = 5),
               legend.title = element_text(size = 5,face="bold"),
               legend.key.size = unit(0.25, "cm"),
               legend.key.spacing.y = unit(0.1,"cm"),
               legend.box.background = element_rect(fill="white",color = "grey50", linewidth = 0.15) ,
               legend.margin = margin(4,4,4,4)
)

# Display and save figure
p
ggsave(plot=p,filename=outpdf,width=5.5,height=4)
