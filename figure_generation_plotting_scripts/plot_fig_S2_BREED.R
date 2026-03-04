library(tidyverse)
library(cowplot)
library(ggrepel)

# Project utilities (data standardization + helper transforms)
source("util.cleaning.R")

# Output figure path
outpdf <- "../../../figures_pdf/fig_S2_BREED.pdf"

## =====================================
## Load all data (only once per session)
## =====================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_* objects in the global environment
}

## =====================================
## Standardize breed naming across input tables
## =====================================
anc <- standardize_breed_names(raw_S15_BREED_ANC, pop)
anc <- anc %>% rename(breed=pop)

breedsizes <- standardize_breed_names(raw_S12_BREED_PHENOS, breed)

# Precompute counts of sequenced dogs by single-breed status (used in axis labels later)
dog_counts_dap <- raw_S14_SEQ_DOGS %>% group_by(single_breed) %>% count()

# Main analysis table: standardized breed field + convert yes/no fields to logical
indata <- standardize_breed_names(raw_S14_SEQ_DOGS, standardized_breed) %>% rename(breed=standardized_breed)
indata <- convert_yesno_to_logical(indata)

## =====================================
## Prepare reference tables used in multiple panels
## =====================================
# Darwin's Ark breed fractions (single vs mixed) for later comparison panel
darwinsArk <- breedsizes %>% select(breed,frac_single_breed=DarwinsArk_frac_single,frac_mixed_breed=DarwinsArk_frac_mixed)

# Keep only size-group annotation per breed; ensure any missing size group gets "Variable"
breedsizes <- anc %>% select(breed) %>% distinct() %>% full_join(breedsizes) %>% replace_na(list(Size.Group="Variable")) %>% select(breed,Size.Group)

## =====================================
## Plot theme helpers (categorical-y barplots; quantitative scatterplots)
## =====================================
theme_categ_y <- function(){
  theme_cowplot(12) %+replace%
    theme(axis.title = element_text(size=7),
          axis.text = element_text(size=6),
          panel.grid.major.x=element_line(color="grey80",linewidth=0.2),
          legend.position = c(0.98,0.02),
          legend.text = element_text(size=5.5),legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.justification.inside = c(1, 0),
          legend.box.just = "right",
          axis.ticks.y=element_blank(),
          strip.text.x=element_text(size=7,hjust=0.5,face="bold",lineheight=1), margin(5,3,5,3)
    )
}
theme_quant <- function(){
  theme_cowplot(12) %+replace%
    theme(axis.title = element_text(size=7),axis.text = element_text(size=6),
            legend.position = c(0.98,0.98),
          legend.text = element_text(size=5.5),legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.justification.inside = c(1, 1),
          legend.box.just = "right",
          legend.background = element_rect(color = "grey20",fill = "white",linewidth = 0.25),
          legend.margin = margin(2,2,2,2),
          strip.text.x=element_text(size=7,hjust=0.5,face="bold",lineheight=1), margin(5,3,5,3)
    )
}

## =====================================
## Panel A: Top single-breed dogs (top N breeds by fraction)
## =====================================
ntot <- length(indata %>% filter(single_breed) %>% pull(dog_id))
nbreeds_plotted <- 10

pd <- indata
pd <- pd %>% filter(single_breed) %>% inner_join(breedsizes) %>% group_by(breed,Size.Group) %>% count() %>% mutate(frac=n/ntot)
pd <- pd %>% ungroup() %>% arrange(-frac) %>% mutate(rank=row_number())

# Order factors so plot sorts by fraction
breeds <- pd %>% arrange(-rank) %>% pull(breed)
pd$breed <- factor(pd$breed,levels=breeds)
pd$Size.Group <- factor(pd$Size.Group,levels=c("Variable","Toy and Small","Medium","Standard","Large","Giant"))

# Shorten one label for better y-axis readability
yaxis <- pd %>% select(breed) %>% mutate(label=str_replace(breed,"cavalier","cav. "))

p <- ggplot(pd %>% filter(rank <= nbreeds_plotted),aes(x=frac,y=breed)) +
  geom_bar(aes(fill=Size.Group),alpha=1,color="#252525",linewidth=0.25,stat="identity",width=0.75)
p <- p + scale_fill_manual(values=c("Toy and Small"="#8dd3c7","Medium"="#ffffb3","Standard"="#bebada","Large"="#fb8072","Giant"="#80b1d3","Variable"="#bdbdbd"))
p <- p + scale_y_discrete("",breaks=yaxis$breed,labels=yaxis$label)
p <- p + scale_x_continuous("fraction of single-breed dogs",breaks=c(0,0.05,0.1),limits=c(0,0.11))
p <- p + theme_categ_y() + theme(axis.text.y = element_text(size=6))
p_breed_top <- p

# Store fractions for later single-vs-mixed comparison
pd_single_vs_mixed <- pd %>% select(breed,frac_single=frac)
rm(pd)

## =====================================
## Panel B: Top breeds contributing to mixed ancestry (top N by total ancestry fraction)
## =====================================
pd <- indata
# Total number of mixed-breed dogs with ancestry annotations present in anc
ntot <- pd %>% filter(!single_breed) %>% select(dog_id) %>% distinct() %>% filter(dog_id %in% anc$dog_id) %>% count() %>% pull(n)

# Expand mixed dogs to (dog_id, breed, pct) rows using ancestry table
pd <- pd %>% filter(!single_breed)  %>% select(dog_id) %>% inner_join(anc)  %>% left_join(breedsizes)

# Convert ancestry percent to overall fraction across all mixed-breed ancestry
tot <- sum(pd$pct)
pd <- pd %>% group_by(breed) %>% summarize(frac=sum(pct)/tot) %>% arrange(desc(frac))
pd <- pd %>% left_join(breedsizes)
pd <- pd %>% ungroup() %>% arrange(-frac) %>% mutate(rank=row_number())
pd_all_anc <- pd

# Order factors so plot sorts by fraction
breeds <- pd %>% arrange(-rank) %>% pull(breed)
pd$breed <- factor(pd$breed,levels=breeds)
pd$Size.Group <- factor(pd$Size.Group,levels=c("Variable","Toy and Small","Medium","Standard","Large","Giant"))

p <- ggplot(pd %>% filter(rank <= nbreeds_plotted),aes(x=frac,y=breed)) +
  geom_bar(aes(fill=Size.Group),alpha=1,color="#252525",linewidth=0.25,stat="identity",width=0.75)
p <- p + scale_fill_manual(values=c("Toy and Small"="#8dd3c7","Medium"="#ffffb3","Standard"="#bebada","Large"="#fb8072","Giant"="#80b1d3","Variable"="#bdbdbd"))
p <- p + scale_y_discrete("")
p <- p + scale_x_continuous("fraction of mixed breed dog ancestry",breaks=c(0,0.05,0.1),limits=c(0,0.11))
p <- p + theme_categ_y() + theme(axis.text.y = element_text(size=6))
p_mixed_top <- p

# Merge single- and mixed- fractions by breed for within-DAP comparison panel
pd_single_vs_mixed <- pd %>% select(breed,frac_mixed=frac) %>% inner_join(pd_single_vs_mixed)

rm(p)

## =====================================
## Panel D: Breed frequency correlation (Darwin's Ark vs DAP single-breed)
## =====================================
ntot_single <- indata %>% filter(single_breed) %>% select(dog_id) %>% distinct() %>% count() %>% pull(n)

# Build breed fractions for DAP single-breed and bind Darwin's Ark fractions
pdSource <- indata %>% filter(single_breed) %>% select(dog_id,breed) %>% group_by(breed) %>% count()
pdSource <- pdSource %>% mutate(source="seqDAP") %>% mutate(frac=n/ntot_single) %>% select(-n)
pdSource <- pdSource %>% bind_rows(darwinsArk %>% select(breed,frac=frac_single_breed) %>% mutate(source="DA") )
pdSource <- pdSource %>% left_join(breedsizes %>% select(breed,Size.Group) %>% distinct())
pdSource <- pdSource %>% filter(!is.na(Size.Group)) %>% distinct()

# Wide table: one row per breed with DA and seqDAP fractions
pd <- pdSource %>% select(source,Size.Group,breed,frac) %>% pivot_wider(names_from=source,values_from=frac)
pd <- pd %>% replace_na(list(DA=0,seqDAP=0,allDAP=0))
pd$Size.Group <- factor(pd$Size.Group,levels=c("Variable","Toy and Small","Medium","Standard","Large","Giant"))
pd <- pd %>% filter(!is.na(Size.Group))

# Spearman correlation + permutation test for p-value
ct <- cor.test(pd$DA, pd$seqDAP,method="spearman")
cor_val <- round(ct$estimate, 3)
p_val   <- format(ct$p.value, scientific = TRUE, digits = 3)
n_perm <- 100000
if (!exists("perm_stats_D", inherits = FALSE) ||
    length(perm_stats_D) != n_perm) {
  perm_stats_D <- replicate(n_perm, {
    cor(pd$DA, sample(pd$seqDAP), method = "spearman", use = "complete.obs")
  })
}
# Two-sided permutation count beyond observed magnitude
b <- sum(abs(perm_stats_D) >= abs(ct$estimate))

# Bias-corrected permutation p-value
p_perm <- (b + 1) / (length(perm_stats_D) + 1)
p_perm <- formatC(signif(p_perm, 2), format = "e", digits = 0)
label_text <- paste0("Spearman's rho = ", cor_val, "\nperm. p < ",p_perm)

p <- ggplot(pd,aes(x=DA,y=seqDAP))
p <- p + geom_abline(linetype=2,linewidth=0.25)
p <- p + geom_point(aes(color=Size.Group,fill=Size.Group,shape=Size.Group),size=1.25)

# Label only the most frequent breeds to reduce clutter
p <- p + ggrepel::geom_text_repel(
  data = pd %>% filter(seqDAP > 0.05|DA>0.05),
  aes(label = str_replace(breed, " ", "\n"), color = Size.Group),
  size = 2,
  lineheight = 0.9,
  hjust = 1,              # right-justify text
  nudge_x = -0.005,        # push labels to the left of points
  direction = "y",        # only repel vertically after nudging
  box.padding = 0.3,
  point.padding = 0.1,
  max.overlaps = Inf,
  segment.size = 0.2,
  show.legend = FALSE
)
# Place correlation annotation in top-left of panel
p <- p +  annotate("text", x = 0, y = Inf, label = label_text, hjust = 0, vjust = 1.5, size = 2.25,lineheight=1)
p <- p + scale_color_manual(values=c("Toy and Small"="#1b9e77","Medium"="#c28f00","Standard"="#7570b3","Large"="#d95f02","Giant"="#1f78b4","Variable"="#4d4d4d"))
p <- p + scale_fill_manual(values=c("Toy and Small"="#1b9e77","Medium"="#c28f00","Standard"="#7570b3","Large"="#d95f02","Giant"="#1f78b4","Variable"="#4d4d4d"))
p <- p + scale_shape_manual(values = c(16, 17, 15, 18, 3, 8))
p <- p + scale_y_continuous(paste0("Dog Aging Project single-breed (N=",(dog_counts_dap %>% filter(single_breed) %>% pull(n))[1],")"),limits=c(0,0.12),breaks=c(0,0.04,0.08,0.12))
p <- p + scale_x_continuous("Darwin's Ark single-breed (N=934)",limits=c(0,0.12),breaks=c(0,0.04,0.08,0.12))
p <- p + theme_quant()
pcompare_single <- p

## =====================================
## Panel C: Within-DAP correlation (single-breed fraction vs mixed-ancestry fraction)
## =====================================
pd <- pd_single_vs_mixed %>% left_join(breedsizes)
pd$Size.Group <- factor(pd$Size.Group,levels=c("Variable","Toy and Small","Medium","Standard","Large","Giant"))
pd <- pd %>% filter(!is.na(Size.Group))

# Spearman correlation + permutation test for p-value
ct <- cor.test(pd$frac_single, pd$frac_mixed,method="spearman")
cor_val <- round(ct$estimate, 3)
p_val   <- format(ct$p.value, scientific = TRUE, digits = 3)
cor_val <- round(ct$estimate, 3)
p_val   <- format(ct$p.value, scientific = TRUE, digits = 3)
n_perm <- 100000
if (!exists("perm_stats_C", inherits = FALSE) ||
    length(perm_stats_C) != n_perm) {
  perm_stats_C <- replicate(n_perm, {
    cor(pd$frac_single, sample(pd$frac_mixed), method = "spearman", use = "complete.obs")
  })
}
# Two-sided permutation count beyond observed magnitude
b <- sum(abs(perm_stats_C) >= abs(ct$estimate))

# Bias-corrected permutation p-value
p_perm <- (b + 1) / (length(perm_stats_C) + 1)
p_perm <- formatC(signif(p_perm, 2), format = "e", digits = 0)
label_text <- paste0("Spearman's rho = ", cor_val, "\nperm. p < ",p_perm)

p <- ggplot(pd,aes(x=frac_mixed,y=frac_single))
p <- p + geom_abline(linetype=2,linewidth=0.25)
p <- p + geom_point(aes(color=Size.Group,fill=Size.Group,shape=Size.Group),size=1.25)

# Label only the most frequent breeds to reduce clutter
p <- p + geom_text_repel(
  data = pd %>% filter(frac_single > 0.05 | frac_mixed > 0.05),
  aes(label = str_replace(breed, " ", "\n"), color = Size.Group),
  size = 2,
  lineheight = 0.9,
  box.padding = 0.3,
  point.padding = 0.1,
  max.overlaps = Inf,
  segment.size = 0.2, show.legend = FALSE
)
# Place correlation annotation in top-left of panel
p <- p +  annotate("text", x = 0, y = Inf, label = label_text, hjust = 0, vjust = 1.5, size = 2.25,lineheight=1)
p <- p + scale_color_manual(values=c("Toy and Small"="#1b9e77","Medium"="#c28f00","Standard"="#7570b3","Large"="#d95f02","Giant"="#1f78b4","Variable"="#4d4d4d"))
p <- p + scale_fill_manual(values=c("Toy and Small"="#1b9e77","Medium"="#c28f00","Standard"="#7570b3","Large"="#d95f02","Giant"="#1f78b4","Variable"="#4d4d4d"))
p <- p + scale_shape_manual(values = c(16, 17, 15, 18, 3, 8))
 p <- p + scale_y_continuous(paste0("Dog Aging Project single-breed (N=",(dog_counts_dap %>% filter(single_breed) %>% pull(n))[1],")"),limits=c(0,0.12),breaks=c(0,0.04,0.08,0.12))
p <- p + scale_x_continuous("Dog Aging Project mixed-breed (N=3485)",limits=c(0,0.12),breaks=c(0,0.04,0.08,0.12))
p <- p + theme_quant()
pcompare_withinDAP <- p
p

## =====================================
## Combine panels and write to PDF
## =====================================
row1 <- plot_grid(p_breed_top,p_mixed_top,labels=LETTERS,nrow=1,label_size=10)
row2 <- plot_grid(pcompare_withinDAP,pcompare_single,labels=c("C","D"),nrow=1,label_size=10)

grid <- plot_grid(row1,row2, nrow=2,rel_heights=c(1,1.3))
ggsave(plot=grid, filename=outpdf, width=6, height=5)
