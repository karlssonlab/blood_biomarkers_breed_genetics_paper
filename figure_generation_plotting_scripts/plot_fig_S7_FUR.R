library(tidyverse)      # Data manipulation (dplyr/tidyr) and plotting (ggplot2)
library(cowplot)        # Plot theming and panel layout helpers (theme_cowplot, plot_grid)
library(ggpubr)         # stat_compare_means for p-value annotations on boxplots
source("util.cleaning.R")  # Project-specific helper functions (loaded for shared pipeline consistency)

# Output figure path
outpdf <- "../../../figures_pdf/fig_S7_FUR.pdf"

# Load precomputed supplementary objects if not already present in the session
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # Expose objects like raw_S6_FUR_REGION, raw_S7_CYST_BREED_FREQ, etc.
}

# Plot theme used across panels (small typography, thin lines, compact legend)
theme_dots <- function(){ 
  theme_cowplot(8) %+replace%
    theme(plot.margin = margin(t=5,r=5,l=5,b = 5),
          axis.line = element_line(linewidth = 0.25),    
          axis.ticks = element_line(linewidth = 0.25),    
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          strip.text = element_text(size = 8,hjust = 0,vjust = 0,face = "bold",margin = margin(l=0,t = 5, b = 5),color = "black"),
          strip.background = element_rect(fill = "white",color = NA),
          panel.grid.major.y=element_line(color="grey80",linewidth=0.2),
          panel.grid.major.x=element_blank(),
          
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.text = element_text(size = 3.5),
          legend.title = element_text(size = 5, hjust = 0), # keep small and inline
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.x = unit(0.04, "cm"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          legend.margin = margin(l=5,r=5), 
          legend.box.background = element_rect(fill="white",color = "grey50", linewidth = 0.15) 
    )
}

# Gene highlight intervals to shade on regional association plots
gene_loc <- tibble(
  region  = c("FGF5", "RSPO2"),
  CHR     = c(32, 13),
  xmin_mb = c(35.474935, 8.824734),
  xmax_mb = c(35.494799, 8.970302)
)

# Known literature variants to mark with vertical lines
known_variants_inc_rare <- tibble(
  CHR    = c(32, 32, 32, 32, 32, 13, 13),
  START  = c(35475217, 35475228, 35494495, 35486607, 35475209,  8823200,  8902413),
  END    = c(35475232, 35475229, 35494495, 35486608, 35475210,  8823200,  8902414),
  NAME   = c("PMID23384345_1", "PMID23384345_2", "PMID16879338",
             "PMID23384345_3", "PMID23384345_4", "PMID713490_1", "PMID713490_2"),
  region = c("FGF5", "FGF5", "FGF5", "FGF5", "FGF5", "RSPO2", "RSPO2"),
  phenotype = c("Fur_length","Fur_length","Fur_length","Fur_length","Fur_length","Furnishings","Furnishings")
)
# Midpoint position used for plotting / matching by basepair coordinate
known_variants_inc_rare <- known_variants_inc_rare %>% mutate(BP=floor((START+END)/2))
# Subset to the main focal variants used as plot reference lines
known_variants <- known_variants_inc_rare %>% filter(START==35494495|START==8823200)

# Prepare association results for the two target regions, keeping LD (R2) and P-values
assoc <- raw_S6_FUR_REGION %>% filter(region %in% c("FGF5","RSPO2")) %>% filter(!is.na(P)) %>% 
  select(region,subset,phenotype,SNP,CHR,BP,P,R2=LD.R2)
# Treat missing LD as 0 so binning/plotting is stable
assoc <- assoc %>% replace_na(list(R2=0))
# Human-readable panel titles per phenotype
assoc <- assoc %>% mutate(title=if_else(phenotype=="furnishings","Furnishings (50 dogs)",
                                                      if_else(phenotype=="fur_length","Fur length (50 dogs)",
                                                              if_else(phenotype=="cystathionine","Cystathionine (924 dogs)",phenotype))))
# Hard plot bounds for each locus (in basepairs)
region_limits <- tibble(
  region  = c("FGF5", "RSPO2"),
  CHR     = c(32, 13),
  min = c(35400000,8600000),
  max = c(35600000,9100000)
)

# Restrict association points to the locus bounds
assoc <- assoc %>% inner_join(region_limits) %>% filter(BP>=min&BP<=max)

# Define LD bins for color/size mapping in regional plots
ld_bins <- tibble(
  minLD = c(0, 0.2, 0.5, 0.8),
  maxLD = c(0.2, 0.5, 0.8,1),
  label = c("0-0.2","0.2-0.5","0.5-0.8",">0.8"),
  size = c(1,1.5,1.5,1.75),
  color = c("#313695", "#84335e", "#ae3243", "#d73027")
)

# Identify top cystathionine GWAS signal per region (used for optional reference/annotation)
top_GWAS <- assoc %>% filter(subset=="gwas"&phenotype=="cystathionine"&P<=5e-8) %>% group_by(region) %>% summarize(P=min(P)) %>% inner_join(assoc)
top_GWAS <- top_GWAS %>% select(region,BP,set=phenotype,P) 
top_GWAS <- top_GWAS %>% inner_join(assoc %>% select(region) %>% distinct())

# Create cut breaks covering [0, 1] in the specified LD intervals
breaks <- c(ld_bins$minLD[1], ld_bins$maxLD)

# Bin LD values for discrete aesthetics (color/size), keeping bins even if empty
assoc <- assoc %>%
  mutate(ld_bin = cut(R2,
                      breaks = breaks,
                      labels = ld_bins$label,
                      include.lowest = TRUE, right = FALSE))

# Force the factor order explicitly
  assoc <- assoc %>% mutate(ld_bin = factor(ld_bin, levels = ld_bins$label))

# Track SNPs that are genome-wide significant for cystathionine (used to set point shape)
  cyst_sig <- assoc %>% filter(subset=="gwas"&phenotype=="cystathionine"&P<=5e-8) %>% select(SNP)

# Regional association plot for GWAS subset (colored by LD bin, cystathionine-significant points shaped)
plot_region_gwas <- function(in_region,plot_start,plot_end) {
    
  
    indata  <- assoc  %>% filter(region == in_region&subset=="gwas") %>% select(title,CHR,SNP,phenotype,BP,P,R2,ld_bin)
    indata <- indata %>% filter(BP>=plot_start&BP<=plot_end)
    # Use the smallest P in the region to set the y-axis limit (rounded up)
    ymax_all <- ceiling(-log10(min(assoc %>% filter(region == in_region) %>% pull(P),na.rm=TRUE)))
                               
    
    xmin <- plot_start/1e6
    xmax <- plot_end/1e6
    
  # Shaded gene interval and (optional) top GWAS signal for the region
  shading <- gene_loc %>% filter(region == in_region)
  vlines <- top_GWAS %>% filter(region == in_region)

  
  # Derive chromosome label for x-axis (assumes single chromosome per region)
  chr_val <- unique(indata$CHR)
  chr_val <- chr_val[!is.na(chr_val)][1]
  xlab <- paste0("chromosome ", chr_val, " (Mb)")
  # Arrange so high-LD and low-P points tend to plot on top
  indata <- indata %>% arrange(R2,-P) 
  # Mark cystathionine-significant SNPs for shape mapping
  indata <- indata %>% mutate(sig_cyst=if_else(SNP %in% cyst_sig$SNP,TRUE,FALSE))
  
  shaded_start_y <- ymax_all-(ymax_all/40)
  shaded_end_y <- ymax_all-(ymax_all/10)
  # Small offset for placing gene label just right of the shaded interval
  increment <- (max(indata$BP/1e6)-min(indata$BP/1e6))/100
  indata$ld_bin <- factor(indata$ld_bin,levels=ld_bins$label)
  
  p <- ggplot(indata, aes(x = BP/1e6, y = -log10(P)))
  # Mark known published variants with vertical lines
  p <- p + geom_vline(aes(xintercept = BP/1e6),color="#a50f15",alpha=0.5,linewidth=0.25,data=known_variants %>% filter(region==in_region) )
  # Shade the gene interval and draw a top segment plus label
  p <- p + annotate("rect",xmin = shading$xmin_mb[[1]], xmax = shading$xmax_mb[[1]],ymin = -Inf, ymax =ymax_all+0.2,alpha = 0.25, fill = "#fcbba1")
  p <- p + annotate("segment",x = shading$xmin_mb[[1]], xend = shading$xmax_mb[[1]],y=ymax_all+0.2, yend = ymax_all+0.2,color="#a50f15")
  p <- p + geom_text(label=shading$region[[1]],color="#a50f15",lineheight=0.8,x=shading$xmax_mb[[1]]+increment,y=ymax_all,vjust=0.5,hjust=0,size=1.75)
  # Genome-wide significance threshold
  p <- p + geom_hline(yintercept = -log10(5e-8),color = "grey30", linetype = 2, linewidth = 0.25)
  # Map LD bins to color/size; shape indicates cystathionine GWAS significance
  p <- p + geom_point(aes(color = ld_bin,shape=sig_cyst,size=ld_bin),alpha = 0.8)
  p <- p + ggtitle(unique(indata$title))
  ld_breaks <- ld_bins$label
  p <- p +
    scale_color_manual(
      values = setNames(ld_bins$color, ld_bins$label),
      breaks = ld_breaks,
      drop = FALSE,
      name = "LD (r²)"
    ) +
    scale_size_manual(
      values = setNames(ld_bins$size, ld_bins$label),
      breaks = ld_breaks,
      drop = FALSE,
      name = "LD (r²)"
    ) 
  p <- p + scale_shape_manual(name=expression(Cystathionine~p < 5 %*% 10^{-8}),values=c(16,17))
  
  
  # Use chromosome label as x-axis title and plot window bounds in Mb
  p <- p + scale_x_continuous(name = xlab,limits=c(xmin,xmax))
  p <- p + scale_y_continuous("-log10(p)",limits=c(0,ymax_all+0.3),breaks=c(0:3)*5,expand = expansion(mult = c(0.05,0.05)))
  # Hide legend in this panel (legends handled at grid level if needed)
  p <- p + theme_dots() + theme(legend.position="none") + guides(color = guide_legend(nrow = 2))
  p
  return(p)
}

# Default region window values used during development (not required by the functions)
in_region <- "RSPO2"
plot_start=8800000
plot_end=9000000

# Regional association plot for the 50-dog precision subset (same styling as GWAS panel)
plot_region_panel <- function(in_region,plot_start,plot_end) {

  indata  <- assoc  %>% filter(region == in_region&subset=="precision_50_dogs")
  indata <- indata %>% filter(BP>=plot_start&BP<=plot_end)
  # Use the smallest P in the region to set the y-axis limit (rounded up)
  ymax_all <- ceiling(-log10(min(assoc %>% filter(region == in_region) %>% pull(P),na.rm=TRUE)))
  
  xmin <- plot_start/1e6
  xmax <- plot_end/1e6
  
  # Shaded gene interval and (optional) top GWAS signal for the region
  shading <- gene_loc %>% filter(region == in_region)
  vlines <- top_GWAS %>% filter(region == in_region)

    # Derive chromosome label for x-axis (assumes single chromosome per region)
    chr_val <- unique(indata$CHR)
  chr_val <- chr_val[!is.na(chr_val)][1]
  xlab <- paste0("chromosome ", chr_val, " (Mb)")
  
  # Mark cystathionine-significant SNPs for shape mapping
  indata <- indata %>% mutate(sig_cyst=if_else(SNP %in% cyst_sig$SNP,TRUE,FALSE))
  # Small offset for placing gene label just right of the shaded interval
  increment <- (max(indata$BP/1e6)-min(indata$BP/1e6))/100
  indata$ld_bin <- factor(indata$ld_bin,levels=ld_bins$label)
  p <- ggplot(indata, aes(x = BP/1e6, y = -log10(P)))
  # Mark known published variants with vertical lines
  p <- p + geom_vline(aes(xintercept = BP/1e6),color="#a50f15",alpha=0.5,linewidth=0.25,data=known_variants %>% filter(region==in_region) )
  # Shade the gene interval and draw a top segment plus label
  p <- p + annotate("rect",xmin = shading$xmin_mb[[1]], xmax = shading$xmax_mb[[1]],ymin = -Inf, ymax =ymax_all+0.2,alpha = 0.25, fill = "#fcbba1")
  p <- p + annotate("segment",x = shading$xmin_mb[[1]], xend = shading$xmax_mb[[1]],y=ymax_all+0.2, yend = ymax_all+0.2,color="#a50f15")
  p <- p + geom_text(label=shading$region[[1]],color="#a50f15",lineheight=0.8,x=shading$xmax_mb[[1]]+increment,y=ymax_all,vjust=0.5,hjust=0,size=1.75)
  # Genome-wide significance threshold
  p <- p + geom_hline(yintercept = -log10(5e-8),color = "grey30", linetype = 2, linewidth = 0.25)
  # Map LD bins to color/size; shape indicates cystathionine GWAS significance
  p <- p + geom_point(aes(color = ld_bin,shape=sig_cyst,size=ld_bin),alpha = 0.8)
  p <- p + ggtitle(unique(indata$title))
  ld_breaks <- ld_bins$label
  p <- p +
    scale_color_manual(
      values = setNames(ld_bins$color, ld_bins$label),
      breaks = ld_breaks,
      drop = FALSE,
      name = "LD (r²)"
    ) +
    scale_size_manual(
      values = setNames(ld_bins$size, ld_bins$label),
      breaks = ld_breaks,
      drop = FALSE,
      name = "LD (r²)"
    ) 
  p <- p + scale_shape_manual(expression(Cyst.~p < 5 %*% 10^{-8}),values=c(16,17))

  
  # Use chromosome label as x-axis title and plot window bounds in Mb
  p <- p + scale_x_continuous(name = xlab,limits=c(xmin,xmax))
  p <- p + scale_y_continuous("-log10(p)",limits=c(0,ymax_all+0.3),breaks=c(0:3)*5,expand = expansion(mult = c(0.05,0.05)))
  p <- p + theme_dots() + guides(color = guide_legend(nrow = 2))

    
  return(p)
}

# Boxplot of breed-level allele frequency by phenotype class for a given SNP/trait
make_maf_boxplot <- function(in_snp, in_trait) {
  pd <- cyst_breed %>% filter(trait == in_trait & SNP == in_snp)
  # Enforce desired phenotype order (only for levels present in the data)
  order <- tibble(phenotype = c("short","no","long","yes"), order = 1:4) %>%
    inner_join(pd %>% distinct(phenotype), by = "phenotype") %>%
    arrange(order)
  pd$phenotype <- factor(pd$phenotype, levels = order$phenotype)
  # Flip allele frequency for this specific SNP so the plotted allele matches the intended direction
  if (in_snp=="chr32:35496522:G:A"){
    pd <- pd %>% mutate(MAF=1-MAF)
  }
  # Build x-axis labels with counts and mean +/- SD per phenotype group
  xlabels <- pd %>%
    group_by(phenotype) %>%
    summarize(n=n(),mean=mean(MAF),sd=sd(MAF)) %>% 
    mutate(label = if_else(
      phenotype == "yes", in_trait,
      if_else(phenotype == "no", paste(phenotype, in_trait),
              if_else(phenotype %in% c("short","long"), paste(phenotype, "fur"), in_trait)
      ))) %>%
    mutate(label = paste0(label, "\n(", n, " breeds)\n",round(mean,2)," +/- ",round(sd,2)))
  
  # Color mapping for phenotype groups
  col_map <- c("no" = "grey30", "yes" = "#a50f15", "medium/long"="#a50f15","short" = "grey30", "long" = "#a50f15")

  
  p <- ggplot(pd, aes(x = phenotype, y = MAF))
  p <- p + geom_boxplot(aes(fill = phenotype, color = phenotype),
                        outlier.shape = NA, width = 0.3, linewidth = 0.25, alpha = 0.25)
  p <- p + geom_jitter(aes(fill = phenotype, color = phenotype),
                       shape = 16, width = 0.1, height = 0, alpha = 0.5)
  
   p <- p + scale_color_manual(name = "phenotype", values = col_map, drop = FALSE)
  p <- p + scale_fill_manual(name = "phenotype",  values = col_map, drop = FALSE)
  
  # Title shown only for furnishings panels (to keep layout cleaner for fur length)
  if(in_trait=="furnishings"){
  p <- p + ggtitle(in_snp)
  } else {
    p <- p + ggtitle("")
  }
  
  # Two-group comparison using Wilcoxon rank-sum test, with formatted p-value label
  p <- p + stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size=2,
    hide.ns = TRUE,
    label.x=1.5,
    label.y=1.05,
    hjust=0.5  # near top of each facet
  )
  p <- p + scale_x_discrete("phenotype",
                            labels = xlabels$label,
                            breaks = xlabels$phenotype)
  p <- p + scale_y_continuous("allele frequency",breaks=c(0,0.5,1),expand = expansion(mult = c(0.1,0.15)))
  
  
  # Boxplot-specific theme tweaks
  p <- p + theme_cowplot(8)
  p <- p + theme(
    plot.title = element_text(size=6,face="bold"),
    axis.line.y = element_line(linewidth = 0.25),
    axis.line.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.text    = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 6, hjust = 0, vjust = 0, face = "bold",
                              margin = margin(t = 2, b = 5), color = "black"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25),
    legend.position = "none"
  )
}

# Build the four regional panels (GWAS + precision for each locus) and arrange into a 2x2 grid
gwasFGF <- plot_region_gwas(in_region  = "FGF5",35450000,35550000)
gwasRSPO2 <- plot_region_gwas(in_region  = "RSPO2",plot_start=8600000,plot_end=9100000)

panelFGF <- plot_region_panel(in_region  = "FGF5",35450000,35550000)
panelRSPO2 <- plot_region_panel(in_region  = "RSPO2",plot_start=8600000,plot_end=9100000)
grid_row1_2 <- plot_grid(gwasRSPO2,gwasFGF,panelRSPO2,panelFGF,nrow=2,labels=LETTERS,label_size = 8,rel_heights = c(1,1),label_x = 0.01,label_y = 1,hjust = 0)

# Prepare breed-level phenotype labels and cystathionine allele-frequency data for boxplots
breed_phenos <- raw_S12_BREED_PHENOS %>% select(breed,furnishings,length,death=`median_age_death.(McMillan.2024)`) %>% distinct()
cyst_breed <- raw_S7_CYST_BREED_FREQ %>% mutate(breed=str_replace_all(CLST,"_"," ")) %>% inner_join(breed_phenos)
cyst_breed <- cyst_breed %>% select(-CLST) %>% pivot_longer(c(furnishings,length)) %>% rename(trait=name,phenotype=value)
# Keep only clean binary phenotype calls used in comparisons
cyst_breed <- cyst_breed %>% filter(!is.na(phenotype)&phenotype!="variable"&phenotype!="medium") 

# Save final multi-panel figure
ggsave(plot=grid_row1_2,filename=outpdf,width=7,height=6)
