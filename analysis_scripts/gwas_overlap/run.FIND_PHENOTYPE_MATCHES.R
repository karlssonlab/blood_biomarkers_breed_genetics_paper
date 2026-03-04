# Load required libraries
library(tidyverse)      # dplyr/tidyr pipelines used throughout for data wrangling

# ----- 1. Define significance thresholds -----
sigP <- 5e-8            # genome-wide significance threshold (human GWAS context)
sugP <- 1e-6             # suggestive significance threshold (used to filter dog regions)

# Google Sheet URL for supplemental tables (kept as reference; script reads local TSVs below)
supp_data_url <- "https://docs.google.com/spreadsheets/d/1-fWTNQvS2NeeVVcO1eoI9pIOiuh_NFO9tI5K_wfuDyA/"

# ----- 2. Load precomputed overlap/matching tables (from local TSVs) -----
# Only load if not already present in the environment
if (!exists("raw_gwas_overlap")  | !exists("raw_matching")){
  raw_gwas_overlap <- as_tibble(read.delim("../../data/human_dog_GWAS_intersect/working.GWAS_RAW_OVERLAP.tsv",na.strings=c("NA",""),header=T))
  raw_matching <- as_tibble(read.delim("../../data/human_dog_GWAS_intersect/working.TRAIT_MATCH.tsv",na.strings=c("NA",""),header=T))
}

## =====================================
## Load all data (RDS containing the supplemental raw_* tables)
## =====================================
# Load the RDS only once; then unpack its elements into the global environment
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates objects like raw_S2_ALL_PHENOS, raw_S3_GWAS_REGIONS, etc.
}

# ----- 3. Define the set of blood-related phenotypes to analyze -----
bloodPhenos <- raw_S2_ALL_PHENOS %>% filter(paper_phenotype_category %in% c("Clinical analytes","Plasma metabolites")) %>% 
  select(phenotype,plot_label,standard_name) %>% mutate(plot_label=tolower(plot_label)) %>% distinct()

# ----- 4. Filter dog GWAS regions to suggestive/significant hits within those phenotypes -----
dog_regions <- raw_S3_GWAS_REGIONS %>% filter(region_P <= sugP) %>% 
  inner_join(bloodPhenos %>% select(phenotype) %>% distinct()) %>% 
  select(phenotype,region,CHR,region_start,region_end)  %>% distinct() 

# ----- 5. Build human-overlap peak table from lifted-over coordinates -----
# Start from the raw lift-over overlap table
lifted <- raw_gwas_overlap 

# Keep only records that successfully lifted-over and have a human P-value; compute -log10(P)
overlap <- lifted %>% filter(lifted_over&!is.na(human_P)) %>% mutate(log_human_P=-log10(human_P))  %>% 
  select(phenotype,region,gwas_catalog_trait_name=human_trait_name,log_human_P,human_chr,gwas_pos,PUBMEDID) %>%
  distinct()

# Reduce lifted table down to just the lift status per dog region (drop human-hit detail columns)
lifted <- lifted %>% select(-dog_SNP,-gwas_hit_name,-human_trait_name,-gwas_catalog_trait_name,-human_P,-PUBMEDID,-gwas_pos,-gwas_pos_end) %>% distinct()

# Rank human hits within each dog phenotype-region by strength of association
overlap <- overlap %>% group_by(phenotype,region) %>% mutate(rank=rank(-log_human_P, ties.method = "min"))

# ----- 6. Collapse PubMed IDs per (phenotype, region, trait) at the max signal -----
pubmedids <- overlap %>% select(phenotype,region,gwas_catalog_trait_name,log_human_P,PUBMEDID) %>% 
  distinct() %>% arrange(-log_human_P) %>% group_by(phenotype,region,gwas_catalog_trait_name,PUBMEDID) %>% 
  summarize(log_human_P=max(log_human_P)) %>% ungroup() 
pubmedids <- pubmedids %>% arrange(-log_human_P) %>% group_by(phenotype,region,gwas_catalog_trait_name) %>% 
  summarize(PUBMEDID=paste(PUBMEDID,collapse="|"),log_human_P=max(log_human_P)) %>% ungroup() %>% distinct()

# Reattach the collapsed PubMed IDs and add a convenient chr/pos identifier for the human peak
overlap <- overlap %>% select(-PUBMEDID) %>% distinct() %>% inner_join(pubmedids) 
overlap <- overlap %>% mutate(gwas_chr_pos=paste0("human_",human_chr,"_",gwas_pos)) %>% 
  select(-human_chr,-gwas_pos) %>% mutate(overlap_human_peak=TRUE) 

# ----- 7. Load/shape matching table between dog phenotype names and GWAS Catalog trait names -----
matching <-  raw_matching %>% select(standard_name,gwas_catalog_trait_name=human_trait_name,match_score=score) %>%  distinct()

# Remove generic matching category not intended for analysis
matching <- matching %>% filter(gwas_catalog_trait_name!="metabolite")

# ----- 8. Create a per-region summary table of lift-over and overlap flags -----
# relationship="many-to-many" allows multiple human hits per (phenotype, region)
all <- dog_regions %>% left_join(lifted) %>% left_join(overlap,relationship = "many-to-many") %>% 
  replace_na(list(lifted_over=FALSE,overlap_human_peak=FALSE)) %>% 
  select(phenotype,region,CHR,region_start,region_end,region_P,lifted_over,overlap_human_peak) %>% distinct()

# ----- 9. Attach phenotype metadata and match scoring; classify match type tiers -----
overlap <- bloodPhenos %>% select(phenotype,standard_name,plot_label) %>% 
  distinct() %>% inner_join(overlap) %>% inner_join(matching)

# For each (phenotype, region), keep the best (lowest) match score across associated traits
overlap <- overlap %>% group_by(phenotype,region) %>% summarize(match_score=min(match_score)) %>% inner_join(overlap) %>% distinct() %>% ungroup()

# Convert match_score to a labeled match type
overlap <- overlap %>% mutate(match_type=if_else(match_score==1,"same",
                                                if_else(match_score<=2,"direct",
                                                if_else(match_score<=3,"indirect",
                                                        if_else(match_score<=4,"correlated","unrelated/unknown")))))

# Ensure log_human_P reflects the strongest human association per (phenotype, region, match tier)
overlap <- overlap %>% group_by(phenotype,region,match_score,match_type) %>% summarize(log_human_P=max(log_human_P)) %>% inner_join(overlap)

# ----- 10. Merge overlap detail back into the region-level table and finalize output columns -----
overlap <- all %>% inner_join(overlap)
all <- all %>% left_join(overlap)
all <- all %>% distinct() %>% left_join(overlap)
all <- all %>% select(phenotype,region,CHR,region_start,region_end,region_P,lifted_over,overlap_human_peak,gwas_catalog_trait_name,match_score,match_type,log_human_P,PUBMEDID) %>% 
  distinct() %>% arrange(region_P) %>% mutate(human_P=10**-log_human_P) %>% select(-log_human_P)
