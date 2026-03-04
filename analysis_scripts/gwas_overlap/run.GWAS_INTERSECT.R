# ================================
# Libraries
# ================================
# Core data manipulation + string/purrr helpers used throughout
library(tidyverse) 
# Reading from Google Sheets (not used in this script)
library(cowplot)

# ================================
# 1) Analysis thresholds
# ================================
sigP <- 5e-8
sugP <- 1e-6 

# ================================
# 2) Input file paths
# ================================
### Can be downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads
in_gwas_catalog <- "../../data/human_dog_GWAS_intersect/gwas_catalog_v1.0-associations_e113_r2025-02-08.tsv"

# ================================
# 3) Load saved supplemental R objects (raw_S2_ALL_PHENOS, raw_S3_GWAS_REGIONS, etc.)
# ================================
if (!exists("dat")) {
  dat <- readRDS("../../data/DAP_supp_data.rds")
  list2env(dat, .GlobalEnv)   # creates raw_* objects in the global environment
}

# ================================
# 4) Load curated list of GWAS Catalog traits to include
# ================================
if (!exists("raw_human_trait_list")){
  raw_human_trait_list <- as_tibble(read.delim("../../data/human_dog_GWAS_intersect/working.HUMAN_TRAITS.tsv",na.strings=c("NA",""),header=T))
}

# ================================
# 5) Define dog blood phenotypes and extract dog GWAS regions ("clumps") for those phenotypes
# ================================
bloodPhenos <- raw_S2_ALL_PHENOS %>% filter(paper_phenotype_category %in% c("Clinical analytes","Plasma metabolites")) %>% select(phenotype,plot_label,paper_phenotype_category) %>% distinct()
bloodPhenos <- bloodPhenos %>% mutate(plot_label=tolower(plot_label))
bloodPhenos <- raw_S3_GWAS_REGIONS  %>% select(phenotype) %>% distinct() %>% inner_join(bloodPhenos)

# All blood phenotype GWAS hits for dog
gwas_dog <- raw_S3_GWAS_REGIONS %>% inner_join(bloodPhenos)
dim(gwas_dog)

# ================================
# 6) Output paths for intermediate BED/TSV files
# ================================
working_dir <- "../../data/human_dog_GWAS_intersect/"
if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

out_dogClumps <- paste0(working_dir, "dog.blood_clumps.tsv")
out_dogClumps_bedFormat <- paste0(working_dir, "dog.blood_clumps.bed")
out_dogClumps_onHuman_lifted_bed <- paste0(working_dir, "dog.blood_clumps.onHuman.mapped.bed")
out_dogClumps_onHuman_unmapped_bed <- paste0(working_dir, "dog.blood_clumps.onHuman.unmapped.bed")
out_dogClumps_onHuman_lifted_500kb_bed <- paste0(working_dir, "dog.blood_clumps.onHuman.mapped.500kb.bed")
out_humanGWAS_all_bed  <- paste0(working_dir, "human.GWAS.25_pubs.sigP.all.bed")
out_humanGWAS_dogClumps.bedtools <- paste0(working_dir, "human.GWAS.25_pubs__dog.blood_clumps.overlap.sugP.bedtools.out")

# ================================
# 7) Process dog GWAS clumps: define regions, attach SNPs within each region, and build locus IDs
# ================================
dog_clumps <- gwas_dog %>%
  select(region, CHR, region_start, region_end, phenotype, region_P) %>%
  distinct()

# Keep SNPs at least as significant as the regional lead P (per-region filter)
dog_snps <- gwas_dog %>%
  filter(P <= region_P) %>%
  select(region, CHR, region_start, region_end, SNP, region_P) %>%
  distinct()

dog_clumps <- dog_clumps %>% 
  left_join(dog_snps, by = c("region","CHR","region_start","region_end","region_P"))

# Standardize column names and create a unique locus identifier string
dog_clumps <- dog_clumps %>% rename(dog_trait = phenotype, dog_region=region, dog_start = region_start, dog_end = region_end, dog_P = region_P) %>%
  mutate(dog_CHR = paste0("chr", CHR),
         dog_input_region = paste0(dog_region,"_",dog_CHR, ":", dog_start, "..", dog_end),
         dog_locus = paste(dog_input_region, dog_trait, sep = "__")) %>%
  distinct()

# Ensure each dog region is at least 1kb (liftOver/bedtools are more stable with non-trivial intervals)
dog_clumps <- dog_clumps %>%
  mutate(wide_start = if_else(dog_end - dog_start < 1000, floor((dog_start + dog_end) / 2) - 500, dog_start),
         wide_end = if_else(dog_end - dog_start < 1000, ceiling((dog_start + dog_end) / 2) + 500, dog_end),
         dog_start = wide_start, dog_end = wide_end) %>%
  select(-wide_start, -wide_end)
all_dog_clumps <- dog_clumps

# ================================
# 8) Write dog clumps to disk (TSV for metadata + BED for liftover)
# ================================
tmp_dog_clumps <- dog_clumps
# Write a human-readable TSV with spaces in headers (restored later with check.names = FALSE)
colnames(tmp_dog_clumps) <- gsub("_", " ", colnames(tmp_dog_clumps))
write.table(tmp_dog_clumps, out_dogClumps, row.names = FALSE, quote = FALSE, sep = "	")

# BED format: chr, start, end, name (no header)
write.table(dog_clumps %>% select(dog_CHR, dog_start, dog_end, dog_locus) %>% distinct(),
            out_dogClumps_bedFormat, row.names = FALSE, quote = FALSE, sep = "	", col.names = FALSE)

# ================================
# 9) LiftOver dog intervals to human coordinates, then widen to 500kb windows
# ================================
system(paste("/Applications/liftOver/liftOver -minMatch=0.1", out_dogClumps_bedFormat, 
             "../../data/human_dog_GWAS_intersect/canFam4ToHg38.over.chain.gz", 
             out_dogClumps_onHuman_lifted_bed, out_dogClumps_onHuman_unmapped_bed))

dog_liftOver_to_human <- as_tibble(read.delim(out_dogClumps_onHuman_lifted_bed, header = FALSE)) %>%
  rename(human_chr = V1, human_start = V2, human_end = V3, dog_locus = V4)

# Expand any lifted interval shorter than 500kb to a 500kb window centered on the interval midpoint
# Also coerce coordinates to non-scientific character strings for downstream BED writing
dog_liftOver_to_human <- dog_liftOver_to_human %>%
  mutate(wide_start = if_else(human_end - human_start < 500000, floor((human_start + human_end) / 2) - 250000, human_start),
         wide_end = if_else(human_end - human_start < 500000, ceiling((human_start + human_end) / 2) + 250000, human_end),
         human_start = str_trim(format(human_start, scientific = FALSE)),
         human_end = str_trim(format(human_end, scientific = FALSE))) %>%
  select(human_chr, human_start, human_end, dog_locus)

# BED for bedtools intersect (no header)
write.table(dog_liftOver_to_human, out_dogClumps_onHuman_lifted_500kb_bed, row.names = FALSE, quote = FALSE, sep = "	", col.names = FALSE)

# ================================
# 10) Process Human GWAS Catalog: parse positions and write significant hits to BED
# ================================
# Convert GWAS Catalog CHR_POS strings into numeric start/end bounds
# - If CHR_POS has multiple positions separated by ';', take min/max
# - Otherwise treat as a single base and create a 1bp interval [pos, pos+1)
get_bounds <- function(x) {
  if (str_detect(x, ";")) {
    nums <- sort(as.numeric(str_split(x, ";")[[1]]))
    c(nums[1], nums[length(nums)])
  } else {
    xnum <- as.numeric(x)
    c(xnum, xnum + 1)
  }
}

# Read the human GWAS catalog (cached in-memory if already loaded)
if (!exists("raw_gwas_human")){
  raw_gwas_human <- read_tsv(in_gwas_catalog,show_col_types = FALSE)
}

# Keep only genome-wide significant hits with valid chromosome and position information
gwas_human <- raw_gwas_human %>% filter(`P-VALUE` <= sigP & !is.na(CHR_ID) & !is.na(CHR_POS)) %>% select(-SNP_ID_CURRENT) 
gwas_human <- gwas_human %>% rename(gwas_catalog_trait_name=`DISEASE/TRAIT`)

# Limit to curated set of GWAS Catalog traits to include
human_trait_list <- raw_human_trait_list %>% filter(include) %>% select(gwas_catalog_trait_name,canonical_name) %>% distinct()
gwas_human <- gwas_human %>% inner_join(human_trait_list %>% select(gwas_catalog_trait_name) %>% distinct(),relationship = "many-to-many")

# Parse CHR_POS into [start,end) interval coordinates
gwas_human <- gwas_human %>%
  mutate(bounds = map(CHR_POS, get_bounds)) %>%
  mutate(
    gwas_pos     = map_dbl(bounds, 1),
    gwas_pos_end = map_dbl(bounds, 2)
  ) %>%
  select(-bounds)

# Build a unique hit name for overlap bookkeeping
gwas_human <- gwas_human %>%
  rename(human_CHR = CHR_ID, human_P = `P-VALUE`) %>%
  mutate(human_CHR = paste0("chr", human_CHR),
         gwas_hit_name = paste(human_CHR, gwas_pos, gwas_pos_end, paste0("PMID",PUBMEDID), str_replace_all(gwas_catalog_trait_name, "\\s+", "_"), sep = "__"))

# BED: chr, start, end, name
bed <- gwas_human %>% select(human_CHR, gwas_pos, gwas_pos_end, gwas_hit_name) %>%
  filter(!is.na(gwas_pos) & !is.na(human_CHR))

write.table(bed, out_humanGWAS_all_bed, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "	")

# ================================
# 11) Find overlaps: dog (lifted, 500kb) vs human GWAS Catalog hits (bedtools)
# ================================
system(paste("bedtools intersect -wao -a", out_dogClumps_onHuman_lifted_500kb_bed, 
             "-b", out_humanGWAS_all_bed, 
             ">", out_humanGWAS_dogClumps.bedtools))

# Parse bedtools output; keep the dog locus and the overlapping human hit name
intersect <- read.delim(out_humanGWAS_dogClumps.bedtools, header = FALSE,na.strings = c(".","","NA")) %>%
  rename(human_chr = V1, dog_start_on_human = V2, dog_end_on_human = V3, dog_locus = V4, gwas_hit_name = V8) %>%
  select(human_chr, dog_start_on_human, dog_end_on_human, dog_locus, gwas_hit_name) %>%
  as_tibble()

# Reload dog clump metadata written earlier (restore underscores in column names)
dog_clumps <- as_tibble(read.delim(out_dogClumps, header = TRUE, check.names = FALSE))
colnames(dog_clumps) <- gsub(" ", "_", colnames(dog_clumps))
dog_clumps <- dog_clumps %>%
  rename(dog_SNP = SNP) %>%
  select(dog_region,dog_locus, dog_SNP, dog_P, dog_trait)

dim(intersect)
# Expand overlap table with dog clump metadata (many-to-many expected)
intersect <- intersect %>% full_join(dog_clumps,relationship = "many-to-many") 
dim(intersect)

# Create a compact table of human GWAS hit metadata keyed by gwas_hit_name
human_clumps <- gwas_human %>%
  select(gwas_hit_name, gwas_catalog_trait_name, human_P, PUBMEDID, human_chr=human_CHR, gwas_pos, gwas_pos_end,SNPS, STUDY) %>%
  distinct()

# Attach human hit metadata to the overlap table
intersect <- intersect %>% left_join(human_clumps, relationship = "many-to-many") %>%
  distinct() 

# ================================
# 12) Build final merged overlap table and write to disk
# ================================
final <- all_dog_clumps %>% select(dog_locus,phenotype=dog_trait,region=dog_region,region_P=dog_P)

# Add lifted-over human coordinates for each dog locus and a boolean flag indicating lift success
final <- final %>% left_join(dog_liftOver_to_human %>% select(human_chr,dog_start_on_human=human_start,	dog_end_on_human=human_end,dog_locus),relationship = "many-to-many") %>% 
  mutate(lifted_over=if_else(is.na(human_chr),FALSE,TRUE)) %>% distinct()
final <- final %>% mutate(dog_start_on_human=as.integer(dog_start_on_human),dog_end_on_human=as.integer(dog_end_on_human))

# Add overlap hits and order by dog/human significance (note: distinct is referenced but not called here)
final <- final %>% left_join(intersect) %>%
  arrange(dog_P, human_P) %>% distinct

# Map GWAS Catalog trait names to canonical names used in this project
final <- final %>% left_join(human_trait_list %>% select(gwas_catalog_trait_name,human_trait_name=canonical_name) %>% distinct(),relationship = "many-to-many")

# Select and order output columns
final <- final %>% select(dog_locus,phenotype,region,dog_SNP,region_P,lifted_over,human_chr,dog_start_on_human,dog_end_on_human,gwas_hit_name,human_trait_name,gwas_catalog_trait_name,human_P,PUBMEDID,gwas_pos,gwas_pos_end) %>% distinct() %>% arrange(region_P,human_P)

# ================================
# 13) Cleanup and remove intermediate files
# ================================
rm(bed)
rm(raw_gwas_human)
rm(gwas_human)
rm(human_clumps)
rm(intersect)

write.table(final,"../../data/human_dog_GWAS_intersect/working.GWAS_RAW_OVERLAP.tsv",row.names=F,quote=F,sep="\t")

file.remove(out_dogClumps)
file.remove(out_dogClumps_bedFormat)
file.remove(out_dogClumps_onHuman_lifted_500kb_bed)
file.remove(out_humanGWAS_all_bed)
file.remove(out_humanGWAS_dogClumps.bedtools)
file.remove(out_dogClumps_onHuman_unmapped_bed)
file.remove(out_dogClumps_onHuman_lifted_bed)
