#!/bin/bash

phenotype_list=$1
ancestry_cutoff=$2
lmerdir=/scratch/vsohrab/dap_lmer
workingdir=${lmerdir}/blood_precision2024_lmer
grm=/scratch/vsohrab/dap_genetic_set_2023/grm/DogAgingProject_2023_N-7627_canfam4_gp-0.70_biallelic_maf-0.01_geno-0.05_hwe-1e-20_midp_keep-fewhet
pheno_data=${workingdir}/input_files/precision2024_16bloodchem_phenotypes.tsv
randeffect_data=${workingdir}/input_files/age_at_sample_collection_chem_966_precision_dogs.tsv
ancestry=${lmerdir}/DAP_GeneticData_2023_ancestry_cleaned.tsv
lmer_script=${lmerdir}/linear_mixed_effects_regressions_age_randEffect_breed_analysis_forCommandline_with_M_ancestryThreshold_variable.R
outname_prefix=BreedLoad_BloodChem_Phenotype_Value
outdir=${workingdir}/individual_lmer_results_chem_all_dogs
phenotype=`sed -n ${SLURM_ARRAY_TASK_ID}p ${phenotype_list}` 

# load R 4.4.0
module load r-4.4.0-gcc-12.1.0

echo "started LMER analysis for ${phenotype} at $(date '+%Y-%m-%d %H:%M:%S')"

# submit lmer script
Rscript --vanilla --no-save ${lmer_script} ${grm} ${ancestry} ${pheno_data} ${randeffect_data} ${outdir} ${outname_prefix} ${phenotype} ${ancestry_cutoff}

echo "LMER analysis ended at $(date '+%Y-%m-%d %H:%M:%S')"
