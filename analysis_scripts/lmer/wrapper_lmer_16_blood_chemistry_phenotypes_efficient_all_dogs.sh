#!/bin/bash

workingdir=/scratch/vsohrab/dap_lmer/blood_precision2024_lmer
scriptdir=${workingdir}
script=${scriptdir}/lmer_efficient_16_blood_chemistry_phenotypes.sh
phenotype_list=${workingdir}/bloodchem_16list.txt
ancestry_cutoff=1

num_runs=`wc -l ${phenotype_list} | awk '{print $1}'`

sbatch -p general -q public --array=1-${num_runs} -t 0-24:00:00 --mem=10G --mail-user=%u@asu.edu --mail-type=FAIL ${script} ${phenotype_list} ${ancestry_cutoff}
