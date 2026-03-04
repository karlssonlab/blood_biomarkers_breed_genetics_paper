# blood_biomarkers_breed_genetics_paper
Analysis and figure generation scripts for the manuscript
### Blood biomarkers and breed genetics of aging in pet dogs

#### The `analysis_scripts` directory contains scripts used for running the following analyses:

* Rscripts to run ANOVA on blood traits can be found in `anova`. 

* Rscript to run linear mixed effects regression modeling (LMER) on blood traits is in `lmer` subdirectory. The Rscript was run on a high-performance computing server (HPC) for each blood trait. An example bash script that submits the Rscript for each blood trait has the prefix lmer_efficient, and the wrapper bash script that submits the bash script for all blood traits has the prefix wrapper_lmer.

* Rscript `run.ANOVA_ON_LMER_SCORES.R` is used to run ANOVA on LMER breed ancestry scores to understand the effects of breed median lifespan and breed average weights on variation in breed ancestry scores of blood traits. 

* Rscript `cox_regression_survival_analysis_20251208.R` contains the code to run survival analysis on individual dogs using longitudinal data from the Dog Aging Project.

* Rscripts needed to assess the degree of overlap between dog and human blood trait genetic associations are available in `gwas_overlap`.

#### The `figure_generation_plotting_scripts` directory includes all Rscripts needed to create the main and supplemental figures of the manuscript.


#### All input data files (including supplemental data tables of the mansucript in R's RDS format) are available in DataDryad at DOI: 10.5061/dryad.sxksn03hq

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
