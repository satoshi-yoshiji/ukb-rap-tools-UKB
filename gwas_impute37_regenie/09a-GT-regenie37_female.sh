#!/bin/sh

# This script runs the step1 analysis for regenie.

# Requirements: 
# (0-4). Please see readme.md
# 5. Must have executed scripts 01, 12, and 03 from 01_prep_gtfile_4_GWAS
# - liftover not necessary as we are staying in grch37

# How to Run:
# Run this shell script using: 
# sh 09a-GT-regenie37.sh
# on the command line on your machine

# Inputs: 
# - /data/gt_genrel_block/ukb22418_allQC_v2_mrg_prun_cohort.bed 
# - /data/gt_genrel_block/ukb22418_allQC_v2_mrg_prun_cohort.bim
# - /data/gt_genrel_block/ukb22418_allQC_v2_mrg_prun_cohort.fam
# - /gwas_cohort_textfiles/phenotypes.v08-04-22.txt
# - /gwas_cohort_textfiles/covariates.v08-04-22.txt

# Output:
# - /data/ap_imp37_gwas/ap37_results_1.loco.gz - Leave One Chromosome Out results 
# - /data/ap_imp37_gwas/ap37_results_pred.list - List of files generated this step 
# - /data/ap_1im37_gwas/ap37_results.log

# Steps:
# 1. run regenie


# Variables
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"
#set this to the exome data field for your release
project="UKB"
data_field="ukb22418"
# working dir and reference text file dir
data_file_dir="/03.sex_stratified_BMI/gt_genrel_block/"
txt_file_dir="/03.sex_stratified_BMI/phenotype_covariate/"
gwas_file_dir="/03.sex_stratified_BMI/ap_imp37_gwas_female/"
pheno_file="female_phenotype.tsv"
covar_file="female_covariate.tsv"

run_regenie_cmd="regenie --step 1 --out ap37_results \
 --bed ${data_field}_allQC_v2_mrg_prun_cohort \
 --phenoFile ${pheno_file} \
 --covarFile ${covar_file} \
 --phenoCol BMI \
 --bsize 1000 --qt --loocv --gz --threads 16 "

# removed options
# --phenoCol AP --covarCol age \
# --covarCol sex.b2 --covarCol bmi \
# --covarCol smoke --covarCol pca{1:9} \

dx run swiss-army-knife -iin="${data_file_dir}/${data_field}_allQC_v2_mrg_prun_cohort.bed"\
   -iin="${data_file_dir}/${data_field}_allQC_v2_mrg_prun_cohort.bim"\
   -iin="${data_file_dir}/${data_field}_allQC_v2_mrg_prun_cohort.fam"\
   -iin="${txt_file_dir}/${pheno_file}" \
   -iin="${txt_file_dir}/${covar_file}" \
   -icmd="${run_regenie_cmd}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16" \
   --destination="${project}:${gwas_file_dir}" --brief --yes;


