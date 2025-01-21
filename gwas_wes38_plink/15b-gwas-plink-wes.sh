#!/bin/bash

# This script runs the logistic regression GWAS using plink2 on WES data.

# Requirements: 
# 0-4 - please refer to readme.md
# 5. Must have executed: 
# - 11b-gwas-s2-wes38-qc-filter.sh

# How to Run:
# Run this shell script using: 
#   sh 15b-gwas-plink-wes.sh  
# on the command line on your own machine

# Inputs:
# Note that you can adjust the output directory by setting the data_file_dir variable
# - /gwas_cohort_textfiles/phenotypes.v08-01-22.txt -  (outside of scope)
# - /gwas_cohort_textfiles/covariates.v08-01-22.txt -  (outside of scope)

# Additional inputs
# for each chromosome, you will run a separate worker
# - /{data_file_dir}/ukb23158_c1_v3.bed - from 11a
# - /{data_file_dir}/ukb23158_c1_v3.bim 
# - /{data_file_dir}/ukb23158_c1_v3.fam 

# Outputs (for each chromosome):
# - /data/ap_wes_gwas/plink/ukb23158_AP_c1_v1.AP.glm.logistic.hybrid - plink results for chromosome 1 
# - /data/ap_wes_gwas/plink/ukb23158_AP_c1_v1.log  - plink log for chromosome 1

# Steps:
# 1. for each chromosome 1-22 and X:
#       - run logistic regression using plink 


exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"
#set this to the exome data field for your release
data_field="ukb23158"
data_file_dir="/data/ap_wes_gwas/"
txt_file_dir="/gwas_cohort_textfiles/"


for i in {1..22}; do

  run_plink_imp="plink2 --bfile ${data_field}_c${i}_b0_v1 --1 \
      --extract WES_c${i}_snps_qc_pass.snplist \
      --pheno phenotypes.v08-01-22.txt --pheno-name AP \
      --covar covariates.v08-01-22.txt --covar-name age,bmi,smoke,pca1-pca6 \
      --logistic sex hide-covar --out ${data_field}_AP_c${i}_v1"


  dx run swiss-army-knife -iin="${exome_file_dir}/${data_field}_c${i}_b0_v1.bed" \
   -iin="${exome_file_dir}/${data_field}_c${i}_b0_v1.bim" \
   -iin="${exome_file_dir}/${data_field}_c${i}_b0_v1.fam" \
   -iin="${data_file_dir}/WES_c${i}_snps_qc_pass.snplist" \
   -iin="${txt_file_dir}/phenotypes.v08-01-22.txt" \
   -iin="${txt_file_dir}/covariates.v08-01-22.txt" \
   -icmd="${run_plink_imp}" --tag="plink" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/data/ap_wes_gwas/plink/" --brief --yes

done

# now run chrX

  run_plink_imp="plink2 --bfile "${data_field}_cX_b0_v1" --1 \
      --extract WES_cX_snps_qc_pass.snplist \
      --pheno phenotypes.v08-01-22.txt --pheno-name AP \
      --covar covariates.v08-01-22.txt --covar-name age,bmi,smoke,pca1-pca6 \
      --logistic sex hide-covar --out ${data_field}_AP_cX_v1"

  dx run swiss-army-knife -iin="${exome_file_dir}/${data_field}_cX_b0_v1.bed" \
   -iin="${exome_file_dir}/${data_field}_cX_b0_v1.bim" \
   -iin="${exome_file_dir}/${data_field}_cX_b0_v1.fam" \
   -iin="${data_file_dir}/WES_cX_snps_qc_pass.snplist" \
   -iin="${txt_file_dir}/phenotypes.v08-01-22.txt" \
   -iin="${txt_file_dir}/covariates.v08-01-22.txt" \
   -icmd="${run_plink_imp}" --tag="plink" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/data/ap_wes_gwas/plink/" --brief --yes

