#!/bin/sh

# Requirements: 
# Please refer to readme.md for more information about the requirements

# How to Run:
# Run this script using: 
# sh 01-GTprep-merge-files-dxfuse.sh 
# on the command line

# Inputs
# - /gwas_cohort_textfiles/phenotypes.v08-04-22.txt - not strictly needed, but swiss-army-knife needs at least one input

# Outputs
# - ukb22418_all_v2_merged.bed  
# - ukb22418_all_v2_merged.bim 
# - ukb22418_all_v2_merged.fam 

# Steps:
# 1. copy plink files to instance 1-22+X
# 2. create a merge list file
# 3. merge via plink2
# 4. delete temporary files


# Variables
# set this to the exome path and data field for your release
project="UKB"
data_field="ukb22418"
# working dir and reference text file dir
data_file_dir="/03.sex_stratified_BMI/gt_genrel_block/"
txt_file_dir="/03.sex_stratified_BMI/"

#cmd to run (use as input with `-icmd="{$run_merge}"`)

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ;\
        cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_cX_* . ;\
        ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt; \
        plink --merge-list files_to_merge.txt --make-bed\
        --out ukb22418_all_v2_merged;\
        rm files_to_merge.txt;\
        rm ukb22418_c*"

dx run swiss-army-knife -iin="${txt_file_dir}/GWAS_all_covariates_sex_stratified_BMI_participant.tsv" \
   -icmd="${run_merge}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/03.sex_stratified_BMI/" --brief --yes 
