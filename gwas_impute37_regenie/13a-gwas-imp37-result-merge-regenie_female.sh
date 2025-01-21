#!/bin/sh

# This script combines all the regenie result files into a single cohort file.

# Requirements: 
# 0-4 - please refer to readme.md
# 5. Must have executed: 
# - all scripts including liftover from 01_prep_gtfile_4_GWAS
# - 09b-step1-regenie.sh
# - 11b-gwas-s2-imp37-qc-filter.sh
# - 12b-gwas-s2-imp37-regenie.sh

# How to Run:
# Run this shell script using: 
#   sh 13b-gwas-imp37-result-merge-regenie.sh
# on the command line on your own machine

# Inputs (for each chromosome):
# Note that you can adjust the output directory by setting the data_file_dir variable
# - /Data/ap_imp37_gwas/assoc.c1_BMI.regenie.gz - regenie results for chromosome 1 
# - /Data/ap_imp37_gwas/assoc.c2_BMI.regenie.gz - regenie results for chromosome 2 
# - /Data/ap_imp37_gwas/assoc.c3_BMI.regenie.gz - regenie results for chromosome 3 
# - /Data/ap_imp37_gwas/assoc.c4_BMI.regenie.gz - regenie results for chromosome 4 
# - etc.

# Output :
# - /data/ap_wes_gwas/assoc.regenie.merged.txt - merged results for all chromosomes in tab-delimited format

# Steps:
# 1. Use dxFUSE to copy the regenie files into the container storage
# 2. unzip regenie files
# 3. add the header back to the top of the merged file with echo command
# 4.  for loop  each .regenie file
#       remove header with tail
#       transform to tab delimited with tr
#       save it into $out_file
# 5. delete regenie files

#data_file_dir="/data/ap_imp37_gwas"
project="UKB"
#data_field="ukb22418"
data_field="ukb22828" # NOTE THE ID IS DIFFERENT FOR THE IMPUTATION DATA
# working dir and reference text file dir
data_file_dir="/03.sex_stratified_BMI/ap_imp37_gwas_female/"
txt_file_dir="/03.sex_stratified_BMI/phenotype_covariate/"
gwas_file_dir="/03.sex_stratified_BMI/ap_imp37_gwas_female/"
pheno_file="female_phenotype.tsv"
covar_file="female_covariate.tsv"

merge_cmd='out_file="assoc.regenie.merged.txt"
     gunzip *.regenie.gz

echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file

files="./*.regenie"
for f in $files
do
   tail -n+2 $f | tr " " "\t" >> $out_file
done

rm *.regenie ' 


dx run swiss-army-knife \
    -iin="/${gwas_file_dir}assoc.c1_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c2_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c3_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c4_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c5_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c6_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c7_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c8_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c9_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c10_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c11_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c12_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c13_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c14_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c15_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c16_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c17_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c18_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c19_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c20_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c21_BMI.regenie.gz" \
    -iin="/${gwas_file_dir}assoc.c22_BMI.regenie.gz" \
   -icmd="${merge_cmd}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:${gwas_file_dir}" --brief --yes 
