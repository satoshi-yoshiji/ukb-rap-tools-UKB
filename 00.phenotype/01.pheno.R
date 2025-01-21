library(tidyverse)
library(vroom)
library(magrittr)
library(readxl)
library(RNOmni)

setwd("/scratch/richards/satoshi.yoshiji/40.UKB_TG_to_HDL/01.phenotype")

# Read the clinical data
pheno0 <- vroom('clinical_data.tsv')

# Step 2: Create a named vector for renaming PC columns
PC_names <- setNames(paste0("p22009_a", 1:20), paste0("PC", 1:20))

# Inspect PC_names to ensure correct mapping
print(PC_names)

# Step 3: Rename the PC columns
pheno_renamed <- pheno0 %>%
  rename(!!!PC_names)

# Verify the renaming
print(colnames(pheno_renamed))
head(pheno_renamed$PC1)  # Should display numeric values from p22009_a1

# Step 4: Select and rename other columns
pheno <- pheno_renamed %>%
  transmute(
    FID = eid,
    IID = eid,
    Sex = p31,
    Genetic_sex = p22001,
    Age_at_recruitment = p21022,
    UK_Biobank_assessment_centre = p54_i0,
    Sex_chromosome_aneuploidy = p22019,
    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
    PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20,
    TG = p30870_i0,
    HDL = p30760_i0
  )

# Step 5: Verify the final dataframe
print(colnames(pheno))
print(nrow(pheno))  # Should output 502239
head(pheno)

# transform Sex and Genetic_sex columns to binary (0/1)
#pheno %<>% mutate(Sex = ifelse(Sex == 'Female', 0, 1)) %>% mutate(Genetic_sex = ifelse(Genetic_sex == 'Female', 0, 1)) # Female = 0 , Male = 1
table(pheno$Sex) # 0 = Female

# one hot encoding for UK_Biobank_assessment_centre
# Step 1: # One-hot encoding using pivot_wider
pheno_encoded <- pheno %>%
  mutate(Flag = 1) %>%
  pivot_wider(names_from = UK_Biobank_assessment_centre, values_from = Flag, values_fill = list(Flag = 0))

# Step 2: Rename columns
new_colnames <- c("centre" %>% paste0(seq_len(length(unique(pheno$UK_Biobank_assessment_centre)))))
names(pheno_encoded)[which(names(pheno_encoded) %in% unique(pheno$UK_Biobank_assessment_centre))] <- new_colnames

# # Genotype_measurement_batch
# #  transform the Genotype_measurement_batch column such that"Batch_" = 0 and "UKBiLEVEAX_" = 1,
# pheno_encoded %<>%
#   mutate(Genotype_measurement_batch = case_when(
#     startsWith(Genotype_measurement_batch, "Batch_") ~ 0,
#     startsWith(Genotype_measurement_batch, "UKBiLEVEAX_") ~ 1,
#     TRUE ~ NA_integer_
#   ))


# ancestry
load('/scratch/richards/tianyuan.lu2/UKB_common/Ancestry_cluster.RData')
cluster %<>% mutate(IID = as.numeric(rownames(cluster)), FID = as.numeric(rownames(cluster)))

pheno2 <- inner_join(pheno_encoded, cluster %>% select(IID, CLUST, Ancestry), by = 'IID')

# Genotype measurement batch
geno_batch <- vroom('data_participant.tsv', delim = '\t')
colnames(geno_batch) <- c('IID', 'genotype_batch')
pheno3 <- inner_join(pheno2, geno_batch, by = c('IID'))
pheno3 %<>%
  mutate(genotype_batch = case_when(
    startsWith(genotype_batch, "Batch_") ~ 0,
    startsWith(genotype_batch, "UKBiLEVEAX_") ~ 1,
    TRUE ~ NA_integer_
  ))
table(pheno3$genotype_batch)

# QC pheno
pheno_qc <- pheno3 %>% 
  filter(Sex == Genetic_sex | is.na(Genetic_sex)) %>% select(-Genetic_sex) %>%  # sex mismatch
  filter(is.na(Sex_chromosome_aneuploidy)) %>% select(-Sex_chromosome_aneuploidy)

# check  
nrow(pheno3)
nrow(pheno_qc)  

# ancestry
table(pheno3$Ancestry)

# remove na
pheno4 <- pheno_qc %>% drop_na()
nrow(pheno4)
table(pheno4$Ancestry)

pheno4 %<>% mutate(Ancestry_merged = case_when(Ancestry == 'White British' ~ 'EUR',
                                              Ancestry == 'Admixed White' ~ 'EUR',
                                              Ancestry == 'African' ~ 'AFR',
                                              Ancestry == 'East Asian' ~ 'EAS',
                                              Ancestry == 'South Asian' ~ 'SAS',
                                              TRUE ~ NA
                                              ))

table(pheno4[c('Sex', 'Ancestry')])
table(pheno4[c('CLUST', 'Ancestry')])
table(pheno4[c('Sex', 'Ancestry_merged')])

# finally, create all required covariates
# For the discovery cohort, association models included the following covariates: age, age2, sex,
# age*sex, age2*sex, batch, UKB centre, UKB genetic array, time between blood sampling and
# measurement and the first 20 genetic principal components (PCs).
pheno4  %<>% mutate(Age2 = Age_at_recruitment*Age_at_recruitment) %>% 
              mutate(TG_to_HDL = TG/HDL)

hist(pheno4$TG*88.57)
max(pheno4$TG*88.57)
hist(pheno4$HDL*38.67)
max(pheno4$HDL*38.67)
hist(pheno4$TG_to_HDL)
max(pheno4$TG_to_HDL)

# relocate columns
pheno4 %<>% relocate(FID, IID,
                      Age_at_recruitment,
                      Age2,
                      Sex,
                      Ancestry)

# remmove duplicated individuals
pheno5 <- pheno4 %>% filter(!duplicated(IID))
nrow(pheno5) # 423557

table(pheno5$Ancestry_merged)
table(pheno5[c('Sex', 'Ancestry_merged')])
table(pheno5[c('CLUST', 'Ancestry_merged')])
table(pheno5[c('Sex', 'Ancestry_merged')])

#check final colnames
colnames(pheno5)

# perform RNT on the final samples
pheno6 <- pheno5 %>% mutate(RNT_TG_to_HDL = RankNorm(TG_to_HDL))

hist(pheno6$RNT_TG_to_HDL)
sd(pheno6$RNT_TG_to_HDL)

# reorder and curate
pheno7 <- pheno6 %>% 
  select(FID, IID, Age_at_recruitment, Age2, Sex, RNT_TG_to_HDL, Ancestry_merged, genotype_batch, matches("centre"), matches("PC")) %>%
  transmute(
    FID, 
    IID,
    Age = Age_at_recruitment,
    Age2,
    Sex,
    RNT_TG_to_HDL,
    Ancestry_merged,
    genotype_batch,
    # All columns matching 'centre' and 'PC' are automatically included from select()
    across(matches("centre")),
    across(matches("PC"))
  )

# remove opted-out individuals
optout <- vroom('/scratch/richards/satoshi.yoshiji/08.UKB/excluded_individuals/w27449_2023-04-25.tsv', col_names = F)
`%!in%` <- Negate(`%in%`)
pheno7 %>% nrow() 
pheno7 %<>% filter(IID %!in% optout$X1)
nrow(pheno7)

# save
# write_tsv(pheno7 %>% filter(Ancestry_merged == 'EUR') %>% select(FID, IID, RNT_TG_to_HDL), file = 'output/EUR_both_pheno.tsv')
# write_tsv(pheno7 %>% filter(Ancestry_merged == 'EUR') %>% select(-RNT_TG_to_HDL, -Ancestry_merged), file = 'output/EUR_both_covar.tsv')

# Define the function to write the files for a specific ancestry and sex
write_pheno_and_covar <- function(data, ancestry, sex_label = "both", output_dir = "output") {
  
  # Create file names based on ancestry and sex
  pheno_file <- paste0(output_dir, "/", ancestry, "_", sex_label, "_pheno.tsv")
  covar_file <- paste0(output_dir, "/", ancestry, "_", sex_label, "_covar.tsv")
  id_file <- paste0(output_dir, "/", ancestry, "_", sex_label, "_FIDIID_withheader.tsv")
  id_noheader_file <- paste0(output_dir, "/", ancestry, "_", sex_label, "_FIDIID_noheader.tsv")
  
  # Filter the data based on ancestry and sex
  filtered_data <- data %>% 
    filter(Ancestry_merged == ancestry) %>% 
    {if(sex_label != "both") filter(., Sex == ifelse(sex_label == "female", 0, 1)) else .}
  
  # Write the phenotype file
  write_tsv(filtered_data %>% select(FID, IID, RNT_TG_to_HDL), file = pheno_file)
  
  # Write the covariate file (excluding RNT_TG_to_HDL and Ancestry_merged)
  write_tsv(filtered_data %>% select(-RNT_TG_to_HDL, -Ancestry_merged), file = covar_file)
  
  # Write the covariate file (excluding RNT_TG_to_HDL and Ancestry_merged)
  write_tsv(filtered_data %>% select(FID, IID), file = id_file)
  
  # Write the covariate file (excluding RNT_TG_to_HDL and Ancestry_merged)
  write_tsv(filtered_data %>% select(FID, IID), file = id_noheader_file, col_names = F)
}

# Apply the function for each ancestry and sex
ancestries <- c("EUR", "AFR", "EAS", "SAS")
sex_labels <- c("both", "female", "male")

for (ancestry in ancestries) {
  for (sex_label in sex_labels) {
    write_pheno_and_covar(pheno7, ancestry, sex_label)
  }
}


