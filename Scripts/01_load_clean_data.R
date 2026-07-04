#!/usr/bin/Rscript

# 01_load_clean_data.R
# Preprocesses raw clinical and mutation data, and saves cleaned objects as RDS/TSV files.

suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(janitor)
  library(data.table)
  library(tidyr)
  library(purrr)
})

setwd("/media/nannu1375/Backpack/Shankara/KNC")

# Create directories
dir.create("Tables", showWarnings = FALSE)
dir.create("Plots/Clinical", recursive = TRUE, showWarnings = FALSE)
dir.create("Plots/Mutational", recursive = TRUE, showWarnings = FALSE)
dir.create("Plots/Survival", recursive = TRUE, showWarnings = FALSE)

message("--- Loading and cleaning clinical datasets ---")

# 3. CPTAC
cptac_clin_patient <- fread("public_data/luad_cptac_2020/data_clinical_patient.txt", skip = "PATIENT_ID") %>% as.data.frame()
cptac_clin_sample <- fread("public_data/luad_cptac_2020/data_clinical_sample.txt", skip = "SAMPLE_ID") %>% as.data.frame()

# Load survival metadata from GDC to enable survival/Cox regression on the CPTAC cohort
cptac_surv <- fread("public_data/luad_cptac_gdc/data_clinical_patient.txt", skip = "PATIENT_ID") %>%
  as.data.frame() %>%
  dplyr::select(PATIENT_ID, OS_STATUS, OS_MONTHS)

cptac_clin_merge <- merge.data.frame(cptac_clin_patient, cptac_clin_sample, by = "PATIENT_ID") %>%
  dplyr::rename(Tumor_Sample_Barcode = SAMPLE_ID) %>%
  merge(y = cptac_surv, by = "PATIENT_ID", all.x = TRUE)

cptac_clean <- cptac_clin_merge %>%
  mutate(
    AGE = as.numeric(AGE),
    STAGE = case_when(
      grepl("^1|^I", STAGE, ignore.case=TRUE) ~ "I",
      grepl("^2|^II", STAGE, ignore.case=TRUE) ~ "II",
      grepl("^3|^III", STAGE, ignore.case=TRUE) ~ "III",
      grepl("^4|^IV", STAGE, ignore.case=TRUE) ~ "IV",
      TRUE ~ "Unknown"
    ),
    STAGE = factor(STAGE, levels = c("I", "II", "III", "IV", "Unknown")),
    T_Stage = factor("TX/Unknown", levels = c("T1", "T2", "T3", "T4", "TX/Unknown")),
    N_Stage = factor("NX/Unknown", levels = c("N0", "N1", "N2", "N3", "NX/Unknown")),
    M_Stage = factor("MX/Unknown", levels = c("M0", "M1", "MX/Unknown")),
    SEX = factor(SEX, levels = c("Male", "Female")),
    Subtype = case_when(
      grepl("acinar", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Acinar",
      grepl("papillary", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Papillary",
      grepl("solid", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Solid",
      grepl("mucinous|colloid", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Mucinous",
      grepl("lepidic", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Lepidic / BAC",
      grepl("mixed", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "Mixed",
      grepl("nos", DOMINANT_HISTOLOGICAL_SUBTYPE, ignore.case = TRUE) ~ "NOS",
      is.na(DOMINANT_HISTOLOGICAL_SUBTYPE) | DOMINANT_HISTOLOGICAL_SUBTYPE == "" | DOMINANT_HISTOLOGICAL_SUBTYPE == "DOMINANT_HISTOLOGICAL_SUBTYPE" ~ "Unknown",
      TRUE ~ "Other"
    ),
    SMOKING_STATUS = case_when(
      grepl("non", SMOKING_STATUS, ignore.case = TRUE) ~ "No",
      grepl("smoker", SMOKING_STATUS, ignore.case = TRUE) ~ "Yes",
      TRUE ~ "Unknown"
    ),
    SMOKING_STATUS = factor(SMOKING_STATUS, levels = c("No", "Yes", "Unknown")),
    PRIOR_DX = factor("Unknown", levels = c("No", "Yes", "Unknown")),
    Purity = as.numeric(TUMOR_PURITY_BYESTIMATE_RNASEQ),
    TMB_NONSYNONYMOUS = as.numeric(TMB_NONSYNONYMOUS)
  )

# Save cleaned clinical dataframes
saveRDS(cptac_clean, "Tables/cptac_clean.rds")

message("--- Loading CPTAC phosphoprotein dataset ---")

# Load phosphoprotein quantification data
cptac_phospho <- fread("public_data/luad_cptac_2020/data_phosphoprotein_quantification.txt", nThread = 8)

# Save the dataset to RDS
saveRDS(cptac_phospho, "Tables/cptac_phospho.rds")

message("Data preprocessing finished successfully!")
