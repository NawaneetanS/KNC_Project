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

# 1. CHINA
china_mut_clin_sample <- fread("public_data/china_pancan_2020/data_clinical_sample.txt", nThread = 3)
china_mut_clin_patient <- fread("public_data/china_pancan_2020/data_clinical_patient.txt", nThread = 3)

china_mut_clin_sample <- china_mut_clin_sample %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(china_mut_clin_sample) <- NULL

china_mut_clin_patient <- china_mut_clin_patient %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(china_mut_clin_patient) <- NULL

china_mut_clin_merge <- merge.data.frame(china_mut_clin_sample, china_mut_clin_patient, by = "PATIENT_ID") %>%
  dplyr::rename(Tumor_Sample_Barcode = SAMPLE_ID)

china_luad <- china_mut_clin_merge %>% filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma")

china_clean <- china_luad %>%
  mutate(
    AGE = as.numeric(`DIAGNOSIS AGE`),
    STAGE = case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE == "0" ~ "0",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "I" ~ "I",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "II" ~ "II",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "III" ~ "III",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "IV" ~ "IV",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "I-II" ~ "I-II",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "III-IV" ~ "III-IV",
      TRUE ~ "Unknown"
    ),
    STAGE = factor(STAGE, levels = c("0", "I", "I-II", "II", "III", "III-IV", "IV", "Unknown")),
    SEX = factor(SEX, levels = c("Male", "Female")),
    SMOKE_STATUS = case_when(
      `SMOKE STATUS` == "Nonsmoker" ~ "Nonsmoker",
      `SMOKE STATUS` == "Smoker" ~ "Smoker",
      TRUE ~ "Unknown"
    ),
    SMOKE_STATUS = factor(SMOKE_STATUS, levels = c("Nonsmoker", "Smoker", "Unknown")),
    Treatment_Cleaned = case_when(
      TREATMENT == "Treatment-naive" ~ "Treatment Naive",
      grepl("Chemotherapy", TREATMENT) & !grepl("Targeted", TREATMENT) & !grepl("Radiation", TREATMENT) ~ "Chemotherapy Only",
      grepl("Targeted_Therapy", TREATMENT) & !grepl("Chemotherapy", TREATMENT) & !grepl("Radiation", TREATMENT) ~ "Targeted Therapy Only",
      TREATMENT == "Unknown" | TREATMENT == "" ~ "Unknown",
      TRUE ~ "Multi-modal / Other"
    ),
    Treatment_Cleaned = factor(Treatment_Cleaned, levels = c("Treatment Naive", "Chemotherapy Only", "Targeted Therapy Only", "Multi-modal / Other", "Unknown")),
    SAMPLE_TYPE = factor(SAMPLE_TYPE, levels = c("Primary", "Metastasis", "Recurrent")),
    SPECIMEN_TYPE = case_when(
      SPECIMEN_TYPE == "Biopsy/Paracentesis" ~ "Biopsy / Paracentesis",
      TRUE ~ SPECIMEN_TYPE
    ),
    SPECIMEN_TYPE = factor(SPECIMEN_TYPE, levels = c("Surgery", "Biopsy / Paracentesis")),
    TUMOR_PURTITY = as.numeric(TUMOR_PURTITY),
    TMB_NONSYNONYMOUS = as.numeric(TMB_NONSYNONYMOUS)
  )

# 2. SINGAPORE
sg_clin_patient <- fread("public_data/singapore_luad_2020/data_clinical_patient.txt", nThread = 3)
sg_clin_sample <- fread("public_data/singapore_luad_2020/data_clinical_sample.txt", nThread = 3)

sg_clin_patient <- sg_clin_patient %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(sg_clin_patient) <- NULL

sg_clin_sample <- sg_clin_sample %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(sg_clin_sample) <- NULL

sg_clin_merge <- merge.data.frame(sg_clin_sample, sg_clin_patient, by = "PATIENT_ID") %>%
  dplyr::rename(Tumor_Sample_Barcode = "SAMPLE_ID")

sg_clean <- sg_clin_merge %>%
  mutate(
    AGE = as.numeric(AGE),
    STAGE = factor(STAGE, levels = c("I", "II", "III", "IV")),
    SEX = factor(SEX, levels = c("Male", "Female")),
    HISTOLOGICAL_GRADE = factor(HISTOLOGICAL_GRADE, levels = c(
      "Well differentiated",
      "Well to Moderately differentiated",
      "Moderately differentiated",
      "Moderately to Poorly differentiated",
      "Poorly differentiated"
    )),
    Subtype = case_when(
      grepl("Acinar", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Acinar",
      grepl("Lepidic", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Lepidic",
      grepl("Papillary", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Papillary",
      grepl("Solid", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Solid",
      grepl("mucinous", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Mucinous",
      grepl("Micropapillary", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "Micropapillary",
      grepl("NOS", ADENOCARCINOMA_SUBTYPE_WHO2015, ignore.case = TRUE) ~ "NOS",
      is.na(ADENOCARCINOMA_SUBTYPE_WHO2015) | ADENOCARCINOMA_SUBTYPE_WHO2015 == "" ~ "Unknown",
      TRUE ~ "Other"
    ),
    PURITY = as.numeric(PURITY),
    TMB_NONSYNONYMOUS = as.numeric(TMB_NONSYNONYMOUS),
    SMOKING_STATUS = factor(SMOKING_STATUS, levels = c("No", "Yes")),
    CHEMOTHERAPY = factor(CHEMOTHERAPY, levels = c("No", "Yes")),
    TKI_TREATMENT = factor(TKI_TREATMENT, levels = c("No", "Yes"))
  )

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

# 4. MSK
msk_clin_patient <- fread("public_data/msk_impact_50k_2026/data_clinical_patient.txt", nThread = 3)
msk_clin_sample <- fread("public_data/msk_impact_50k_2026/data_clinical_sample.txt", nThread = 3)

msk_clin_patient <- msk_clin_patient %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(msk_clin_patient) <- NULL

msk_clin_sample <- msk_clin_sample %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()
rownames(msk_clin_sample) <- NULL

msk_clin_merge <- merge.data.frame(msk_clin_patient, msk_clin_sample, by = "PATIENT_ID") %>%
  dplyr::rename(Tumor_Sample_Barcode = SAMPLE_ID)

msk_luad <- msk_clin_merge %>% filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma")

msk_clean <- msk_luad %>%
  mutate(
    AGE = as.numeric(AGE_AT_DX),
    SEX = factor(SEX, levels = c("Male", "Female")),
    Sample_Type = factor(SAMPLE_TYPE, levels = c("Primary", "Metastasis", "Local Recurrence")),
    Met_Site = case_when(
      SAMPLE_TYPE == "Primary" ~ "Primary Tumor (N/A)",
      grepl("Lymph Node|LN", METASTATIC_SITE, ignore.case=TRUE) ~ "Lymph Node",
      grepl("Brain", METASTATIC_SITE, ignore.case=TRUE) ~ "Brain",
      grepl("Bone", METASTATIC_SITE, ignore.case=TRUE) ~ "Bone",
      grepl("Liver", METASTATIC_SITE, ignore.case=TRUE) ~ "Liver",
      grepl("Pleura", METASTATIC_SITE, ignore.case=TRUE) ~ "Pleura",
      grepl("Adrenal", METASTATIC_SITE, ignore.case=TRUE) ~ "Adrenal",
      grepl("Soft tissue|Soft Tissue", METASTATIC_SITE, ignore.case=TRUE) ~ "Soft Tissue",
      is.na(METASTATIC_SITE) | METASTATIC_SITE == "" | METASTATIC_SITE == "Not Applicable" | METASTATIC_SITE == "Not available" ~ "Unknown/Other",
      TRUE ~ "Other Metastasis"
    ),
    GENE_PANEL = factor(GENE_PANEL),
    MSI_SCORE = as.numeric(MSI_SCORE),
    MSI_TYPE = factor(MSI_TYPE, levels = c("Stable", "Indeterminate", "Instable", "Do not report")),
    TMB_SCORE = as.numeric(TMB_SCORE),
    TMB_SCORE = ifelse(TMB_SCORE < 0, NA_real_, TMB_SCORE),
    TUMOR_PURITY = as.numeric(TUMOR_PURITY),
    TUMOR_PURITY = ifelse(TUMOR_PURITY > 100, NA_real_, TUMOR_PURITY)
  )

# Save cleaned clinical dataframes
saveRDS(china_clean, "Tables/china_clean.rds")
saveRDS(sg_clean, "Tables/sg_clean.rds")
saveRDS(cptac_clean, "Tables/cptac_clean.rds")
saveRDS(msk_clean, "Tables/msk_clean.rds")

message("--- Loading CPTAC phosphoprotein dataset ---")

# Load phosphoprotein quantification data
cptac_phospho <- fread("public_data/luad_cptac_2020/data_phosphoprotein_quantification.txt", nThread = 8)

# Save the dataset to RDS
saveRDS(cptac_phospho, "Tables/cptac_phospho.rds")

message("Data preprocessing finished successfully!")
