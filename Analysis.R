#!usr/bin/Rscript

# Load necessary libraries for data manipulation, genomics analysis, and visualization
suppressPackageStartupMessages({
  library(dplyr)        # Data manipulation
  library(maftools)     # Analysis of Mutation Annotation Format (MAF) files
  library(biomaRt)      # Interface to BioMart databases
  library(janitor)      # Data cleaning and rounding
  library(survminer)    # Survival analysis visualization
  library(survival)     # Survival analysis
  library(data.table)   # Fast data reading and manipulation
  library(purrr)        # Functional programming tools
  library(tidyr)        # Data tidying
  library(ggplot2)      # Data visualization
  library(scales)       # Scale functions for visualization
  library(ggpubr)       # Publication-ready plots
  library(patchwork)    # Combining multiple plots
  library(ggsci)        # Scientific journal color palettes
  library(tidyverse)
})

# Set working directory to the project root
setwd("/mnt/Linux_storage/KNC")

# Create output directories if they do not exist
dir.create("Plots/Clinical", recursive = TRUE, showWarnings = FALSE)
dir.create("Plots/Mutational", recursive = TRUE, showWarnings = FALSE)
dir.create("Plots/Survival", recursive = TRUE, showWarnings = FALSE)


## China mutation clinical data - Reading, cleaning and merging 2 clinical files.

# Load clinical sample and patient data for the China Pan-Cancer cohort
china_mut_clin_sample <- fread("public_data/china_pancan_2020/data_clinical_sample.txt",
                    nThread = 3)

china_mut_clin_patient <- fread("public_data/china_pancan_2020/data_clinical_patient.txt",
                               nThread = 3)

# Clean sample data: remove metadata rows, set proper column names, and convert to data frame
china_mut_clin_sample <- china_mut_clin_sample %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(china_mut_clin_sample) <- NULL

# Clean patient data: remove metadata rows, set proper column names, and convert to data frame
china_mut_clin_patient <- china_mut_clin_patient %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(china_mut_clin_patient) <- NULL

# Merge sample and patient clinical data on PATIENT_ID
china_mut_clin_merge <- merge.data.frame(china_mut_clin_sample, china_mut_clin_patient, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

# Rename SAMPLE_ID to Tumor_Sample_Barcode for compatibility with maftools
china_mut_clin_merge <- china_mut_clin_merge %>% 
  rename(Tumor_Sample_Barcode = SAMPLE_ID)

## Singapore clinical data

# Load clinical data for the Singapore LUAD cohort
sg_clin_patient <- fread("public_data/singapore_luad_2020/data_clinical_patient.txt",
                 nThread = 3)

sg_clin_sample <- fread("public_data/singapore_luad_2020/data_clinical_sample.txt",
                        nThread = 3)

# Clean patient clinical data
sg_clin_patient <- sg_clin_patient %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(sg_clin_patient) <- NULL

# Clean sample clinical data
sg_clin_sample <- sg_clin_sample %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(sg_clin_sample) <- NULL

# Merge sample and patient data for Singapore cohort
sg_clin_merge <- merge.data.frame(sg_clin_sample, sg_clin_patient, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

# Standardize sample column name
sg_clin_merge <- sg_clin_merge %>% 
  rename(Tumor_Sample_Barcode = SAMPLE_ID)

## MSK clinical data

# Load clinical data for the MSK Impact cohort
msk_clin_patient <- fread("public_data/msk_impact_50k_2026/data_clinical_patient.txt",
                  nThread = 3)

msk_clin_sample <- fread("public_data/msk_impact_50k_2026/data_clinical_sample.txt",
                         nThread = 3)

# Clean patient clinical data
msk_clin_patient <- msk_clin_patient %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(msk_clin_patient) <- NULL

# Clean sample clinical data
msk_clin_sample <- msk_clin_sample %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(msk_clin_sample) <- NULL

# Merge patient and sample data for MSK cohort
msk_clin_merge <- merge.data.frame(msk_clin_patient, msk_clin_sample, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

# Standardize sample column name
msk_clin_merge <- msk_clin_merge %>% 
  rename(Tumor_Sample_Barcode = SAMPLE_ID)

## TCGA LUAD clinical data

# Load clinical data for the MSK Impact cohort
tcga_clin_patient <- fread("public_data/TCGA/luad_tcga_gdc/data_clinical_patient.txt",
                          nThread = 3)

tcga_clin_sample <- fread("public_data/TCGA/luad_tcga_gdc/data_clinical_sample.txt",
                         nThread = 3)

# Clean patient clinical data
tcga_clin_patient <- tcga_clin_patient %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(tcga_clin_patient) <- NULL

# Clean sample clinical data
tcga_clin_sample <- tcga_clin_sample %>% 
  slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame()

rownames(tcga_clin_sample) <- NULL

# Merge patient and sample data for MSK cohort
tcga_clin_merge <- merge.data.frame(tcga_clin_patient, tcga_clin_sample, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

# Standardize sample column name
tcga_clin_merge <- tcga_clin_merge %>% 
  rename(Tumor_Sample_Barcode = SAMPLE_ID)

# Remove intermediate clinical data objects to free up memory
rm(china_mut_clin_sample,
   china_mut_clin_patient, 
   china_rna_clin_sample, 
   china_rna_clin_patient, 
   sg_clin_sample,
   sg_clin_patient,
   msk_clin_sample,
   msk_clin_patient,
   tcga_clin_patient,
   tcga_clin_sample
)

## ==============================================
#            Clinical dashboard
# ==============================================

# Define a generic clinical dashboard function
create_clinical_dashboard <- function(df, cohort_name, plot_vars, outfile) {
  plot_list <- list()
  
  for (var_info in plot_vars) {
    col_name <- var_info$col
    var_label <- var_info$label
    var_type <- var_info$type
    
    # Filter out NAs for this plot
    plot_df <- df %>% 
      filter(!is.na(.data[[col_name]]), as.character(.data[[col_name]]) != "", as.character(.data[[col_name]]) != "NA")
    
    # Skip if no data
    if (nrow(plot_df) == 0) next
    
    if (var_type == "continuous") {
      # Make sure it is numeric
      plot_df[[col_name]] <- as.numeric(plot_df[[col_name]])
      plot_df <- plot_df %>% filter(!is.na(.data[[col_name]]))
      
      fill_col <- if (!is.null(var_info$fill)) var_info$fill else "#2C7FB8"
      
      if (col_name == "AGE") {
        # Keep density curve for age
        p <- ggplot(plot_df, aes(x = .data[[col_name]])) +
          geom_histogram(aes(y = after_stat(density)), fill = fill_col, color = "white", alpha = 0.6, bins = 15) +
          geom_density(color = "#1D3557", linewidth = 1.2, fill = fill_col, alpha = 0.1) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = var_label, y = "Density")
      } else {
        # Standard frequency histogram without density line
        p <- ggplot(plot_df, aes(x = .data[[col_name]])) +
          geom_histogram(fill = fill_col, color = "white", alpha = 0.85, bins = 15) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = var_label, y = "Patient Count")
      }
      
      p <- p + theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        axis.line = element_line(color = "grey80"),
        plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
      )
    } else if (var_type == "categorical") {
      # If not already a factor, convert to factor and order by frequency
      if (!is.factor(plot_df[[col_name]])) {
        plot_df[[col_name]] <- factor(plot_df[[col_name]])
        plot_df[[col_name]] <- reorder(plot_df[[col_name]], plot_df[[col_name]], FUN = length)
      }
      
      cat_counts <- plot_df %>% count(.data[[col_name]], name = "n")
      num_cats <- nrow(cat_counts)
      
      if (num_cats > 6) {
        # Horizontal bar plot for many categories to prevent label overlapping
        p <- ggplot(plot_df, aes(y = .data[[col_name]], fill = .data[[col_name]])) +
          geom_bar(width = 0.7, alpha = 0.85, color = "white", linewidth = 0.3) +
          geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.2, fontface = "bold", size = 3) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
          scale_y_discrete(labels = scales::label_wrap(15)) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = "Patient Count", y = NULL) +
          theme(
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey95"),
            axis.line = element_line(color = "grey80"),
            plot.margin = margin(t = 15, r = 15, b = 15, l = 10)
          )
      } else {
        # Vertical bar plot for few categories
        p <- ggplot(plot_df, aes(x = .data[[col_name]], fill = .data[[col_name]])) +
          geom_bar(width = 0.6, alpha = 0.85, color = "white", linewidth = 0.3) +
          geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, fontface = "bold", size = 3) +
          scale_x_discrete(labels = scales::label_wrap(12)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = NULL, y = "Patient Count") +
          theme(
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey95"),
            axis.line = element_line(color = "grey80"),
            plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
          )
      }
      
      # Apply custom palette if provided
      if (!is.null(var_info$palette)) {
        p <- p + scale_fill_manual(values = var_info$palette, breaks = names(var_info$palette))
      } else {
        if (num_cats <= 4) {
          p <- p + scale_fill_manual(values = c("#7FCDBB", "#41B6C4", "#1D91C0", "#081D58", "#E2E8F0"))
        } else {
          p <- p + scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8)
        }
      }
    }
    
    plot_list[[col_name]] <- p
  }
  
  num_plots <- length(plot_list)
  num_cols <- min(3, num_plots)
  num_rows <- ceiling(num_plots / num_cols)
  
  # Combine using patchwork
  dashboard <- wrap_plots(plot_list, ncol = num_cols) +
    plot_annotation(
      title = paste(cohort_name, "LUAD Clinical & Pathological Dashboard"),
      subtitle = paste0("Total Cohort Size: N = ", nrow(df), " Patients"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 12, b = 4)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#4A5568", margin = margin(b = 12))
      )
    )
  
  png(filename = outfile, width = 4.5 * num_cols, height = 3.5 * (num_rows + 0.2), units = "in", res = 600)
  print(dashboard)
  dev.off()
}

# Subset each cohort to Lung Adenocarcinoma
china_luad <- china_mut_clin_merge %>% filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma")
sg_luad    <- sg_clin_merge # Singapore is all LUAD
tcga_luad  <- tcga_clin_merge # TCGA is all LUAD
msk_luad   <- msk_clin_merge %>% filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma")

# 1. SINGAPORE CLINICAL & PATHOLOGICAL DASHBOARD
sg_clean <- sg_luad %>%
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

sg_vars <- list(
  list(col = "AGE", type = "continuous", label = "Age (Years)"),
  list(col = "SEX", type = "categorical", label = "Gender", palette = c("Male" = "#2C7FB8", "Female" = "#D53F8C")),
  list(col = "STAGE", type = "categorical", label = "Pathological Stage", palette = c("I" = "#7FCDBB", "II" = "#41B6C4", "III" = "#1D91C0", "IV" = "#081D58")),
  list(col = "HISTOLOGICAL_GRADE", type = "categorical", label = "Histological Grade"),
  list(col = "SMOKING_STATUS", type = "categorical", label = "Smoking Status", palette = c("No" = "#319795", "Yes" = "#ED8936")),
  list(col = "Subtype", type = "categorical", label = "Adenocarcinoma Subtype"),
  list(col = "CHEMOTHERAPY", type = "categorical", label = "Received Chemotherapy", palette = c("No" = "#A0AEC0", "Yes" = "#319795")),
  list(col = "TKI_TREATMENT", type = "categorical", label = "Received TKI Treatment", palette = c("No" = "#A0AEC0", "Yes" = "#805AD5")),
  list(col = "PURITY", type = "continuous", label = "Tumor Purity", fill = "#319795"),
  list(col = "TMB_NONSYNONYMOUS", type = "continuous", label = "TMB (Nonsynonymous)", fill = "#D95F02")
)

create_clinical_dashboard(sg_clean, "Singapore", sg_vars, "Plots/Clinical/Singapore_clinical_dashboard.png")

# 2. TCGA CLINICAL & PATHOLOGICAL DASHBOARD
tcga_clean <- tcga_luad %>%
  mutate(
    AGE = as.numeric(AGE),
    STAGE = case_when(
      grepl("^Stage I$|^Stage IA$|^Stage IB$", PATH_STAGE, ignore.case=TRUE) ~ "I",
      grepl("^Stage II$|^Stage IIA$|^Stage IIB$", PATH_STAGE, ignore.case=TRUE) ~ "II",
      grepl("^Stage III$|^Stage IIIA$|^Stage IIIB$", PATH_STAGE, ignore.case=TRUE) ~ "III",
      grepl("^Stage IV$", PATH_STAGE, ignore.case=TRUE) ~ "IV",
      TRUE ~ "Unknown"
    ),
    STAGE = factor(STAGE, levels = c("I", "II", "III", "IV", "Unknown")),
    T_Stage = case_when(
      grepl("^T1|^T1a|^T1b", PATH_T_STAGE, ignore.case=TRUE) ~ "T1",
      grepl("^T2|^T2a|^T2b", PATH_T_STAGE, ignore.case=TRUE) ~ "T2",
      grepl("^T3", PATH_T_STAGE, ignore.case=TRUE) ~ "T3",
      grepl("^T4", PATH_T_STAGE, ignore.case=TRUE) ~ "T4",
      TRUE ~ "TX/Unknown"
    ),
    T_Stage = factor(T_Stage, levels = c("T1", "T2", "T3", "T4", "TX/Unknown")),
    N_Stage = case_when(
      grepl("^N0", PATH_N_STAGE, ignore.case=TRUE) ~ "N0",
      grepl("^N1", PATH_N_STAGE, ignore.case=TRUE) ~ "N1",
      grepl("^N2", PATH_N_STAGE, ignore.case=TRUE) ~ "N2",
      grepl("^N3", PATH_N_STAGE, ignore.case=TRUE) ~ "N3",
      TRUE ~ "NX/Unknown"
    ),
    N_Stage = factor(N_Stage, levels = c("N0", "N1", "N2", "N3", "NX/Unknown")),
    M_Stage = case_when(
      grepl("^M0", PATH_M_STAGE, ignore.case=TRUE) ~ "M0",
      grepl("^M1|^M1a|^M1b", PATH_M_STAGE, ignore.case=TRUE) ~ "M1",
      TRUE ~ "MX/Unknown"
    ),
    M_Stage = factor(M_Stage, levels = c("M0", "M1", "MX/Unknown")),
    SEX = factor(SEX, levels = c("Male", "Female")),
    Subtype = case_when(
      grepl("Acinar", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Acinar",
      grepl("Papillary", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Papillary",
      grepl("Solid", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Solid",
      grepl("mucinous|colloid", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Mucinous",
      grepl("Bronchioloalveolar Carcinoma Nonmucinous", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Lepidic / BAC",
      grepl("Mixed", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "Mixed",
      grepl("NOS", PRIMARY_DIAGNOSIS, ignore.case = TRUE) ~ "NOS",
      is.na(PRIMARY_DIAGNOSIS) | PRIMARY_DIAGNOSIS == "" | PRIMARY_DIAGNOSIS == "PRIMARY_DIAGNOSIS" ~ "Unknown",
      TRUE ~ "Other"
    ),
    PRIOR_TREATMENT = case_when(
      grepl("true", PRIOR_TREATMENT, ignore.case=TRUE) ~ "Yes",
      grepl("false", PRIOR_TREATMENT, ignore.case=TRUE) ~ "No",
      TRUE ~ "Unknown"
    ),
    PRIOR_TREATMENT = factor(PRIOR_TREATMENT, levels = c("No", "Yes", "Unknown")),
    PRIOR_MALIGNANCY = case_when(
      grepl("true", PRIOR_MALIGNANCY, ignore.case=TRUE) ~ "Yes",
      grepl("false", PRIOR_MALIGNANCY, ignore.case=TRUE) ~ "No",
      TRUE ~ "Unknown"
    ),
    PRIOR_MALIGNANCY = factor(PRIOR_MALIGNANCY, levels = c("No", "Yes", "Unknown")),
    VITAL_STATUS = case_when(
      VITAL_STATUS == "Alive" ~ "Alive",
      VITAL_STATUS == "Dead" ~ "Deceased",
      TRUE ~ "Unknown"
    ),
    VITAL_STATUS = factor(VITAL_STATUS, levels = c("Alive", "Deceased", "Unknown")),
    TMB_NONSYNONYMOUS = as.numeric(TMB_NONSYNONYMOUS),
    SMOKING_PACK_YEARS = as.numeric(SMOKING_PACK_YEARS)
  )

tcga_vars <- list(
  list(col = "AGE", type = "continuous", label = "Age (Years)"),
  list(col = "SEX", type = "categorical", label = "Gender", palette = c("Male" = "#2C7FB8", "Female" = "#D53F8C")),
  list(col = "STAGE", type = "categorical", label = "Pathological Stage", palette = c("I" = "#7FCDBB", "II" = "#41B6C4", "III" = "#1D91C0", "IV" = "#081D58", "Unknown" = "#E2E8F0")),
  list(col = "T_Stage", type = "categorical", label = "T Stage", palette = c("T1" = "#C7E9B4", "T2" = "#7FCDBB", "T3" = "#41B6C4", "T4" = "#1D91C0", "TX/Unknown" = "#E2E8F0")),
  list(col = "N_Stage", type = "categorical", label = "N Stage", palette = c("N0" = "#C7E9B4", "N1" = "#7FCDBB", "N2" = "#41B6C4", "N3" = "#1D91C0", "NX/Unknown" = "#E2E8F0")),
  list(col = "M_Stage", type = "categorical", label = "M Stage", palette = c("M0" = "#C7E9B4", "M1" = "#1D91C0", "MX/Unknown" = "#E2E8F0")),
  list(col = "Subtype", type = "categorical", label = "Adenocarcinoma Subtype"),
  list(col = "PRIOR_TREATMENT", type = "categorical", label = "Prior Treatment History", palette = c("No" = "#A0AEC0", "Yes" = "#319795", "Unknown" = "#E2E8F0")),
  list(col = "PRIOR_MALIGNANCY", type = "categorical", label = "Prior Malignancy History", palette = c("No" = "#A0AEC0", "Yes" = "#D53F8C", "Unknown" = "#E2E8F0")),
  list(col = "VITAL_STATUS", type = "categorical", label = "Vital Status", palette = c("Alive" = "#38A169", "Deceased" = "#E53E3E", "Unknown" = "#E2E8F0")),
  list(col = "TMB_NONSYNONYMOUS", type = "continuous", label = "TMB (Nonsynonymous)", fill = "#D95F02"),
  list(col = "SMOKING_PACK_YEARS", type = "continuous", label = "Smoking Pack Years", fill = "#319795")
)

create_clinical_dashboard(tcga_clean, "TCGA", tcga_vars, "Plots/Clinical/TCGA_clinical_dashboard.png")

# 3. CHINA CLINICAL & PATHOLOGICAL DASHBOARD
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
 
china_vars <- list(
  list(col = "AGE", type = "continuous", label = "Age (Years)"),
  list(col = "SEX", type = "categorical", label = "Gender", palette = c("Male" = "#2C7FB8", "Female" = "#D53F8C")),
  list(col = "STAGE", type = "categorical", label = "Pathological Stage", palette = c("0" = "#D95F02", "I" = "#7FCDBB", "I-II" = "#41B6C4", "II" = "#1D91C0", "III" = "#225EA8", "III-IV" = "#253494", "IV" = "#081D58", "Unknown" = "#E2E8F0")),
  list(col = "SMOKE_STATUS", type = "categorical", label = "Smoking Status", palette = c("Nonsmoker" = "#319795", "Smoker" = "#ED8936", "Unknown" = "#E2E8F0")),
  list(col = "Treatment_Cleaned", type = "categorical", label = "Treatment Profile", palette = c("Treatment Naive" = "#319795", "Chemotherapy Only" = "#2C7FB8", "Targeted Therapy Only" = "#805AD5", "Multi-modal / Other" = "#ED8936", "Unknown" = "#E2E8F0")),
  list(col = "SAMPLE_TYPE", type = "categorical", label = "Sample Type", palette = c("Primary" = "#319795", "Metastasis" = "#ED8936", "Recurrent" = "#6B46C1")),
  list(col = "SPECIMEN_TYPE", type = "categorical", label = "Specimen Type", palette = c("Surgery" = "#319795", "Biopsy / Paracentesis" = "#ED8936")),
  list(col = "TUMOR_PURTITY", type = "continuous", label = "Tumor Purity", fill = "#319795"),
  list(col = "TMB_NONSYNONYMOUS", type = "continuous", label = "TMB (Nonsynonymous)", fill = "#D95F02")
)

create_clinical_dashboard(china_clean, "China", china_vars, "Plots/Clinical/China_clinical_dashboard.png")

# 4. MSK CLINICAL & PATHOLOGICAL DASHBOARD
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
    TUMOR_PURITY = ifelse(TUMOR_PURITY > 100, NA_real_, TUMOR_PURITY),
    FACETS_PURITY = as.numeric(FACETS_PURITY),
    FACETS_PURITY = ifelse(FACETS_PURITY > 1, NA_real_, FACETS_PURITY),
    FACETS_PLOIDY = as.numeric(FACETS_PLOIDY),
    FACETS_WGD = factor(FACETS_WGD, levels = c("FALSE", "TRUE"))
  )

msk_vars <- list(
  list(col = "AGE", type = "continuous", label = "Age (Years)"),
  list(col = "SEX", type = "categorical", label = "Gender", palette = c("Male" = "#2C7FB8", "Female" = "#D53F8C")),
  list(col = "Sample_Type", type = "categorical", label = "Sample Type", palette = c("Primary" = "#319795", "Metastasis" = "#ED8936", "Local Recurrence" = "#6B46C1")),
  list(col = "Met_Site", type = "categorical", label = "Metastatic Site Profile", palette = c(
    "Primary Tumor (N/A)" = "#A0AEC0", "Lymph Node" = "#2C7FB8", "Brain" = "#E53E3E",
    "Bone" = "#ED8936", "Liver" = "#319795", "Pleura" = "#805AD5",
    "Adrenal" = "#D53F8C", "Soft Tissue" = "#38A169", "Other Metastasis" = "#CBD5E0",
    "Unknown/Other" = "#E2E8F0"
  )),
  list(col = "GENE_PANEL", type = "categorical", label = "MSK Gene Panel"),
  list(col = "MSI_SCORE", type = "continuous", label = "MSI Score", fill = "#319795"),
  list(col = "MSI_TYPE", type = "categorical", label = "MSI Status", palette = c("Stable" = "#38A169", "Indeterminate" = "#ED8936", "Instable" = "#E53E3E", "Do not report" = "#A0AEC0")),
  list(col = "TMB_SCORE", type = "continuous", label = "TMB Score", fill = "#D95F02"),
  list(col = "TUMOR_PURITY", type = "continuous", label = "Tumor Purity", fill = "#319795"),
  list(col = "FACETS_PURITY", type = "continuous", label = "FACETS Purity", fill = "#4299E1"),
  list(col = "FACETS_PLOIDY", type = "continuous", label = "FACETS Ploidy", fill = "#4FD1C5")
)

create_clinical_dashboard(msk_clean, "MSK", msk_vars, "Plots/Clinical/MSK_clinical_dashboard.png")


##==============================================
#             Mutational Analysis
#===============================================

## Retrieve mutation data.

# Read MAF (Mutation Annotation Format) data for China, Singapore, and MSK cohorts
# and associate them with their respective clinical data
china_mut <- read.maf("public_data/china_pancan_2020/data_mutations.txt", clinicalData = china_mut_clin_merge)

sg_mut <- read.maf("public_data/singapore_luad_2020/data_mutations.txt", clinicalData = sg_clin_merge)

msk_mut <- read.maf("public_data/msk_impact_50k_2026/data_mutations.txt", clinicalData = msk_clin_merge)

tcga_mut <- read.maf("public_data/TCGA/luad_tcga_gdc/data_mutations.txt", clinicalData = tcga_clin_merge)

## Subset Lung cancer only

# Filter China mutation data to include only LUAD
lung_china_tsb <- china_mut_clin_merge %>%
  filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma") %>% 
  pull(Tumor_Sample_Barcode)
china_mut <- subsetMaf(maf = china_mut, tsb = lung_china_tsb)

# Filter MSK mutation data to include only LUAD
lung_msk_tsb <- msk_clin_merge %>%
  filter(CANCER_TYPE_DETAILED == "Lung Adenocarcinoma") %>% 
  pull(Tumor_Sample_Barcode)
msk_mut <- subsetMaf(maf = msk_mut, tsb = lung_msk_tsb)

# KNC percentage analysis
# Create a list of mutation data for all cohorts (KEAP1, NFE2L2, CUL3 analysis)
cohorts <- list(
  MSK_LUAD = msk_mut@data,
  China_LUAD = china_mut@data,
  TCGA_LUAD = tcga_mut@data,
  Singapore_LUAD = sg_mut@data
) 

# Define a function to calculate the number and percentage of patients with KNC mutations
# KNC mutations are defined as mutations in KEAP1, NFE2L2, or CUL3 genes
knc_percent <- function(file) {
  
  # Calculate total number of unique patients (samples) in the cohort
  total_mutations <- n_distinct(
    file$Tumor_Sample_Barcode
  )
  
  # Calculate number of unique patients with mutations in KNC genes
  knc_mutations <- n_distinct(
    file$Tumor_Sample_Barcode[
      file$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")
    ]
  )
  
  # Return results as a data frame
  data.frame(
    total_mutations = total_mutations,
    knc_mutations = knc_mutations,
    knc_percent = 100 * knc_mutations / total_mutations
  )
}

# Apply the knc_percent function to each cohort and combine results into a single data frame
knc_dist <- map_dfr(cohorts, 
                    knc_percent,
                    .id = "Cohorts")

  
# ==============================================
#            KNC Stacked Bar Plot
# ==============================================

# 1. Reshape the data for a 100% percentage stack visualization
plot_data <- knc_dist %>%
  mutate(
    Cohorts = factor(
      Cohorts,
      levels = c(
        "China_LUAD",
        "Singapore_LUAD",
        "TCGA_LUAD",
        "MSK_LUAD"
      )
    ),
    Non_KNC = 100 - knc_percent,
    KNC = knc_percent
  ) %>%
  pivot_longer(
    cols = c(Non_KNC, KNC),
    names_to = "Patient_Group",
    values_to = "Percentage"
  ) %>%
  mutate(
    Patient_Group = factor(
      Patient_Group,
      levels = c("KNC", "Non_KNC")
    )
  )

# 2. Generate and save the stacked bar plot
png(filename = "Plots/Mutational/KNC_percentage.png", width = 12, height = 8, units = "in", res = 600)

ggplot(plot_data, aes(x = Cohorts, y = Percentage, fill = Patient_Group)) +
  # Create the stacked bars
  geom_col(width = 0.4, alpha = 0.95, color = "white", linewidth = 0.4) +
  
  # Add text labels for KNC percentage and raw counts (e.g., "15.2% (100/658)")
  geom_text(
    aes(label = ifelse(Patient_Group == "KNC", sprintf("%.1f%%\n(%d/%d)",
                                                       Percentage,
                                                       knc_mutations,
                                                       total_mutations), "")),
    position = position_stack(vjust = 1.0), 
    vjust = -0.5,                          
    color = "black",                       
    fontface = "bold",
    size = 5
  ) +
  
  # Manually set colors for KNC and Non-KNC groups
  scale_fill_manual(
    values = c("Non_KNC" = "#E2E8F0", "KNC" = "#2C7FB8"),
    labels = c("Non_KNC" = "Other Patients", "KNC" = "KNC Patients")
  ) +
  
  # Format the Y-axis (0-100%) with extra space at the top for labels
  scale_y_continuous(
    limits = c(0, 115), # Expanded slightly past 100 so the top labels don't get clipped
    expand = c(0, 0),
    breaks = seq(0, 115, by = 10),
    labels = function(x) paste0(x, "%")
  ) +
  
  # Set plot titles and axis labels
  labs(
    title = "KNC Relative Prevalence Across Cohorts",
    x = NULL,
    y = "Percentage of Patients",
    fill = NULL
  ) +
  
  # Legend formatting
  guides(fill = guide_legend(reverse = FALSE)) +
  
  # Apply a clean publication-ready theme
  theme_pubr(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 15), hjust = 0.5), # Center-aligned title
    axis.title.y = element_text(face = "bold", margin = margin(r = 15), color = "#444444"),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.6),
    axis.line.x = element_line(color = "black")
  )

dev.off()

# ==============================================
#            KNC survival
# ==============================================

knc_survival_plot <- function(
    maf,
    time_col,
    status_col,
    cohort_name,
    outfile = NULL
){
  
  clin <- maf@clinical.data
  
  # Create KNC status
  knc_pax <- unique(
    maf@data$Tumor_Sample_Barcode[
      maf@data$Hugo_Symbol %in%
        c("KEAP1", "NFE2L2", "CUL3")
    ]
  )
  
  clin$KNC_status <- ifelse(
    clin$Tumor_Sample_Barcode %in% knc_pax,
    "KNC Mutated",
    "KNC Wild-Type"
  )
  
  # Convert survival status
  if (is.character(clin[[status_col]])) {
    
    clin[[status_col]] <- ifelse(
      clin[[status_col]] %in%
        c("DECEASED", "Dead", "DEAD", "1", "1:DECEASED"),
      1,
      0
    )
  }
  
  # Clean survival data
  clin <- clin %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[status_col]]),
      .data[[status_col]] != ""
    )
  
  clin[[time_col]] <- as.numeric(clin[[time_col]])
  clin[[status_col]] <- as.numeric(clin[[status_col]])
  
  # Survival fit
  # 1. Construct the formula string and convert to formula
  form_str <- paste0("Surv(", time_col, ", ", status_col, ") ~ KNC_status")
  surv_form <- as.formula(form_str)
  
  fit <- survfit(surv_form, data = clin)
  
  fit$call$formula <- surv_form
  
  # Log-rank test
  pval <- survdiff(surv_form, data = clin)
  
  # Plot
  surv_plot <- ggsurvplot(
    fit,
    data = clin,
    
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    
    title = paste(
      "Overall Survival by KNC Mutation Status in",
      cohort_name,
      "cohort"
    ),
    
    font.title = c(16, "bold"),
    
    xlab = "Overall Survival (Months)",
    ylab = "Survival Probability",
    
    legend.title = NULL,
    legend.labs = c(
      "KNC Mutated",
      "KNC Wild-Type"
    ),
    
    legend = "top",
    
    palette = c(
      "#2C7FB8",
      "#D95F02"
    ),
    
    surv.median.line = "hv",
    
    censor = TRUE,
    censor.shape = 124,
    censor.size = 3,
    
    ggtheme = theme_pubr(base_size = 14) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        )
      ),
    
    conf.int.alpha = 0.15
  )
  
  if (!is.null(outfile)) {
    
    png(
      outfile,
      width = 12,
      height = 9,
      units = "in",
      res = 600
    )
    
    print(surv_plot)
    dev.off()
  }
  
  return(
    list(
      fit = fit,
      survdiff = pval,
      plot = surv_plot,
      data = clin
    )
  )
}

msk_os <- knc_survival_plot(
  maf = msk_mut,
  time_col = "OS_MONTHS",
  status_col = "OS_STATUS",
  cohort_name = "MSK",
  outfile = "Plots/Survival/MSK_KNC_survival.png"
)

sg_os <- knc_survival_plot(
  maf = sg_mut,
  time_col = "OS_MONTHS",
  status_col = "OS_STATUS",
  cohort_name = "Singapore",
  outfile = "Plots/Survival/SG_KNC_survival.png"
)

tcga_os <- knc_survival_plot(
  maf = tcga_mut,
  time_col = "OS_MONTHS",
  status_col = "OS_STATUS",
  cohort_name = "TCGA",
  outfile = "Plots/Survival/TCGA_KNC_survival.png"
)
# ==============================================
#            KNC co-mutational analysis
# ==============================================

## Plot Co-occurence plot for all cohorts

png("Plots/Mutational/MSK_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
  msk_com <- somaticInteractions(msk_mut,
                                 topMar = 5.5,
                                 leftMar = 5.5)
dev.off()

write.csv(msk_com, "Tables/MSK_co-mutation.csv")

png("Plots/Mutational/TCGA_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
  tcga_com <- somaticInteractions(tcga_mut,
                                  topMar = 5.5,
                                  leftMar = 5.5)
dev.off()

write.csv(tcga_com, "Tables/TCGA_co-mutation.csv")

png("Plots/Mutational/China_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
  china_com <- somaticInteractions(china_mut,
                                 topMar = 5.5,
                                 leftMar = 5.5)
dev.off()

write.csv(china_com, "Tables/China_co-mutation.csv")

png("Plots/Mutational/Singapore_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
  sg_com <- somaticInteractions(sg_mut,
                                  topMar = 5.5,
                                  leftMar = 5.5)
dev.off()

write.csv(sg_com, "Tables/Singapore_co-mutation.csv")

## Calculate number of pax with combinations of knc gene mutations

mutationMatrix <- function(maf) {
  maf@data %>%
    distinct(
      Hugo_Symbol,
      Tumor_Sample_Barcode
    ) %>%
    mutate(value = 1) %>%
    pivot_wider(
      names_from = Tumor_Sample_Barcode,
      values_from = value,
      values_fill = 0
    ) %>% 
    t() %>% 
    row_to_names(row_number = 1)
}

msk_mm <- as.data.frame(mutationMatrix(maf = msk_mut))

keap_nrf2_msk_mm <- msk_mm %>% 
  dplyr::filter(KEAP1 == "1" & NFE2L2 == "1") %>% 
  dplyr::select(KEAP1, NFE2L2)

keap_cul3_msk_mm <- msk_mm %>% 
  dplyr::filter(KEAP1 == "1" & CUL3 == "1") %>% 
  dplyr::select(KEAP1, CUL3)

cul3_nrf2_msk_mm <- msk_mm %>% 
  dplyr::filter(CUL3 == "1" & NFE2L2 == "1") %>% 
  dplyr::select(CUL3, NFE2L2)

knc_msk_mm <- msk_mm %>% 
  dplyr::filter(KEAP1 == "1" & NFE2L2 == "1" & CUL3 == "1") %>% 
  dplyr::select(KEAP1, NFE2L2, CUL3)

## Lolipop plots for knc genes
plot_lollipop <- function(maf,
                          gene,
                          cohort_name,
                          outdir = "Plots",
                          top_n_labels = 5,
                          width = 15,
                          height = 5){
  
  # Find recurrent positions
  gene_pos <- maf@data %>%
    filter(Hugo_Symbol == gene) %>%
    mutate(
      pos = as.numeric(
        stringr::str_extract(HGVSp_Short, "\\d+")
      )
    ) %>%
    filter(!is.na(pos)) %>%
    count(pos, sort = TRUE)
  
  top_pos <- head(gene_pos$pos, top_n_labels)
  
  png(
    file.path(
      outdir,
      paste0(cohort_name, "_", gene, "_lollipop.png")
    ),
    width = width,
    height = height,
    units = "in",
    res = 600
  )
  
  par(oma = c(0, 0, 4, 0))
  
  lollipopPlot(
    maf = maf,
    gene = gene,
    showMutationRate = FALSE,
    titleSize = c(0.001, 0.001)
  )
  
  mtext(
    paste0(
      gene,
      " Mutational Landscape in ",
      cohort_name,
      " Cohort"
    ),
    side = 3,
    outer = TRUE,
    line = 1,
    cex = 1.5,
    font = 2
  )
  
  dev.off()
  
  invisible(top_pos)
}

lapply(
  c("KEAP1", "NFE2L2", "CUL3"),
  function(g)
    plot_lollipop(
      maf = msk_mut,
      gene = g,
      cohort_name = "MSK",
      outdir = "Plots/Mutational"
    )
)

lapply(
  c("KEAP1", "NFE2L2", "CUL3"),
  function(g)
    plot_lollipop(
      maf = tcga_mut,
      gene = g,
      cohort_name = "TCGA",
      outdir = "Plots/Mutational"
    )
)

lapply(
  c("KEAP1", "NFE2L2", "CUL3"),
  function(g)
    plot_lollipop(
      maf = sg_mut,
      gene = g,
      cohort_name = "Singapore",
      outdir = "Plots/Mutational"
    )
)

lapply(
  c("KEAP1", "NFE2L2", "CUL3"),
  function(g)
    plot_lollipop(
      maf = china_mut,
      gene = g,
      cohort_name = "China",
      outdir = "Plots/Mutational"
    )
)

# ==============================================
#            KNC variant analysis
# ==============================================

## Subset the maf files to keep only knc variants

knc_variant_survival <- function(
    maf,
    genes,
    cohort_name,
    time_col = "OS_MONTHS",
    status_col = "OS_STATUS",
    outdir = "Plots",
    top_n = 10
){
  
  # ============================
  # Subset MAF
  # ============================
  
  maf_sub <- subsetMaf(
    maf = maf,
    genes = genes
  )
  
  # ============================
  # Variant frequencies
  # ============================
  
  variant_freq <- maf_sub@data %>%
    filter(!is.na(HGVSp_Short)) %>%
    count(
      Hugo_Symbol,
      HGVSp_Short,
      name = "Frequency",
      sort = TRUE
    )
  
  # ============================
  # Top recurrent variants
  # Must be recurrent (>1)
  # Then keep top N
  # ============================
  
  top_variants <- variant_freq %>%
    filter(Frequency > 1) %>%
    slice_max(
      Frequency,
      n = top_n,
      with_ties = TRUE
    )
  
  if(nrow(top_variants) == 0){
    
    message(
      cohort_name,
      ": No recurrent variants found."
    )
    
    return(NULL)
    
  }
  
  # ============================
  # Patients carrying top variants
  # ============================
  
  top_patients <- unique(
    maf_sub@data$Tumor_Sample_Barcode[
      maf_sub@data$HGVSp_Short %in%
        top_variants$HGVSp_Short
    ]
  )
  
  # ============================
  # Clinical data
  # ============================
  
  surv_df <- maf_sub@clinical.data
  
  surv_df$var_prevalence <- ifelse(
    surv_df$Tumor_Sample_Barcode %in% top_patients,
    "Top recurrent",
    "Other variants"
  )
  
  cat("\n")
  cat("Cohort:", cohort_name, "\n")
  
  cat("Before filtering:", nrow(surv_df), "\n")
  
  print(
    summary(surv_df[[time_col]])
  )
  
  print(
    table(
      surv_df[[status_col]],
      useNA = "always"
    )
  )
  
  surv_df <- surv_df[
    !is.na(surv_df[[time_col]]) &
      !is.na(surv_df[[status_col]]) &
      surv_df[[status_col]] != "",
  ]
  
  cat("After filtering:", nrow(surv_df), "\n")
  
  surv_df$TIME <- as.numeric(
    surv_df[[time_col]]
  )
  
  surv_df$STATUS <- ifelse(
    surv_df[[status_col]] %in%
      c(
        "DECEASED",
        "Dead",
        "dead",
        1
      ),
    1,
    0
  )
  
  surv_df$var_prevalence <- factor(
    surv_df$var_prevalence,
    levels = c(
      "Other variants",
      "Top recurrent"
    )
  )
  
  cat("\n")
  cat(cohort_name, "\n")
  print(table(surv_df$var_prevalence))
  cat("\n")
  
  if(length(unique(surv_df$var_prevalence)) < 2){
    
    message(
      cohort_name,
      ": only one group present."
    )
    
    return(NULL)
    
  }
  
  # ============================
  # KM fit
  # ============================
  
  fit <- survfit(
    Surv(TIME, STATUS) ~ var_prevalence,
    data = surv_df
  )
  
  # ============================
  # Plot
  # ============================
  
  p <- ggsurvplot(
    fit,
    data = surv_df,
    
    pval = TRUE,
    pval.method = TRUE,
    
    conf.int = TRUE,
    
    risk.table = TRUE,
    risk.table.height = 0.25,
    
    title = paste0(
      cohort_name,
      ": Survival by KNC Variant Recurrence"
    ),
    
    legend.title = NULL,
    
    legend.labs = c(
      "Other variants",
      "Top recurrent variants"
    ),
    
    ggtheme = theme_pubr(base_size = 14)
  )
  
  png(
    file.path(
      outdir,
      paste0(
        cohort_name,
        "_variant_survival.png"
      )
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 600
  )
  
  print(p)
  
  dev.off()
  
  return(
    list(
      maf_subset = maf_sub,
      variant_frequency = variant_freq,
      top_variants = top_variants,
      top_patients = top_patients,
      survival_data = surv_df,
      fit = fit,
      plot = p
    )
  )
  
}

cohorts <- list(
  MSK = msk_mut,
  China = china_mut,
  Singapore = sg_mut,
  TCGA = tcga_mut
)

results <- lapply(
  names(cohorts),
  function(x)
    knc_variant_survival(
      maf = cohorts[[x]],
      genes = c("KEAP1","NFE2L2","CUL3"),
      cohort_name = x,
      outdir = "Plots/Survival"
    )
)

names(results) <- names(cohorts)