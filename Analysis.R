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
})

# Set working directory to the project root
setwd("/mnt/Linux_storage/KNC")


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
  TCGA_LUAD = TCGA_mut,
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
png(filename = "Plots/KNC_percentage.png", width = 12, height = 8, units = "in", res = 600)

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

msk_mut@clinical.data$OS_STATUS <- ifelse(msk_mut@clinical.data$OS_STATUS == "DECEASED", 1, 0)
knc_pax <- unique(
  msk_mut@data$Tumor_Sample_Barcode[
    msk_mut@data$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")
  ]
)
msk_surv_df <- msk_mut@clinical.data
msk_surv_df$KNC_status <- ifelse(
  msk_surv_df$Tumor_Sample_Barcode %in% knc_pax,
  "KNC Mutated",
  "No KNC mutation"
)
msk_surv_df <- msk_surv_df %>%
  filter(
    !is.na(OS_MONTHS),
    !is.na(OS_STATUS),
    OS_STATUS != ""
  )
msk_surv_df$OS_MONTHS <- as.numeric(msk_surv_df$OS_MONTHS)
msk_fit <- survfit(
  Surv(event = OS_STATUS, time = OS_MONTHS) ~ KNC_status,
  data = msk_surv_df
)
msk_p_val <- survdiff(
  Surv(OS_MONTHS, OS_STATUS) ~ KNC_status,
  data = msk_surv_df
)
png("Plots/MSK_KNC_mutVSwt_survival.png", width = 12, height = 9, units = "in",
    res = 600)
ggsurvplot(
  msk_fit,
  data = msk_surv_df,
  # Statistics
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  # Risk table
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  # Titles
  title = "Overall Survival by KNC Mutation Status in MSK cohort",
  font.title = c(16, "bold"),
  xlab = "Overall Survival (Months)",
  ylab = "Survival Probability",
  # Legend
  legend.title = NULL,
  legend.labs = c("KNC Mutated", "KNC Wild-Type"),
  legend = "top",
  # Colors
  palette = c("#2C7FB8", "#D95F02"),
  # Median survival lines
  surv.median.line = "hv",
  # Censor marks
  censor = TRUE,
  censor.shape = 124,
  censor.size = 3,
  # Theme
  ggtheme = theme_pubr(base_size = 14) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    ),
  # Confidence interval transparency
  conf.int.alpha = 0.15
)
dev.off()

# ==============================================
#            KNC co-mutational analysis
# ==============================================

com_msk <- somaticInteractions(msk_mut)

keap1_com_msk <- com_msk %>% 
  filter(gene1 %in% c("KEAP1", "NFE2L2", "CUL3") | gene2 %in% c("KEAP1", "NFE2L2", "CUL3"))
 