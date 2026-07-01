#!/usr/bin/Rscript

# ==============================================================================
# PIPELINE: KNC Mutant vs. Non-KNC Differential Expression & Survival Analysis
#
# Description:
#   This script identifies Differentially Expressed Genes (DEGs) and proteins
#   associated with KNC pathway mutations (KEAP1, NFE2L2, CUL3) using TCGA LUAD
#   (RNA-seq) and CPTAC LUAD (Proteomics) datasets. It then fits univariate and
#   multivariate Cox proportional hazard models to establish a prognostic gene
#   signature.
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)            # Differential expression analysis of RNA-seq count data
  library(dplyr)            # Data manipulation grammar
  library(tidyverse)        # Collection of packages for data science
  library(AnnotationDbi)    # Interface to query annotation data
  library(EnhancedVolcano)  # Highly-customizable volcano plot generation
  library(biomaRt)          # Query BioMart databases (Ensembl)
  library(maftools)         # Analyze and visualize Mutation Annotation Format (MAF) data
  library(janitor)          # Data cleaning and table formatting
  library(limma)            # Linear models for microarray and proteomic analysis
  library(DESeq2)           # Differential gene expression & variance stabilization
  library(survival)         # Survival analysis (Cox regression, Kaplan-Meier)
})

# Set working directory to project root
setwd("/media/nannu1375/Backpack/Shankara/KNC")

# ==============================================================================
# 1. IDENTIFY KNC-MUTATED PATIENTS IN TCGA LUAD
# ==============================================================================

## Read the TCGA rds file and take KNC patients
# Load TCGA mutation data (MAF) for lung adenocarcinoma (LUAD)
tcgaRDS <- readRDS("Tables/tcga_maf_luad.rds")

# Extract Tumor Sample Barcodes for patients carrying mutations in key KNC pathway genes
knc_pax <- tcgaRDS@data$Tumor_Sample_Barcode[
  tcgaRDS@data$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")
]

# Convert sample barcodes (e.g., "TCGA-XX-XXXX-01A-...") to 12-character Patient IDs (e.g., "TCGA-XX-XXXX")
knc_pax <- substr(knc_pax, 1, 12)
knc_pax <- unique(knc_pax)

# ==============================================================================
# 2. LOAD & PREPROCESS TCGA RNA-SEQ COUNT DATA
# ==============================================================================

## Read RNA-seq counts
# Read raw RNA-seq counts using fast data.table reader, filtering for protein-coding genes
raw_counts <- data.table::fread(
  "public_data/TCGA/luad_tcga_pan_can_atlas_2018/raw_counts_LUAD.csv",
  sep = ",",
  nThread = 8
) %>%
  filter(gene_type == "protein_coding") %>% 
  dplyr::select(-gene_type)

## Keep only primary tumour samples (sample type = 01)
# Extract column names (representing TCGA barcodes) excluding the first column (Gene)
sample_cols <- colnames(raw_counts)[-1]

# Filter for primary solid tumor samples (indicated by "01" as characters 14-15 in barcode)
tumor_cols <- sample_cols[
  substr(sample_cols, 14, 15) == "01"
]

# Subset counts to keep only the primary tumor samples
raw_counts <- raw_counts %>%
  dplyr::select(
    Gene,
    all_of(tumor_cols)
  )

## Convert TCGA barcodes to patient IDs
# Truncate barcodes to 12-character Patient IDs to enable matching with clinical and mutation data
colnames(raw_counts)[-1] <- substr(
  colnames(raw_counts)[-1],
  1,
  12
)

## Remove duplicate tumour aliquots
# Remove duplicate tumor aliquots (keeping only the first column for each patient)
raw_counts <- raw_counts[
  ,
  !duplicated(names(raw_counts)),
  with = FALSE
]

## Collapse duplicate gene symbols
# Group by Gene symbol and sum counts to collapse duplicate gene features
raw_counts <- raw_counts %>%
  group_by(Gene) %>%
  summarise(
    across(
      everything(),
      sum
    ),
    .groups = "drop"
  )

## Move genes to rownames
# Convert to standard data.frame and move Gene names to row index (required by edgeR/DESeq2)
raw_counts <- as.data.frame(raw_counts)

rownames(raw_counts) <- raw_counts$Gene

raw_counts$Gene <- NULL

# ==============================================================================
# 3. DIFFERENTIAL GENE EXPRESSION (edgeR)
# ==============================================================================

## Create sample information
# Initialize metadata sheet mapping samples to groups (KNC vs. Non-KNC)
sample_info <- data.frame(
  Sample = colnames(raw_counts)
)

sample_info$Group <- ifelse(
  sample_info$Sample %in% knc_pax,
  "KNC",
  "Non-KNC"
)

# Set "Non-KNC" as the baseline/reference group for statistical modeling
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Non-KNC", "KNC")
)

sample_info$Group <- as.factor(sample_info$Group)

# Create DGE list for edgeR analysis
dge <- DGEList(
  counts = raw_counts,
  group = sample_info$Group
)

# Filter out low counts per million genes (CPM) to improve statistical power
keep <- filterByExpr(
  dge,
  sample_info$Group
)

# Keep only high CPM genes and recalculate library sizes
dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

# Normalise (TMM) to account for library composition differences
dge <- calcNormFactors(dge)

# Set the comparison groups and create a design matrix (~ Group)
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Non-KNC", "KNC")
)

design <- model.matrix(
  ~ Group,
  data = sample_info
)

# Estimate dispersion using Empirical Bayes
dge <- estimateDisp(
  dge,
  design
)

# Fit quasi-likelihood (QL) negative binomial generalized linear model
fit <- glmQLFit(
  dge,
  design
)

# Differential expression: Conduct QL F-test to compare KNC group against the Non-KNC baseline (coefficient 2)
res <- glmQLFTest(
  fit,
  coef = 2
)

# Extract full edgeR statistics table
deg <- topTags(
  res,
  n = Inf
)$table

# Filter for significantly differentially expressed genes (FDR < 0.05 and |log2FC| > 1.5)
sig_deg <- deg %>%
  filter(
    FDR < 0.05 &
      abs(logFC) > 1.5
  )

# ==============================================================================
# 4. PLOTTING & SAVE DEG OUTPUTS
# ==============================================================================

# Save high-resolution Volcano Plot representing significant DEGs
png("Plots/TCGA_KNC_volcano.png",
    width = 10,
    height = 10,
    units = "in",
    res = 600)

EnhancedVolcano(
  deg,
  lab = rownames(deg),
  x = "logFC",
  y = "FDR"
)

dev.off()

# Save significant DEGs metadata to CSV
write.csv(sig_deg, "Tables/Significant_genes.csv")

# ==============================================================================
# 5. CPTAC LUAD PROTEOMICS & limma ANALYSIS
# ==============================================================================

## Get protein data
# Load mass spectrometry protein quantification data, rename index column
prot_dat <- read.delim("public_data/luad_cptac_2020/data_protein_quantification.txt") %>%
  dplyr::rename(protein = prot) %>%
  column_to_rownames(var = "protein")

## Take mutation data and filter for KNC mutated vs non mutated
# Load CPTAC patient clinical characteristics, skipping descriptive header rows
cptac_patient <- read.delim("public_data/luad_cptac_2020/data_clinical_patient.txt") %>%
  dplyr::slice(-c(1,2,3)) %>%
  janitor::row_to_names(row_number = 1) %>%
  as.data.frame()

# Load CPTAC sample annotations
cptac_sample <- read.delim("public_data/luad_cptac_2020/data_clinical_sample.txt") %>%
  dplyr::slice(-c(1,2,3)) %>%
  janitor::row_to_names(row_number = 1) %>%
  as.data.frame()

# Merge clinical patient and sample files to link PATIENT_ID with SAMPLE_ID (Tumor_Sample_Barcode)
cptac_clin <- merge.data.frame(cptac_patient, cptac_sample, by = "PATIENT_ID") %>%
  dplyr::rename(Tumor_Sample_Barcode = SAMPLE_ID)

# Read mutations and attach clinical annotations using maftools
cptac_maf <- read.maf("public_data/luad_cptac_2020/data_mutations.txt",
                      clinicalData = cptac_clin)

# Extract CPTAC sample IDs for KNC-mutated patients
cptac_knc <- unique(
  cptac_maf@data$Tumor_Sample_Barcode[
    cptac_maf@data$Hugo_Symbol %in% c(
      "KEAP1",
      "NFE2L2",
      "CUL3"
    )
  ]
)

## Create a samplesheet of tsb and knc status for limma
# Build metadata mapping CPTAC sample columns to KNC status
prot_meta <- data.frame(
  Sample = colnames(prot_dat[, 1:ncol(prot_dat)])
)

# Standardize delimiter format: replace dots (.) in sample names with hyphens (-) to match MAF IDs
prot_meta$Sample <- gsub(
  pattern = "\\.",
  replacement = "-",
  x = prot_meta$Sample
)

# Assign groups and set reference factor
prot_meta$Group <- ifelse(
  prot_meta$Sample %in% cptac_knc,
  "KNC",
  "Non-KNC"
)

prot_meta$Group <- factor(
  prot_meta$Group,
  levels = c("Non-KNC", "KNC")
)

# Cleaning up data for limma
# Align protein data column names with the standardized metadata names
colnames(prot_dat) <- gsub(
  "\\.",
  "-",
  colnames(prot_dat)
)

prot_dat <- as.matrix(prot_dat)

# Create Design matrix for limma
design <- model.matrix(
  ~ Group,
  data = prot_meta
)

# Fit limma linear model for each protein feature
fit <- lmFit(
  prot_dat,
  design
)

# Moderate standard errors via Empirical Bayes method
fit <- eBayes(fit)

# Get protein info (adjusted p-value < 0.05)
prot_limma <- topTable(
  fit = fit,
  coef = "GroupKNC",
  number = Inf,
  p.value = 0.05
)

# Filter for differential proteins with a log2 fold-change threshold > 1.5
top_prot_limma <- prot_limma %>%
  filter(abs(logFC) > 1.5)

## Filter DEGs in top_prot_limma
# Intersect significant mRNA transcripts (TCGA DEGs) with significant proteins (CPTAC limma)
# to obtain a robust, cross-validated signature (concordant at both transcription and translation level)
final_deg <- sig_deg %>%
  filter(rownames(sig_deg) %in% rownames(top_prot_limma))

gene_sig <- rownames(final_deg)

# ==============================================================================
# 6. NORMALIZATION (DESeq2 VST)
# ==============================================================================

## Run VST normalisation on the raw counts
# Set up DESeq2 object from raw counts for normalization
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_info,
  design = ~ Group
)

# Perform Variance Stabilizing Transformation (vst) to homoscedasticize count data
vsd <- vst(dds, blind = FALSE)

vst_mat <- assay(vsd)

## Filter normalised data for the 14 genes
# Subset using the 14 protein-concordant genes signature
vst_14 <- vst_mat[
  rownames(vst_mat) %in% gene_sig,
]

expr_df <- as.data.frame(t(vst_14))

## Filter normalised data for the 14 genes
# NOTE: In the original pipeline design, this second filtering step overwrites vst_14 and expr_df 
# using all raw TCGA DEGs (rownames(sig_deg)) instead of the protein-cross-validated gene_sig.
vst_14 <- vst_mat[
  rownames(vst_mat) %in% rownames(sig_deg), 
  ]

expr_df <- as.data.frame(t(vst_14))

# ==============================================================================
# 7. CLINICAL OUTCOME & SURVIVAL (OS) MERGING
# ==============================================================================

## Combine OS data to this matrix
# Load TCGA patient survival metadata, skip descriptive headers, clean column names
os_df <- read.delim("public_data/TCGA/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt") %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame() %>% 
  dplyr::select(PATIENT_ID, OS_STATUS, OS_MONTHS) %>% 
  dplyr::filter(OS_MONTHS != "",
                !is.na(OS_MONTHS)) %>% 
  # Convert status string to binary event marker (1 = deceased, 0 = living)
  dplyr::mutate(OS_STATUS = ifelse(OS_STATUS == "1:DECEASED", 1, 0))

# Clean expression patient barcodes and merge with survival metrics
cox_df <- expr_df %>%
  rownames_to_column(var = "PATIENT_ID") %>%
  mutate(
    # Remove sample code suffix "-01A" to get patient ID matching clinical sheet
    PATIENT_ID = gsub("\\-01A$", "", PATIENT_ID)
  ) %>%
  merge(
    y = os_df,
    by = "PATIENT_ID"
  ) %>% 
  dplyr::mutate(
    OS_MONTHS = as.numeric(OS_MONTHS),
    OS_STATUS = as.numeric(OS_STATUS)
  )

# ==============================================================================
# 8. UNIVARIATE COX REGRESSION SURVIVAL ANALYSIS
# ==============================================================================

## Univariate COX regression
# Identify all gene columns in the merged dataset (excluding sample ID and clinical event outcomes)
genes <- colnames(cox_df)[
  !(colnames(cox_df) %in%
      c("PATIENT_ID", "OS_STATUS", "OS_MONTHS"))
]

# Run a univariate Cox Proportional Hazards model for each gene individually
cox_results <- lapply(
  genes,
  function(gene){
    
    fit <- coxph(
      as.formula(
        paste0(
          "Surv(OS_MONTHS, OS_STATUS) ~ `",
          gene,
          "`"
        )
      ),
      data = cox_df
    )
    s <- summary(fit)
    
    # Extract Hazard Ratio (HR), confidence limits, and Wald test p-value
    data.frame(
      Gene = gene,
      HR = s$coef[,"exp(coef)"],
      Lower95 = s$conf.int[,"lower .95"],
      Upper95 = s$conf.int[,"upper .95"],
      PValue = s$coef[,"Pr(>|z|)"]
    )
  }
)

# Combine list elements into a single data frame
cox_results <- do.call(rbind, cox_results)

# Control FDR using the Benjamini-Hochberg (BH) adjustment procedure
cox_results$FDR <- p.adjust(
  cox_results$PValue,
  method = "BH"
)

# Print sorted results by raw P-value
cox_results[order(cox_results$PValue), ]

# Filter to retain only prognostic genes (PValue < 0.05)
cox_results <- cox_results %>% 
  dplyr::filter(PValue < 0.05)

## Get final gene signature
# Extract the names of the final prognostic genes
final_gene_sig <- cox_results$Gene

## Convert univariate cox results to z score
# Subset patient clinical outcomes and the prognostic genes
cox_df_z <- cox_df[, c("PATIENT_ID", final_gene_sig, "OS_STATUS", "OS_MONTHS")]

# Perform Z-score scaling on prognostic genes to standardize expression variance (improves multi-cox parameter comparisons)
cox_df_z[, final_gene_sig] <- scale(cox_df_z[, final_gene_sig])

# ==============================================================================
# 9. MULTIVARIATE COX REGRESSION
# ==============================================================================

## Multivariate cox regression.
# Construct model formula containing all prognostic genes: Surv(time, event) ~ gene1 + gene2 + ...
cox_formula <- as.formula(
  paste(
    "Surv(OS_MONTHS, OS_STATUS) ~",
    paste(
      paste0("`", final_gene_sig, "`"),
      collapse = " + "
    )
  )
)

# Fit multivariate Cox model to assess independent prognostic significance
fit_multi <- coxph(
  formula = cox_formula,
  data = cox_df_z
)

# Structure and format coefficients, hazard ratios, and Wald stats for the report
multi_results <- data.frame(
  Gene = rownames(summary(fit_multi)$coefficients),
  Coefficient = summary(fit_multi)$coefficients[, "coef"],
  HR = summary(fit_multi)$coefficients[, "exp(coef)"],
  SE = summary(fit_multi)$coefficients[, "se(coef)"],
  Z = summary(fit_multi)$coefficients[, "z"],
  PValue = summary(fit_multi)$coefficients[, "Pr(>|z|)"],
  Lower95 = summary(fit_multi)$conf.int[, "lower .95"],
  Upper95 = summary(fit_multi)$conf.int[, "upper .95"]
)

# Identify genes that remain statistically significant under joint modeling (P < 0.05)
multi_results <- multi_results %>% 
  filter(PValue < 0.05)

# Note: The risk score prediction can be extracted for clinical stratification (low/high-risk patients)
# using the linear predictor of the fitted model:
# Code for prediction for each patient
# cox_df$risk_score <- predict(
#   fit_multi,
#   type = "lp"
# )