#!/usr/bin/Rscript

# ==============================================================================
# PIPELINE: pNRF2 vs. non_pNRF2 Differential Expression & Survival Analysis
#
# Description:
#   This script identifies Differentially Expressed Genes (DEGs) and proteins
#   associated with phospho-NRF2 abundance using CPTAC LUAD
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
  library(org.Hs.eg.db)
})

# Set working directory to project root
setwd("/media/nannu1375/Backpack/Shankara/KNC")

# ==============================================================================
# 1. IDENTIFY pNRF2-POSITIVE PATIENTS IN CPTAC LUAD
# ==============================================================================

# Load CPTAC phosphoprotein dataset
cptac_phospho <- readRDS("Tables/cptac_phospho.rds")

# Extract the NFE2L2 (NRF2) phosphosites
nrf2_phospho <- cptac_phospho %>%
  filter(NAME %in% c("NFE2L2_S215s", "NFE2L2_S433s")) %>%
  as.data.frame()

# Extract patient columns (excluding metadata columns 1 to 4)
patient_cols <- colnames(nrf2_phospho)[-c(1,2,3,4)]

# Determine if each patient has any detected pNRF2 (non-NA in S215s or S433s)
val_s215 <- as.numeric(nrf2_phospho[nrf2_phospho$NAME == "NFE2L2_S215s", patient_cols])
val_s433 <- as.numeric(nrf2_phospho[nrf2_phospho$NAME == "NFE2L2_S433s", patient_cols])

# Identify patient IDs that are pNRF2 positive
pnrf2_pax <- patient_cols[!is.na(val_s215) | !is.na(val_s433)]

# ==============================================================================
# 2. LOAD & PREPROCESS CPTAC RNA-SEQ COUNT DATA
# ==============================================================================

## Read RNA-seq raw counts from GDC dataset
# Read raw read counts using fast data.table reader
raw_counts <- data.table::fread(
  "public_data/luad_cptac_gdc/data_mrna_seq_read_counts.txt",
  nThread = 8
)

## Convert Entrez IDs to numeric
raw_counts$Entrez_Gene_Id <- as.numeric(raw_counts$Entrez_Gene_Id)

## Obtain official HGNC symbols from org.Hs.eg.db
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(unique(raw_counts$Entrez_Gene_Id)),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
) %>%
  dplyr::rename(
    Entrez_Gene_Id = ENTREZID,
    Gene = SYMBOL
  ) %>%
  dplyr::mutate(Entrez_Gene_Id = as.numeric(Entrez_Gene_Id)) %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::distinct(Entrez_Gene_Id, .keep_all = TRUE)

## Add gene symbols
raw_counts <- raw_counts %>%
  dplyr::left_join(gene_map, by = "Entrez_Gene_Id") %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::select(Gene, everything())

# Standardize GDC column names to Patient IDs (first 9 characters)
colnames(raw_counts) <- sapply(colnames(raw_counts), function(col) {
  if (col == "Gene") return(col)
  return(substr(col, 1, 9))
})

# Filter columns to keep only valid patient IDs from the phosphoprotein cohort
valid_cols <- colnames(raw_counts)[-1] %in% patient_cols
raw_counts <- raw_counts[, c(TRUE, valid_cols), with = FALSE]

## Remove duplicate patient columns (keep first aliquot for each patient)
duplicate_cols <- duplicated(colnames(raw_counts))
raw_counts <- raw_counts[, !duplicate_cols, with = FALSE]

# Group by Gene symbol and sum counts for duplicate symbol mappings
raw_counts <- raw_counts %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    across(
      everything(),
      sum
    ),
    .groups = "drop"
  )

## Move genes to rownames
raw_counts <- as.data.frame(raw_counts)
rownames(raw_counts) <- raw_counts$Gene
raw_counts$Gene <- NULL

# Set NA values in counts to 0
raw_counts[is.na(raw_counts)] <- 0

# ==============================================================================
# 3. DIFFERENTIAL GENE EXPRESSION (edgeR)
# ==============================================================================

## Create sample information
# Initialize metadata sheet mapping samples to groups (pNRF2 vs. non_pNRF2)
sample_info <- data.frame(
  Sample = colnames(raw_counts)
)

sample_info$Group <- ifelse(
  sample_info$Sample %in% pnrf2_pax,
  "pNRF2",
  "non_pNRF2"
)

# Set "non_pNRF2" as the baseline/reference group for statistical modeling
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("non_pNRF2", "pNRF2")
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
  levels = c("non_pNRF2", "pNRF2")
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

# Differential expression: Conduct QL F-test to compare pNRF2 group against the non_pNRF2 baseline (coefficient 2)
res <- glmQLFTest(
  fit,
  coef = 2
)

# Extract full edgeR statistics table
deg <- topTags(
  res,
  n = Inf
)$table

# Filter for significantly differentially expressed genes (FDR < 0.05 and |log2FC| >= 0.5)
sig_deg <- deg %>%
  filter(
    FDR < 0.05 &
      abs(logFC) >= 0.5
  )

# ==============================================================================
# 4. PLOTTING & SAVE DEG OUTPUTS
# ==============================================================================

# Save high-resolution Volcano Plot representing significant DEGs
png("Plots/CPTAC_pNRF2_volcano.png",
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
# 5. NORMALIZATION (DESeq2 VST)
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

vst_sig <- vst_mat[rownames(vst_mat) %in% rownames(sig_deg),
                   ]

expr_df <- as.data.frame(t(vst_sig))

# ==============================================================================
# 6. CLINICAL OUTCOME & SURVIVAL (OS) MERGING
# ==============================================================================

## Combine OS data to this matrix
# Load CPTAC patient survival metadata, skip descriptive headers, clean column names
os_df <- read.delim("public_data/luad_cptac_gdc/data_clinical_patient.txt") %>% 
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

# Write univariate cox results to Table
write.csv(cox_results, file = "Tables/cptac_univariate_cox_results.csv", row.names = FALSE)

## Get final gene signature
# Extract the names of the final prognostic genes
final_gene_sig <- cox_results$Gene

## Convert univariate cox results to z score
# Subset patient clinical outcomes and the prognostic genes
cox_df_z <- cox_df[, c("PATIENT_ID", final_gene_sig, "OS_STATUS", "OS_MONTHS")]

# Perform Z-score scaling on prognostic genes to standardize expression variance (improves multi-cox parameter comparisons)
cox_df_z[, final_gene_sig] <- scale(cox_df_z[, final_gene_sig])

# ==============================================================================
# 8. MULTIVARIATE COX REGRESSION
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

# Write multivariate cox results to Table
write.csv(multi_results, file = "Tables/cptac_multivariate_cox_results.csv", row.names = FALSE)

# Note: The risk score prediction can be extracted for clinical stratification (low/high-risk patients)
# using the linear predictor of the fitted model:
# Code for prediction for each patient
# cox_df$risk_score <- predict(
#   fit_multi,
#   type = "lp"
# )