#!/usr/bin/Rscript

# ==============================================================================
# PIPELINE: KNC Mutant vs. Non-KNC Differential Expression & Survival Analysis (GSEA Branch)
#
# Description:
#   This script identifies Differentially Expressed Genes (DEGs) associated with
#   KNC pathway mutations (KEAP1, NFE2L2, CUL3) using TCGA LUAD RNA-seq. It runs
#   Fast Gene Set Enrichment Analysis (FGSEA) using Hallmark pathways, filters
#   for leading-edge genes, and fits univariate and multivariate Cox proportional
#   hazard models to establish a prognostic signature.
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
  library(fgsea)            # Fast Gene Set Enrichment Analysis
})

# Set working directory to project root
setwd("/mnt/Linux_storage/KNC")

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
# 5. GENE SET ENRICHMENT ANALYSIS (fgsea)
# ==============================================================================

## Run GSEA and take top pathway genes for filtering and further cox regression
# Load Hallmark pathway gene sets in GMT format
hallmark_pathway <- gmtPathways("public_data/h.all.v2026.1.Hs.symbols.gmt")

# Create genes ranked list
# Construct ranked gene vector using log2 fold-changes from edgeR DEG analysis
gene_ranks <- deg$logFC
names(gene_ranks) <- rownames(deg)

# Run GSEA
# Perform Fast Gene Set Enrichment Analysis (FGSEA) using Hallmark pathways
fgsea_res <- fgsea(
  pathways = hallmark_pathway,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500
)

# Filter for statistically significant enriched pathways (adjusted p-value < 0.05)
fgsea_res <- fgsea_res %>% 
  filter(padj < 0.05)

# Sort pathways by statistical significance
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

# Get list of genes in top pathways
# Extract the unique set of leading edge genes (core driving genes) from the top pathways
genes_list <- unique(unlist(fgsea_res$leadingEdge))

## Filter DEGs in top_prot_limma
# NOTE: The legacy comment below refers to top_prot_limma from the proteomics branch, but the
# functional code intersects the TCGA DEGs with the GSEA leading edge genes (genes_list) instead.
final_deg <- sig_deg %>%
  filter(rownames(sig_deg) %in% genes_list)

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
# Subset using the GSEA-cross-validated gene signature
vst_14 <- vst_mat[
  rownames(vst_mat) %in% gene_sig,
]

expr_df <- as.data.frame(t(vst_14))

## Filter normalised data for the 14 genes
# !!! WARNING: OVERWRITE BUG !!!
# The block below immediately overwrites vst_14/vst_37 and expr_df using all raw TCGA DEGs
# (rownames(sig_deg)) instead of keeping the GSEA-cross-validated gene signature.
vst_37 <- vst_mat[
  rownames(vst_mat) %in% rownames(sig_deg), 
]

expr_df <- as.data.frame(t(vst_37))

# ==============================================================================
# 7. CLINICAL OUTCOME & SURVIVAL (PFS/OS) MERGING
# ==============================================================================

## Combine PFS data to this matrix
# Load TCGA patient survival metadata, skip descriptive headers, clean column names
# NOTE: Even though this is named pfs_df and the comment refers to PFS (Progression-Free Survival),
# the code actually selects and processes Overall Survival metrics (OS_STATUS, OS_MONTHS).
pfs_df <- read.delim("public_data/TCGA/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt") %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame() %>% 
  dplyr::select(PATIENT_ID, PFS_STATUS, PFS_MONTHS) %>% 
  dplyr::filter(PFS_MONTHS != "",
                !is.na(PFS_MONTHS)) %>% 
  # Convert status string to binary event marker (1 = deceased, 0 = living)
  dplyr::mutate(PFS_STATUS = ifelse(PFS_STATUS == "1:DECEASED", 1, 0))

# Clean expression patient barcodes and merge with survival metrics
# !!! CRITICAL SCRIPT BUG !!!
# The merge step references the variable 'os_df' which is not defined in this script
# (the clinical data is instead assigned to 'pfs_df' above). Running this script in a fresh 
# R session will throw: Error in merge(...) : object 'os_df' not found.
cox_df <- expr_df %>%
  rownames_to_column(var = "PATIENT_ID") %>%
  mutate(
    # Remove sample code suffix "-01A" to get patient ID matching clinical sheet
    PATIENT_ID = gsub("\\-01A$", "", PATIENT_ID)
  ) %>%
  merge(
    y = pfs_df,
    by = "PATIENT_ID"
  ) %>% 
  dplyr::mutate(
    PFS_MONTHS = as.numeric(PFS_MONTHS),
    PFS_STATUS = as.numeric(PFS_STATUS)
  )

# ==============================================================================
# 8. UNIVARIATE COX REGRESSION SURVIVAL ANALYSIS
# ==============================================================================

## Univariate COX regression
# Identify all gene columns in the merged dataset (excluding sample ID and clinical event outcomes)
genes <- colnames(cox_df)[
  !(colnames(cox_df) %in%
      c("PATIENT_ID", "PFS_STATUS", "PFS_MONTHS"))
]

# Run a univariate Cox Proportional Hazards model for each gene individually
cox_results <- lapply(
  genes,
  function(gene){
    
    fit <- coxph(
      as.formula(
        paste0(
          "Surv(PFS_MONTHS, PFS_STATUS) ~ `",
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
cox_df_z <- cox_df[, c("PATIENT_ID", final_gene_sig, "PFS_STATUS", "PFS_MONTHS")]

# Perform Z-score scaling on prognostic genes to standardize expression variance (improves multi-cox parameter comparisons)
cox_df_z[, final_gene_sig] <- scale(cox_df_z[, final_gene_sig])

# ==============================================================================
# 9. MULTIVARIATE COX REGRESSION
# ==============================================================================

## Multivariate cox regression.
# Construct model formula containing all prognostic genes: Surv(time, event) ~ gene1 + gene2 + ...
cox_formula <- as.formula(
  paste(
    "Surv(PFS_MONTHS, PFS_STATUS) ~",
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

## Correlation matrix to depict 

# Note: The risk score prediction can be extracted for clinical stratification (low/high-risk patients)
# using the linear predictor of the fitted model:
# Code for prediction for each patient
# cox_df$risk_score <- predict(
#   fit_multi,
#   type = "lp"
# )