#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyverse)
  library(AnnotationDbi)
  library(EnhancedVolcano)
  library(biomaRt)
  library(maftools)
  library(janitor)
  library(limma)
  library(DESeq2)
  library(survival)
  library(fgsea)
})

setwd("/mnt/Linux_storage/KNC")

## Read the TCGA rds file and take KNC patients
tcgaRDS <- readRDS("Tables/tcga_maf_luad.rds")

knc_pax <- tcgaRDS@data$Tumor_Sample_Barcode[
  tcgaRDS@data$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")
]

# Convert to patient IDs
knc_pax <- substr(knc_pax, 1, 12)
knc_pax <- unique(knc_pax)

## Read RNA-seq counts
raw_counts <- data.table::fread(
  "public_data/TCGA/luad_tcga_pan_can_atlas_2018/raw_counts_LUAD.csv",
  sep = ",",
  nThread = 8
) %>%
  filter(gene_type == "protein_coding") %>% 
  dplyr::select(-gene_type)

## Keep only primary tumour samples (sample type = 01)
sample_cols <- colnames(raw_counts)[-1]

tumor_cols <- sample_cols[
  substr(sample_cols, 14, 15) == "01"
]

raw_counts <- raw_counts %>%
  dplyr::select(
    Gene,
    all_of(tumor_cols)
  )

## Convert TCGA barcodes to patient IDs
colnames(raw_counts)[-1] <- substr(
  colnames(raw_counts)[-1],
  1,
  12
)

## Remove duplicate tumour aliquots
raw_counts <- raw_counts[
  ,
  !duplicated(names(raw_counts)),
  with = FALSE
]

## Collapse duplicate gene symbols
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
raw_counts <- as.data.frame(raw_counts)

rownames(raw_counts) <- raw_counts$Gene

raw_counts$Gene <- NULL

## Create sample information
sample_info <- data.frame(
  Sample = colnames(raw_counts)
)

sample_info$Group <- ifelse(
  sample_info$Sample %in% knc_pax,
  "KNC",
  "Non-KNC"
)

sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Non-KNC", "KNC")
)

sample_info$Group <- as.factor(sample_info$Group)

# Create DGE list
dge <- DGEList(
  counts = raw_counts,
  group = sample_info$Group
)

# Filter out low counts per million genes (CPM)
keep <- filterByExpr(
  dge,
  sample_info$Group
)

# Keep only high CPM genes
dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

# Normalise (TMM)
dge <- calcNormFactors(dge)

# Set the comparison groups and create a design matrix
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Non-KNC", "KNC")
)

design <- model.matrix(
  ~ Group,
  data = sample_info
)

# Estimate dispersion
dge <- estimateDisp(
  dge,
  design
)

# Fit model
fit <- glmQLFit(
  dge,
  design
)

# Differential expression
res <- glmQLFTest(
  fit,
  coef = 2
)

# Extract results
deg <- topTags(
  res,
  n = Inf
)$table

# Significant deg
sig_deg <- deg %>%
  filter(
    FDR < 0.05 &
      abs(logFC) > 1.5
  )

# Volcano
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

# CSV
write.csv(sig_deg, "Tables/Significant_genes.csv")

## Run GSEA and take top pathway genes for filtering and further cox regression
hallmark_pathway <- gmtPathways("public_data/h.all.v2026.1.Hs.symbols.gmt")

# Create genes ranked list
gene_ranks <- deg$logFC
names(gene_ranks) <- rownames(deg)

# Run GSEA
fgsea_res <- fgsea(
  pathways = hallmark_pathway,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500
)

fgsea_res <- fgsea_res %>% 
  filter(padj < 0.05)

fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

# Get list of genes in top pathways
genes_list <- unique(unlist(fgsea_res$leadingEdge))

## Filter DEGs in top_prot_limma
final_deg <- sig_deg %>%
  filter(rownames(sig_deg) %in% genes_list)

gene_sig <- rownames(final_deg)


## Run VST normalisation on the raw counts
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_info,
  design = ~ Group
)

vsd <- vst(dds, blind = FALSE)

vst_mat <- assay(vsd)

## Filter normalised data for the 14 genes
vst_14 <- vst_mat[
  rownames(vst_mat) %in% gene_sig,
]

expr_df <- as.data.frame(t(vst_14))

## Filter normalised data for the 14 genes

vst_37 <- vst_mat[
  rownames(vst_mat) %in% rownames(sig_deg), 
]

expr_df <- as.data.frame(t(vst_37))

## Combine PFS data to this matrix
pfs_df <- read.delim("public_data/TCGA/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt") %>% 
  dplyr::slice(-c(1,2,3)) %>% 
  janitor::row_to_names(row_number = 1) %>% 
  as.data.frame() %>% 
  dplyr::select(PATIENT_ID, OS_STATUS, OS_MONTHS) %>% 
  dplyr::filter(OS_MONTHS != "",
                !is.na(OS_MONTHS)) %>% 
  dplyr::mutate(OS_STATUS = ifelse(OS_STATUS == "1:DECEASED", 1, 0))

cox_df <- expr_df %>%
  rownames_to_column(var = "PATIENT_ID") %>%
  mutate(
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


## Univariate COX regression
genes <- colnames(cox_df)[
  !(colnames(cox_df) %in%
      c("PATIENT_ID", "OS_STATUS", "OS_MONTHS"))
]

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
    
    data.frame(
      Gene = gene,
      HR = s$coef[,"exp(coef)"],
      Lower95 = s$conf.int[,"lower .95"],
      Upper95 = s$conf.int[,"upper .95"],
      PValue = s$coef[,"Pr(>|z|)"]
    )
  }
)

cox_results <- do.call(rbind, cox_results)

cox_results$FDR <- p.adjust(
  cox_results$PValue,
  method = "BH"
)

cox_results[order(cox_results$PValue), ]

cox_results <- cox_results %>% 
  dplyr::filter(PValue < 0.05)

## Get final gene signature
final_gene_sig <- cox_results$Gene

## Convert univariate cox results to z score
cox_df_z <- cox_df[, c("PATIENT_ID", final_gene_sig, "OS_STATUS", "OS_MONTHS")]

cox_df_z[, final_gene_sig] <- scale(cox_df_z[, final_gene_sig])

## Multivariate cox regression.
cox_formula <- as.formula(
  paste(
    "Surv(OS_MONTHS, OS_STATUS) ~",
    paste(
      paste0("`", final_gene_sig, "`"),
      collapse = " + "
    )
  )
)

fit_multi <- coxph(
  formula = cox_formula,
  data = cox_df_z
)

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

multi_results <- multi_results %>% 
  filter(PValue < 0.05)

## Correlation matrix to depict 

# Code for prediction for each patient
# cox_df$risk_score <- predict(
#   fit_multi,
#   type = "lp"
# )