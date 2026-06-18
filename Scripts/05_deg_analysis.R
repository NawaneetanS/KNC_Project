#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyverse)
  library(AnnotationDbi)
  library(EnhancedVolcano)
  library(biomaRt)
})

setwd("/mnt/Linux_storage/KNC")

## Read the TCGA rds file and take knc patients tsb
tcgaRDS <- readRDS("Tables/tcga_maf_luad.rds")

knc_pax <- tcgaRDS@data$Tumor_Sample_Barcode[tcgaRDS@data$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")]

## Read the RNA seq data
raw_counts <- data.table::fread("public_data/TCGA/luad_tcga_gdc/data_mrna_seq_read_counts.txt", sep = "\t", nThread = 8)

# Convert entrez to hugo
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- raw_counts$Entrez_Gene_Id

# Fetch HUGO symbols using getBM
hugo <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "hgnc_symbol"),
  values = genes,
  mart = mart
)

# Join the column of hgnc with the raw counts
raw_counts <- raw_counts %>% 
  left_join(
    hugo,
    by = c("Entrez_Gene_Id" = "entrezgene_id")
  ) %>% 
  dplyr::select(hgnc_symbol, everything(), -Entrez_Gene_Id) %>% 
  dplyr::filter(!is.na(hgnc_symbol),
                hgnc_symbol != "",
                rowSums(across(2:ncol(.))) > 0) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(across(everything(), sum), 
            .groups = "drop" ) %>% 
  column_to_rownames(var = "hgnc_symbol")

# Prepare files for edgeR
sample_info <- data.frame(
  Sample = colnames(raw_counts)
)

sample_info$Group <- ifelse(
  sample_info$Sample %in% knc_pax,
  "KNC",
  "Non-KNC"
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
write.csv(sig_deg, "Tables/Significant_genes.csv", sep = ",")
