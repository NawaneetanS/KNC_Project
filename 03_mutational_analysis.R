#!/usr/bin/env Rscript

# 03_mutational_analysis.R
# Performs mutational analyses, including prevalence bar plots, co-mutations, and lollipop plots.

suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(ggpubr)
  library(janitor)
})

setwd("/mnt/Linux_storage/KNC")

message("--- Loading MAF datasets ---")
china_mut <- readRDS("Tables/china_maf_luad.rds")
sg_mut    <- readRDS("Tables/sg_maf_luad.rds")
msk_mut   <- readRDS("Tables/msk_maf_luad.rds")
tcga_mut  <- readRDS("Tables/tcga_maf_luad.rds")

# KNC percentage analysis
cohorts <- list(
  MSK_LUAD = msk_mut@data,
  China_LUAD = china_mut@data,
  TCGA_LUAD = tcga_mut@data,
  Singapore_LUAD = sg_mut@data
) 

knc_percent <- function(file) {
  total_mutations <- n_distinct(file$Tumor_Sample_Barcode)
  knc_mutations <- n_distinct(file$Tumor_Sample_Barcode[file$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")])
  data.frame(
    total_mutations = total_mutations,
    knc_mutations = knc_mutations,
    knc_percent = 100 * knc_mutations / total_mutations
  )
}

knc_dist <- map_dfr(cohorts, knc_percent, .id = "Cohorts")

# ==============================================
#            KNC Stacked Bar Plot
# ==============================================
plot_data <- knc_dist %>%
  mutate(
    Cohorts = factor(Cohorts, levels = c("China_LUAD", "Singapore_LUAD", "TCGA_LUAD", "MSK_LUAD")),
    Non_KNC = 100 - knc_percent,
    KNC = knc_percent
  ) %>%
  pivot_longer(
    cols = c(Non_KNC, KNC),
    names_to = "Patient_Group",
    values_to = "Percentage"
  ) %>%
  mutate(
    Patient_Group = factor(Patient_Group, levels = c("KNC", "Non_KNC"))
  )

png(filename = "Plots/Mutational/KNC_percentage.png", width = 12, height = 8, units = "in", res = 600)

ggplot(plot_data, aes(x = Cohorts, y = Percentage, fill = Patient_Group)) +
  geom_col(width = 0.4, alpha = 0.95, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(Patient_Group == "KNC", sprintf("%.1f%%\n(%d/%d)", Percentage, knc_mutations, total_mutations), "")),
    position = position_stack(vjust = 1.0), 
    vjust = -0.5,                          
    color = "black",                       
    fontface = "bold",
    size = 5
  ) +
  scale_fill_manual(
    values = c("Non_KNC" = "#E2E8F0", "KNC" = "#2C7FB8"),
    labels = c("Non_KNC" = "Other Patients", "KNC" = "KNC Patients")
  ) +
  scale_y_continuous(
    limits = c(0, 115), 
    expand = c(0, 0),
    breaks = seq(0, 115, by = 10),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "KNC Relative Prevalence Across Cohorts",
    x = NULL,
    y = "Percentage of Patients",
    fill = NULL
  ) +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_pubr(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 15), hjust = 0.5),
    axis.title.y = element_text(face = "bold", margin = margin(r = 15), color = "#444444"),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.6),
    axis.line.x = element_line(color = "black")
  )

dev.off()

# ==============================================
#            KNC co-mutational analysis
# ==============================================
message("--- Running Somatic Interactions ---")

png("Plots/Mutational/MSK_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
msk_com <- somaticInteractions(msk_mut, topMar = 5.5, leftMar = 5.5)
dev.off()
write.csv(msk_com, "Tables/MSK_co-mutation.csv")

png("Plots/Mutational/TCGA_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
tcga_com <- somaticInteractions(tcga_mut, topMar = 5.5, leftMar = 5.5)
dev.off()
write.csv(tcga_com, "Tables/TCGA_co-mutation.csv")

png("Plots/Mutational/China_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
china_com <- somaticInteractions(china_mut, topMar = 5.5, leftMar = 5.5)
dev.off()
write.csv(china_com, "Tables/China_co-mutation.csv")

png("Plots/Mutational/Singapore_co-mutations.png", width = 13, height = 9, units = "in", res = 600)
sg_com <- somaticInteractions(sg_mut, topMar = 5.5, leftMar = 5.5)
dev.off()
write.csv(sg_com, "Tables/Singapore_co-mutation.csv")

# Mutation matrix combinations
mutationMatrix <- function(maf) {
  maf@data %>%
    distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = Tumor_Sample_Barcode, values_from = value, values_fill = 0) %>% 
    t() %>% 
    row_to_names(row_number = 1)
}

msk_mm <- as.data.frame(mutationMatrix(maf = msk_mut))
write.csv(msk_mm, "Tables/MSK_mutation_matrix.csv")

# ==============================================
#            Lollipop plots for KNC genes
# ==============================================
message("--- Generating Lollipop Plots ---")

plot_lollipop <- function(maf, gene, cohort_name, outdir = "Plots", top_n_labels = 5, width = 15, height = 5){
  gene_pos <- maf@data %>%
    filter(Hugo_Symbol == gene) %>%
    mutate(pos = as.numeric(stringr::str_extract(HGVSp_Short, "\\d+"))) %>%
    filter(!is.na(pos)) %>%
    count(pos, sort = TRUE)
  
  top_pos <- head(gene_pos$pos, top_n_labels)
  
  png(
    file.path(outdir, paste0(cohort_name, "_", gene, "_lollipop.png")),
    width = width, height = height, units = "in", res = 600
  )
  par(oma = c(0, 0, 4, 0))
  
  lollipopPlot(maf = maf, gene = gene, showMutationRate = FALSE, titleSize = c(0.001, 0.001))
  
  mtext(
    paste0(gene, " Mutational Landscape in ", cohort_name, " Cohort"),
    side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2
  )
  dev.off()
  invisible(top_pos)
}

lapply(c("KEAP1", "NFE2L2", "CUL3"), function(g) plot_lollipop(msk_mut, g, "MSK", "Plots/Mutational"))
lapply(c("KEAP1", "NFE2L2", "CUL3"), function(g) plot_lollipop(tcga_mut, g, "TCGA", "Plots/Mutational"))
lapply(c("KEAP1", "NFE2L2", "CUL3"), function(g) plot_lollipop(sg_mut, g, "Singapore", "Plots/Mutational"))
lapply(c("KEAP1", "NFE2L2", "CUL3"), function(g) plot_lollipop(china_mut, g, "China", "Plots/Mutational"))

message("Mutational analysis finished successfully!")
