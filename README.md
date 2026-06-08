# KNC Project: Genomic and Clinical Analysis of KNC Mutations

This project investigates the prevalence and clinical significance of mutations in the **KNC** gene triad—**KEAP1**, **NFE2L2**, and **CUL3**—across diverse lung adenocarcinoma (LUAD) cohorts. These genes are central components of the Nrf2 signaling pathway, which plays a critical role in cellular antioxidant defense and is frequently dysregulated in cancer.

## Research Objectives

The primary goal of this analysis is to characterize the landscape of KNC mutations and their impact on patient outcomes across global populations. Key objectives include:

- **Prevalence Mapping:** Quantifying the mutation frequency of KEAP1, NFE2L2, and CUL3 across four major lung cancer cohorts:
  - **China LUAD** (China Pan-Cancer 2020)
  - **Singapore LUAD** (Singapore LUAD 2020)
  - **TCGA-LUAD** (The Cancer Genome Atlas)
  - **MSK-LUAD** (MSK-IMPACT 50k 2026)
- **Cohort Comparison:** Identifying variations in KNC mutation relative prevalence across different geographical and ethnic populations.
- **Survival Analysis:** Evaluating the impact of KNC mutation status on overall survival (OS), specifically within the large-scale MSK-IMPACT cohort, to determine if these mutations correlate with poorer clinical outcomes.

## Key Findings (Visualized)

The analysis generates publication-quality visualizations to summarize these findings:

1.  **KNC Relative Prevalence (`Plots/KNC_percentage.png`):** A comparative analysis showing the percentage of patients harboring KNC mutations in each cohort. This visualization highlights how KNC pathway alterations vary by cohort, providing insights into population-specific genomic landscapes.
2.  **Survival Impact (`Plots/MSK_KNC_mutVSwt_survival.png`):** Kaplan-Meier survival curves comparing KNC-mutated versus wild-type patients in the MSK cohort. The plot includes risk tables and p-values to assess the statistical significance of survival differences.

## Analytical Workflow

The analysis is implemented in R (`Analysis.R`) and follows a rigorous pipeline:
1.  **Data Harmonization:** Cleaning and merging clinical sample and patient data from multiple sources.
2.  **Genomic Integration:** Processing Mutation Annotation Format (MAF) and TSV files to identify non-synonymous mutations in KNC genes.
3.  **Statistical Modeling:** Calculating prevalence rates and performing survival analysis using the `survival` and `survminer` packages.
4.  **Visualization:** Generating professional-grade plots using `ggplot2` and `ggpubr`.

---

*Note: This repository contains the analytical code and resulting visualizations. Raw genomic and clinical datasets are stored locally in the `public_data/` directory and are not tracked in this repository.*
