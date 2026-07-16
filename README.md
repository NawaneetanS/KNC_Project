# KNC Project: Genomic, Clinical, and Transcriptomic Analysis of KNC Mutations in Lung Adenocarcinoma

This repository contains the complete analytical pipeline and visual outputs for investigating the prevalence, clinical significance, and transcriptomic impact of mutations in the **KNC** gene triad—**KEAP1**, **NFE2L2**, and **CUL3**—across diverse Lung Adenocarcinoma (LUAD) cohorts. These genes regulate the Nrf2 signaling pathway, which controls cellular antioxidant responses and is frequently dysregulated in cancer, leading to therapy resistance and poor clinical outcomes.

---

## Research Objectives
1. **Prevalence Mapping:** Quantify the mutational frequency of *KEAP1*, *NFE2L2*, and *CUL3* across four major cohorts:
   - **China LUAD** (China Pan-Cancer 2020)
   - **Singapore LUAD** (Singapore LUAD 2020)
   - **TCGA-LUAD** (The Cancer Genome Atlas)
   - **MSK-LUAD** (MSK-IMPACT 50k 2026)
2. **Clinical Profiling:** Analyze cohort demographics, staging, smoking status, tumor purity, and Tumor Mutational Burden (TMB).
3. **Survival Analysis:** Estimate the impact of KNC mutation status and specific recurrent variants on Overall Survival (OS).
4. **Differential Expression & Prognostic Signature:** Identify Differentially Expressed Genes (DEGs) and proteins associated with KNC mutations to build a prognostic gene signature.
5. **Spatial QC:** Perform quality control and preprocessing on Visium spatial transcriptomics datasets.

---

## Project Structure

```
KNC/
├── Articles/                  # Relevant literature and reports
│   ├── download.pdf
│   └── knc_gene-sig.pdf
├── Plots/                     # Generated publication-quality figures
│   ├── Clinical/              # Cohort demographic & pathological dashboards
│   ├── Genes/                 # Gene correlation plots
│   ├── Mutational/            # Prevalence histograms, co-mutation heatmaps, lollipop plots
│   ├── Spatial/               # Spatial transcriptomics QC plots
│   ├── Survival/              # Kaplan-Meier survival curves
│   └── TCGA_KNC_volcano.png   # Volcano plot for TCGA differential expression
├── Scripts/                   # R and Python scripts for the analysis
│   ├── 01_load_clean_data.R   # Raw data cleaning and object exporting
│   ├── 02_clinical_dashboards.R # Dashboard visualization code
│   ├── 03_mutational_analysis.R # Lollipop, bar plots, co-mutations
│   ├── 04_survival_analysis.R # Kaplan-Meier survival & variant analysis
│   ├── 05_deg_analysis.R      # edgeR/DESeq2 and Cox hazard modeling
│   ├── 06_spatial_analysis.ipynb # Spatial transcriptomics QC pipeline
│   ├── run_pipeline.R         # Master pipeline orchestrator script
│   └── deg_cox_analysis.txt   # Text log output of the DEG/Cox analysis
├── Tables/                    # Preprocessed datasets and statistical summaries
├── environment.yml            # Conda environment definition
├── requirements.txt           # Python packages and requirements
└── README.md                  # This file
```

---

## Analytical Workflow

The pipeline is organized into modular scripts under the [Scripts/](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts) directory:

1. **Step 1: Data Preprocessing & Formatting**
   - **Script:** [01_load_clean_data.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/01_load_clean_data.R)
   - **Details:** Harmonizes demographic metrics (age, stage, sex, smoking history, treatment profile) across cohorts and processes Mutation Annotation Format (MAF) files, saving cleaned objects into the `Tables/` directory.

2. **Step 2: Demographic & Pathological Dashboards**
   - **Script:** [02_clinical_dashboards.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/02_clinical_dashboards.R)
   - **Details:** Compiles comprehensive dashboards visualizing age distributions, stage, gender, smoking status, and molecular features like TMB.
   - **Outputs:**
     - [China clinical dashboard](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Clinical/China_clinical_dashboard.png)
     - [MSK clinical dashboard](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Clinical/MSK_clinical_dashboard.png)
     - [Singapore clinical dashboard](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Clinical/Singapore_clinical_dashboard.png)
     - [TCGA clinical dashboard](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Clinical/TCGA_clinical_dashboard.png)

3. **Step 3: Mutational Profiling & Lollipop Plots**
   - **Script:** [03_mutational_analysis.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/03_mutational_analysis.R)
   - **Details:** Quantifies relative prevalence, maps mutation spots along protein domains (lollipop plots), and executes somatic interaction analysis to identify co-mutation or mutual exclusivity trends.
   - **Outputs:**
     - [KNC relative prevalence histogram](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Mutational/KNC_percentage.png)
     - [MSK co-mutations heatmap](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Mutational/MSK_co-mutations.png)
     - [TCGA co-mutations heatmap](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Mutational/TCGA_co-mutations.png)
     - Lollipop plots for *KEAP1*, *NFE2L2*, and *CUL3* across all cohorts (saved in [Plots/Mutational/](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Mutational/))

4. **Step 4: Survival Analysis**
   - **Script:** [04_survival_analysis.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/04_survival_analysis.R)
   - **Details:** Generates Kaplan-Meier curves comparing Overall Survival (OS) between KNC-mutated and wild-type patients. It also looks at outcomes for patients harboring top recurrent variants.
   - **Outputs:**
     - [MSK cohort KNC survival KM plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Survival/MSK_KNC_survival.png)
     - [TCGA cohort KNC survival KM plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Survival/TCGA_KNC_survival.png)
     - [Singapore cohort KNC survival KM plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Survival/SG_KNC_survival.png)
     - [MSK variant-level survival plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Survival/MSK_variant_survival.png)
     - [TCGA variant-level survival plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Survival/TCGA_variant_survival.png)

5. **Step 5: Differential Gene Expression & Cox Proportional Hazards Model**
   - **Script:** [05_deg_analysis.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/05_deg_analysis.R)
   - **Details:** Compares KNC-mutated vs. wild-type tumors. Performs RNA-seq DEG analysis in TCGA-LUAD (edgeR/DESeq2) and proteomics analysis in CPTAC-LUAD (limma) to identify differentially expressed genes/proteins. Fits univariate and multivariate Cox models to identify prognostic genes and outputs correlation plots.
   - **Outputs:**
     - [TCGA KNC DEGs Volcano Plot](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/TCGA_KNC_volcano.png)
     - [Univariate Cox Genes Correlation Heatmap](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Genes/Correlation_plot_uni-cox_genes.png)
     - Tables of univariate and multivariate results (saved in `Tables/`)

6. **Step 6: Spatial Transcriptomics Preprocessing & Quality Control**
   - **Script:** [06_spatial_analysis.ipynb](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/06_spatial_analysis.ipynb)
   - **Details:** Runs a python-based pipeline using `scanpy` and `squidpy` to ingest Visium spatial coordinates and H&E images, check dimensions, annotate mitochondrial genes, and filter out low-quality spots/bins.
   - **Outputs:** QC reports and distribution plots located in [Plots/Spatial/QC/](file:///media/nannu1375/Backpack/Shankara/KNC/Plots/Spatial/QC/)

---

## Execution Guide

### Environment Setup
To replicate the environment and dependencies:
- **R Environment:** Check `environment.yml` for Conda setup, which includes necessary bioinformatic packages (`maftools`, `edgeR`, `DESeq2`, `survival`, `survminer`, `EnhancedVolcano`, `pheatmap`, etc.).
- **Python Environment:** Install standard dependencies using:
  ```bash
  pip install -r requirements.txt
  ```

### Running the R Analysis Pipeline
The R analysis workflow is orchestrated via the master script [run_pipeline.R](file:///media/nannu1375/Backpack/Shankara/KNC/Scripts/run_pipeline.R). 

To run the **entire R pipeline** (Steps 1 to 5) sequentially:
```bash
Rscript Scripts/run_pipeline.R
```

To run a **specific step** of the analysis:
- **Preprocessing only:** `Rscript Scripts/run_pipeline.R clean`
- **Clinical dashboards:** `Rscript Scripts/run_pipeline.R dashboard`
- **Mutational profiling:** `Rscript Scripts/run_pipeline.R mutation`
- **Survival analysis:** `Rscript Scripts/run_pipeline.R survival`
- **DEG & Cox regression:** `Rscript Scripts/run_pipeline.R deg`

---

*Note: Raw genomic and clinical datasets from cohorts are kept locally in `public_data/` and are not tracked under version control.*

