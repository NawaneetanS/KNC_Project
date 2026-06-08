# KNC Project: Genomic Analysis of KNC Mutations

This project involves the genomic analysis of KNC mutations across various cancer cohorts, including MSK-IMPACT, TCGA-LUAD, and others. The analysis focuses on mutation prevalence, survival outcomes, and clinical correlations.

## Project Structure

- `Analysis.R`: The main R script containing the data cleaning, analysis, and visualization logic.
- `Plots/`: Directory containing generated plots and visualizations.
  - `KNC_percentage.png`: Visualization of KNC mutation percentages across different cohorts.
  - `MSK_KNC_mutVSwt_survival.png`: Survival analysis comparing KNC mutant vs. wild-type cases.
- `public_data/`: (Not tracked) Directory containing the raw datasets used for analysis (China Pan-Cancer, MSK-IMPACT, Singapore LUAD, and TCGA).

## Datasets

The analysis utilizes the following datasets (stored locally in `public_data/`):
- China Pan-Cancer 2020
- MSK-IMPACT 50k 2026
- Singapore LUAD 2020
- TCGA-LUAD

## Getting Started

### Prerequisites

The analysis is performed in R. The following libraries are required:
- `dplyr`, `maftools`, `biomaRt`, `janitor`, `survminer`, `survival`, `data.table`, `purrr`, `tidyr`, `ggplot2`, `scales`, `ggpubr`, `patchwork`, `ggsci`

### Usage

1. Clone the repository.
2. Ensure the raw datasets are placed in the `public_data/` directory (not included in the repo).
3. Run `Analysis.R`.

## License
Refer to the licenses provided with the individual datasets in the `public_data/` directory.
