#!/usr/bin/Rscript

# run_pipeline.R
# Orchestrates the execution of KNC cohort analysis steps.

setwd("/mnt/Linux_storage/KNC")

args <- commandArgs(trailingOnly = TRUE)

run_step <- function(script_name, desc) {
  message("\n==================================================")
  message(" Running: ", desc, " (", script_name, ")")
  message("==================================================")
  
  # Execute the script in a separate R session
  status <- system2("Rscript", script_name)
  if (status != 0) {
    stop("Error: ", script_name, " failed with status ", status)
  }
}

if (length(args) == 0) {
  # Run entire pipeline
  run_step("Scripts/01_load_clean_data.R", "Step 1: Data Preprocessing & RDS Export")
  run_step("Scripts/02_clinical_dashboards.R", "Step 2: Clinical & Pathological Dashboards")
  run_step("Scripts/03_mutational_analysis.R", "Step 3: Mutational Profiling & Lollipop Plots")
  run_step("Scripts/04_survival_analysis.R", "Step 4: Overall & Variant Survival Analysis")
  message("\nPipeline completed successfully!")
} else {
  step <- args[1]
  if (step == "1" || step == "clean") {
    run_step("Scripts/01_load_clean_data.R", "Step 1: Data Preprocessing & RDS Export")
  } else if (step == "2" || step == "dashboard") {
    run_step("Scripts/02_clinical_dashboards.R", "Step 2: Clinical & Pathological Dashboards")
  } else if (step == "3" || step == "mutation") {
    run_step("Scripts/03_mutational_analysis.R", "Step 3: Mutational Profiling & Lollipop Plots")
  } else if (step == "4" || step == "survival") {
    run_step("Scripts/04_survival_analysis.R", "Step 4: Overall & Variant Survival Analysis")
  } else if (step == "5" || step == "deg") {
    run_step("Scripts/05_deg_analysis.R", "Step 5: DEG analysis and cox regressions")
  }
    else {
    message("Invalid argument. Usage:")
    message("  Rscript Scripts/run_pipeline.R             (Runs entire pipeline)")
    message("  Rscript Scripts/run_pipeline.R clean       (Runs data cleaning/preprocessing only)")
    message("  Rscript Scripts/run_pipeline.R dashboard   (Runs clinical dashboards only)")
    message("  Rscript Scripts/run_pipeline.R mutation    (Runs mutational analysis only)")
    message("  Rscript Scripts/run_pipeline.R survival    (Runs survival analysis only)")
  }
}
