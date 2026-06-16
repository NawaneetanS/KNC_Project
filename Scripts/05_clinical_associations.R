#!/usr/bin/env Rscript

# 05_clinical_associations.R
# Calculates statistical associations (Fisher's exact test) between KNC pathway mutations and clinical characteristics.

suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(purrr)
  library(tidyr)
})

setwd("/mnt/Linux_storage/KNC")

# Define genes in KNC pathway (KEAP1, NFE2L2, CUL3)
knc_genes <- c("KEAP1", "NFE2L2", "CUL3")

run_association <- function(cohort_name, clin_file, maf_file) {
  # Load cleaned clinical df
  clin <- readRDS(clin_file)
  
  # Load filtered MAF
  maf <- readRDS(maf_file)
  
  # Find barcodes with KNC mutations
  knc_mut_barcodes <- unique(
    maf@data$Tumor_Sample_Barcode[maf@data$Hugo_Symbol %in% knc_genes]
  )
  
  # Annotate clinical data with mutation status
  clin$knc_mut <- ifelse(
    clin$Tumor_Sample_Barcode %in% knc_mut_barcodes,
    "Mutated",
    "Wild-Type"
  )
  clin$knc_mut <- factor(clin$knc_mut, levels = c("Mutated", "Wild-Type"))
  
  # List to store results
  results_list <- list()
  
  # Overall KNC mutation prevalence
  total_mut <- sum(clin$knc_mut == "Mutated")
  total_pts <- nrow(clin)
  results_list[[length(results_list) + 1]] <- data.frame(
    Characteristic = "All cases",
    Subgroup = "All",
    Mutated_Total = paste0(total_mut, "/", total_pts, " (", round(100 * total_mut / total_pts, 1), "%)"),
    p_value = ""
  )
  
  # Helper function to run Fisher's exact test and format results
  test_char <- function(var_name, var_label) {
    if (!var_name %in% colnames(clin)) {
      return(NULL)
    }
    
    # Filter out missing/unknown/NA categories for the association test
    temp_df <- clin %>% 
      filter(!is.na(.data[[var_name]]), as.character(.data[[var_name]]) != "", as.character(.data[[var_name]]) != "Unknown", as.character(.data[[var_name]]) != "NA")
    
    if (nrow(temp_df) == 0) return(NULL)
    
    # Create contingency table
    tbl <- table(temp_df[[var_name]], temp_df$knc_mut)
    
    # If dimensions are not correct, skip
    if (nrow(tbl) < 2) return(NULL)
    
    # Two-sided Fisher's Exact Test
    ft <- fisher.test(tbl)
    p_val <- ft$p.value
    p_str <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
    p_str <- ifelse(p_val >= 0.05, "NS", p_str)
    
    subgroups <- rownames(tbl)
    rows <- list()
    for (i in 1:length(subgroups)) {
      sub <- subgroups[i]
      mut <- tbl[i, "Mutated"]
      wt <- tbl[i, "Wild-Type"]
      tot <- mut + wt
      pct <- ifelse(tot > 0, round(100 * mut / tot, 1), 0)
      
      rows[[i]] <- data.frame(
        Characteristic = ifelse(i == 1, var_label, ""),
        Subgroup = sub,
        Mutated_Total = paste0(mut, "/", tot, " (", pct, "%)"),
        p_value = ifelse(i == 1, p_str, "")
      )
    }
    do.call(rbind, rows)
  }
  
  # 1. Sex
  sex_res <- test_char("SEX", "Sex")
  if (!is.null(sex_res)) results_list[[length(results_list) + 1]] <- sex_res
  
  # 2. Age (< 65 vs >= 65)
  if ("AGE" %in% colnames(clin)) {
    clin$Age_Group <- ifelse(clin$AGE < 65, "< 65", ">= 65")
    clin$Age_Group <- factor(clin$Age_Group, levels = c("< 65", ">= 65"))
    age_res <- test_char("Age_Group", "Age (years)")
    if (!is.null(age_res)) results_list[[length(results_list) + 1]] <- age_res
  }
  
  # 3. AJCC Stage (Grouped Stage I-II vs III-IV)
  if ("STAGE" %in% colnames(clin)) {
    # Individual stage
    stage_res <- test_char("STAGE", "AJCC Stage (Individual)")
    if (!is.null(stage_res)) results_list[[length(results_list) + 1]] <- stage_res
    
    # Grouped stage
    clin$Stage_Grouped <- case_when(
      clin$STAGE %in% c("I", "II", "0", "I-II") ~ "Stage I-II",
      clin$STAGE %in% c("III", "IV", "III-IV") ~ "Stage III-IV",
      TRUE ~ NA_character_
    )
    stage_grp_res <- test_char("Stage_Grouped", "AJCC Stage (Grouped)")
    if (!is.null(stage_grp_res)) results_list[[length(results_list) + 1]] <- stage_grp_res
  }
  
  # 4. Smoking Status
  if ("SMOKING_STATUS" %in% colnames(clin)) {
    smoke_res <- test_char("SMOKING_STATUS", "Smoking History")
    if (!is.null(smoke_res)) results_list[[length(results_list) + 1]] <- smoke_res
  } else if ("SMOKE_STATUS" %in% colnames(clin)) {
    smoke_res <- test_char("SMOKE_STATUS", "Smoking History")
    if (!is.null(smoke_res)) results_list[[length(smoke_res) + 1]] <- smoke_res
  }
  
  # 5. Sample Type
  if ("Sample_Type" %in% colnames(clin)) {
    st_res <- test_char("Sample_Type", "Sample Type")
    if (!is.null(st_res)) results_list[[length(results_list) + 1]] <- st_res
  } else if ("SAMPLE_TYPE" %in% colnames(clin)) {
    st_res <- test_char("SAMPLE_TYPE", "Sample Type")
    if (!is.null(st_res)) results_list[[length(st_res) + 1]] <- st_res
  }
  
  # Combine results
  final_df <- do.call(rbind, results_list)
  
  # Save to tables folder (ignored by git)
  write.csv(final_df, paste0("Tables/", cohort_name, "_knc_associations.csv"), row.names = FALSE)
  
  invisible(final_df)
}

# Run for all four cohorts
message("--- Calculating KNC clinical associations ---")
run_association("Singapore", "Tables/sg_clean.rds", "Tables/sg_maf_luad.rds")
run_association("China", "Tables/china_clean.rds", "Tables/china_maf_luad.rds")
run_association("TCGA", "Tables/tcga_clean.rds", "Tables/tcga_maf_luad.rds")
run_association("MSK", "Tables/msk_clean.rds", "Tables/msk_maf_luad.rds")

message("Clinical associations analysis completed successfully!")
