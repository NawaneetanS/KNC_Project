#!/usr/bin/env Rscript

# 04_survival_analysis.R
# Conducts overall survival and variant-level survival analysis by cohort.

suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(survival)
  library(survminer)
  library(ggpubr)
})

setwd("/mnt/Linux_storage/KNC")

message("--- Loading MAF datasets ---")
china_mut <- readRDS("Tables/china_maf_luad.rds")
sg_mut    <- readRDS("Tables/sg_maf_luad.rds")
msk_mut   <- readRDS("Tables/msk_maf_luad.rds")
tcga_mut  <- readRDS("Tables/tcga_maf_luad.rds")

# ==============================================
#            Overall KNC survival functions
# ==============================================
knc_survival_plot <- function(maf, time_col, status_col, cohort_name, outfile = NULL) {
  clin <- maf@clinical.data
  
  # Create KNC status
  knc_pax <- unique(maf@data$Tumor_Sample_Barcode[maf@data$Hugo_Symbol %in% c("KEAP1", "NFE2L2", "CUL3")])
  
  clin$KNC_status <- ifelse(
    clin$Tumor_Sample_Barcode %in% knc_pax,
    "KNC Mutated",
    "KNC Wild-Type"
  )
  
  # Convert survival status
  if (is.character(clin[[status_col]])) {
    clin[[status_col]] <- ifelse(
      clin[[status_col]] %in% c("DECEASED", "Dead", "DEAD", "1", "1:DECEASED"),
      1, 0
    )
  }
  
  # Clean survival data
  clin <- clin %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[status_col]]),
      .data[[status_col]] != ""
    )
  
  clin[[time_col]] <- as.numeric(clin[[time_col]])
  clin[[status_col]] <- as.numeric(clin[[status_col]])
  
  form_str <- paste0("Surv(", time_col, ", ", status_col, ") ~ KNC_status")
  surv_form <- as.formula(form_str)
  
  fit <- survfit(surv_form, data = clin)
  fit$call$formula <- surv_form
  
  # Log-rank test
  pval <- survdiff(surv_form, data = clin)
  
  # Plot
  surv_plot <- ggsurvplot(
    fit,
    data = clin,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    title = paste("Overall Survival by KNC Mutation Status in", cohort_name, "cohort"),
    font.title = c(16, "bold"),
    xlab = "Overall Survival (Months)",
    ylab = "Survival Probability",
    legend.title = NULL,
    legend.labs = c("KNC Mutated", "KNC Wild-Type"),
    legend = "top",
    palette = c("#2C7FB8", "#D95F02"),
    surv.median.line = "hv",
    censor = TRUE,
    censor.shape = 124,
    censor.size = 3,
    ggtheme = theme_pubr(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)),
    conf.int.alpha = 0.15
  )
  
  if (!is.null(outfile)) {
    png(outfile, width = 12, height = 9, units = "in", res = 600)
    print(surv_plot)
    dev.off()
  }
  
  return(list(fit = fit, survdiff = pval, plot = surv_plot, data = clin))
}

message("--- Running Overall KNC Survival ---")
msk_os  <- knc_survival_plot(msk_mut, "OS_MONTHS", "OS_STATUS", "MSK", "Plots/Survival/MSK_KNC_survival.png")
sg_os   <- knc_survival_plot(sg_mut,  "OS_MONTHS", "OS_STATUS", "Singapore", "Plots/Survival/SG_KNC_survival.png")
tcga_os <- knc_survival_plot(tcga_mut, "OS_MONTHS", "OS_STATUS", "TCGA", "Plots/Survival/TCGA_KNC_survival.png")

# ==============================================
#            KNC variant survival
# ==============================================
message("--- Running Recurrent Variant Survival ---")

knc_variant_survival <- function(maf, genes, cohort_name, time_col = "OS_MONTHS", status_col = "OS_STATUS", outdir = "Plots", top_n = 10) {
  maf_sub <- subsetMaf(maf = maf, genes = genes)
  
  variant_freq <- maf_sub@data %>%
    filter(!is.na(HGVSp_Short)) %>%
    count(Hugo_Symbol, HGVSp_Short, name = "Frequency", sort = TRUE)
  
  top_variants <- variant_freq %>%
    filter(Frequency > 1) %>%
    slice_max(Frequency, n = top_n, with_ties = TRUE)
  
  if (nrow(top_variants) == 0) {
    message(cohort_name, ": No recurrent variants found.")
    return(NULL)
  }
  
  top_patients <- unique(maf_sub@data$Tumor_Sample_Barcode[maf_sub@data$HGVSp_Short %in% top_variants$HGVSp_Short])
  
  surv_df <- maf_sub@clinical.data
  surv_df$var_prevalence <- ifelse(
    surv_df$Tumor_Sample_Barcode %in% top_patients,
    "Top recurrent",
    "Other variants"
  )
  
  cat("\n")
  cat("Cohort:", cohort_name, "\n")
  cat("Before filtering:", nrow(surv_df), "\n")
  print(summary(surv_df[[time_col]]))
  print(table(surv_df[[status_col]], useNA = "always"))
  
  surv_df <- surv_df[
    !is.na(surv_df[[time_col]]) &
      !is.na(surv_df[[status_col]]) &
      surv_df[[status_col]] != "",
  ]
  
  cat("After filtering:", nrow(surv_df), "\n")
  
  surv_df$TIME <- as.numeric(surv_df[[time_col]])
  surv_df$STATUS <- ifelse(
    surv_df[[status_col]] %in% c("DECEASED", "Dead", "dead", 1),
    1, 0
  )
  
  surv_df$var_prevalence <- factor(surv_df$var_prevalence, levels = c("Other variants", "Top recurrent"))
  
  cat("\n")
  cat(cohort_name, "\n")
  print(table(surv_df$var_prevalence))
  cat("\n")
  
  if (length(unique(surv_df$var_prevalence)) < 2) {
    message(cohort_name, ": only one group present.")
    return(NULL)
  }
  
  fit <- survfit(Surv(TIME, STATUS) ~ var_prevalence, data = surv_df)
  
  p <- ggsurvplot(
    fit, data = surv_df, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
    risk.table = TRUE, risk.table.height = 0.25,
    title = paste0(cohort_name, ": Survival by KNC Variant Recurrence"),
    legend.title = NULL,
    legend.labs = c("Other variants", "Top recurrent variants"),
    ggtheme = theme_pubr(base_size = 14)
  )
  
  png(
    file.path(outdir, paste0(cohort_name, "_variant_survival.png")),
    width = 8, height = 6, units = "in", res = 600
  )
  print(p)
  dev.off()
  
  return(list(
    maf_subset = maf_sub, variant_frequency = variant_freq,
    top_variants = top_variants, top_patients = top_patients,
    survival_data = surv_df, fit = fit, plot = p
  ))
}

cohorts_list <- list(
  MSK = msk_mut,
  China = china_mut,
  Singapore = sg_mut,
  TCGA = tcga_mut
)

results <- lapply(
  names(cohorts_list),
  function(x)
    knc_variant_survival(
      maf = cohorts_list[[x]],
      genes = c("KEAP1", "NFE2L2", "CUL3"),
      cohort_name = x,
      outdir = "Plots/Survival"
    )
)

names(results) <- names(cohorts_list)
message("Survival analysis finished successfully!")
