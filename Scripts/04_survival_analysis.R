#!/usr/bin/Rscript

# ==============================================================================
# PIPELINE: CPTAC pNRF2 vs. non-pNRF2 Overall Survival Analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

setwd("/media/nannu1375/Backpack/Shankara/KNC")

message("--- Loading CPTAC clinical and phosphoprotein datasets ---")
cptac_clean <- readRDS("Tables/cptac_clean.rds")
cptac_phospho <- readRDS("Tables/cptac_phospho.rds")

# Extract the NFE2L2 (NRF2) phosphosites
nrf2_phospho <- cptac_phospho %>%
  filter(NAME %in% c("NFE2L2_S215s", "NFE2L2_S433s")) %>%
  as.data.frame()

# Extract patient columns (excluding metadata columns 1 to 4)
patient_cols <- colnames(nrf2_phospho)[-c(1,2,3,4)]

# Determine if each patient has positive average pNRF2 expression (> 0)
val_s215 <- as.numeric(nrf2_phospho[nrf2_phospho$NAME == "NFE2L2_S215s", patient_cols])
val_s433 <- as.numeric(nrf2_phospho[nrf2_phospho$NAME == "NFE2L2_S433s", patient_cols])

# Compute average of both sites, treating NAs (undetected) as below average (<= 0)
avg_nrf2 <- colMeans(rbind(val_s215, val_s433), na.rm = TRUE)
avg_nrf2[is.na(avg_nrf2)] <- -Inf  # Map patients with NA for both sites to negative group

# Set up grouping dataframe
group_df <- data.frame(
  PATIENT_ID = patient_cols,
  pNRF2_Status = ifelse(avg_nrf2 > 0, "pNRF2", "non pNRF2"),
  stringsAsFactors = FALSE
)

# Merge survival data from cptac_clean
surv_df <- cptac_clean %>%
  dplyr::select(PATIENT_ID, OS_STATUS, OS_MONTHS) %>%
  dplyr::inner_join(group_df, by = "PATIENT_ID") %>%
  dplyr::filter(!is.na(OS_MONTHS), OS_MONTHS != "", !is.na(OS_STATUS), OS_STATUS != "") %>%
  dplyr::mutate(
    TIME = as.numeric(OS_MONTHS),
    STATUS = ifelse(grepl("DECEASED", OS_STATUS, ignore.case = TRUE), 1, 0)
  )

# Convert to factor with "non pNRF2" as the reference level
surv_df$pNRF2_Status <- factor(surv_df$pNRF2_Status, levels = c("non pNRF2", "pNRF2"))

cat("pNRF2 Group Counts:\n")
print(table(surv_df$pNRF2_Status))

# Check if we have at least 2 groups with patients to plot survival curves
if (length(unique(surv_df$pNRF2_Status)) < 2) {
  stop("Cannot perform survival analysis: less than 2 groups present in the survival dataset.")
}

message("--- Fitting survival curve for pNRF2 vs. non-pNRF2 ---")
fit <- survfit(Surv(TIME, STATUS) ~ pNRF2_Status, data = surv_df)

# Fit Cox Proportional Hazards model to get Hazard Ratio
fit_cox <- coxph(Surv(TIME, STATUS) ~ pNRF2_Status, data = surv_df)
sum_cox <- summary(fit_cox)
hr <- sum_cox$conf.int[1]
hr_lower <- sum_cox$conf.int[3]
hr_upper <- sum_cox$conf.int[4]
p_val <- sum_cox$sctest[3]

# Create custom annotation text depicting p-value and Hazard Ratio (with 95% CI)
annotation_text <- paste0(
  "Log-rank p = ", round(p_val, 2), "\n",
  "HR = ", round(hr, 2), " (95% CI: ", round(hr_lower, 2), "-", round(hr_upper, 2), ")"
)

# Plot Kaplan-Meier curve
p <- ggsurvplot(
  fit,
  data = surv_df,
  pval = annotation_text,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  title = "CPTAC LUAD: Overall Survival by pNRF2 Status",
  legend.title = "pNRF2 Status",
  legend.labs = c("non pNRF2", "pNRF2"),
  palette = c("#ED8936", "#319795"),
  ggtheme = theme_pubr(base_size = 14)
)

# Save the plot
png("Plots/Survival/CPTAC_pNRF2_survival.png", width = 8, height = 6, units = "in", res = 600)
print(p)
dev.off()

message("Survival analysis completed successfully!")
