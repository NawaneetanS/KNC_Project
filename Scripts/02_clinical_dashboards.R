#!/usr/bin/Rscript

# 02_clinical_dashboards.R
# Loads cleaned clinical data and generates dashboards for each cohort.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(janitor)
})

setwd("/media/nannu1375/Backpack/Shankara/KNC")

# Define a generic clinical dashboard function
create_clinical_dashboard <- function(df, cohort_name, plot_vars, outfile) {
  plot_list <- list()
  
  for (var_info in plot_vars) {
    col_name <- var_info$col
    var_label <- var_info$label
    var_type <- var_info$type
    
    # Filter out NAs for this plot
    plot_df <- df %>% 
      filter(!is.na(.data[[col_name]]), as.character(.data[[col_name]]) != "", as.character(.data[[col_name]]) != "NA")
    
    # Skip if no data
    if (nrow(plot_df) == 0) next
    
    if (var_type == "continuous") {
      # Make sure it is numeric
      plot_df[[col_name]] <- as.numeric(plot_df[[col_name]])
      plot_df <- plot_df %>% filter(!is.na(.data[[col_name]]))
      
      fill_col <- if (!is.null(var_info$fill)) var_info$fill else "#2C7FB8"
      
      if (col_name == "AGE") {
        # Keep density curve for age
        p <- ggplot(plot_df, aes(x = .data[[col_name]])) +
          geom_histogram(aes(y = after_stat(density)), fill = fill_col, color = "white", alpha = 0.6, bins = 15) +
          geom_density(color = "#1D3557", linewidth = 1.2, fill = fill_col, alpha = 0.1) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = var_label, y = "Density")
      } else {
        # Standard frequency histogram without density line
        p <- ggplot(plot_df, aes(x = .data[[col_name]])) +
          geom_histogram(fill = fill_col, color = "white", alpha = 0.85, bins = 15) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = var_label, y = "Patient Count")
      }
      
      p <- p + theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        axis.line = element_line(color = "grey80"),
        plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
      )
      
    } else if (var_type == "categorical") {
      # If not already a factor, convert to factor and order by frequency
      if (!is.factor(plot_df[[col_name]])) {
        plot_df[[col_name]] <- factor(plot_df[[col_name]])
        plot_df[[col_name]] <- reorder(plot_df[[col_name]], plot_df[[col_name]], FUN = length)
      }
      
      cat_counts <- plot_df %>% count(.data[[col_name]], name = "n")
      num_cats <- nrow(cat_counts)
      
      if (num_cats > 6) {
        # Horizontal bar plot for many categories to prevent label overlapping
        p <- ggplot(plot_df, aes(y = .data[[col_name]], fill = .data[[col_name]])) +
          geom_bar(width = 0.7, alpha = 0.85, color = "white", linewidth = 0.3) +
          geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.2, fontface = "bold", size = 3) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
          scale_y_discrete(labels = scales::label_wrap(15)) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = "Patient Count", y = NULL) +
          theme(
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey95"),
            axis.line = element_line(color = "grey80"),
            plot.margin = margin(t = 15, r = 15, b = 15, l = 10)
          )
      } else {
        # Vertical bar plot for few categories
        p <- ggplot(plot_df, aes(x = .data[[col_name]], fill = .data[[col_name]])) +
          geom_bar(width = 0.6, alpha = 0.85, color = "white", linewidth = 0.3) +
          geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, fontface = "bold", size = 3) +
          scale_x_discrete(labels = scales::label_wrap(12)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
          theme_minimal(base_size = 11) +
          labs(title = var_label, x = NULL, y = "Patient Count") +
          theme(
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b=8)),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey95"),
            axis.line = element_line(color = "grey80"),
            plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
          )
      }
      
      # Apply custom palette if provided
      if (!is.null(var_info$palette)) {
        p <- p + scale_fill_manual(values = var_info$palette, breaks = names(var_info$palette))
      } else {
        if (num_cats <= 4) {
          p <- p + scale_fill_manual(values = c("#7FCDBB", "#41B6C4", "#1D91C0", "#081D58", "#E2E8F0"))
        } else {
          p <- p + scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8)
        }
      }
    }
    
    plot_list[[col_name]] <- p
  }
  
  num_plots <- length(plot_list)
  num_cols <- min(3, num_plots)
  num_rows <- ceiling(num_plots / num_cols)
  
  # Combine using patchwork
  dashboard <- wrap_plots(plot_list, ncol = num_cols) +
    plot_annotation(
      title = paste(cohort_name, "LUAD Clinical & Pathological Dashboard"),
      subtitle = paste0("Total Cohort Size: N = ", nrow(df), " Patients"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 12, b = 4)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#4A5568", margin = margin(b = 12))
      )
    )
  
  png(filename = outfile, width = 4.5 * num_cols, height = 3.5 * (num_rows + 0.2), units = "in", res = 600)
  print(dashboard)
  dev.off()
}

message("--- Generating Clinical Dashboards ---")

# Load intermediate clinical data
cptac_clean  <- readRDS("Tables/cptac_clean.rds")

# 1. CPTAC
cptac_vars <- list(
  list(col = "AGE", type = "continuous", label = "Age (Years)"),
  list(col = "SEX", type = "categorical", label = "Gender", palette = c("Male" = "#2C7FB8", "Female" = "#D53F8C")),
  list(col = "STAGE", type = "categorical", label = "Pathological Stage", palette = c("I" = "#7FCDBB", "II" = "#41B6C4", "III" = "#1D91C0", "IV" = "#081D58", "Unknown" = "#E2E8F0")),
  list(col = "SMOKING_STATUS", type = "categorical", label = "Smoking Status", palette = c("No" = "#319795", "Yes" = "#ED8936", "Unknown" = "#E2E8F0")),
  list(col = "Subtype", type = "categorical", label = "Adenocarcinoma Subtype"),
  list(col = "Purity", type = "continuous", label = "Tumor Purity", fill = "#319795"),
  list(col = "TMB_NONSYNONYMOUS", type = "continuous", label = "TMB (Nonsynonymous)", fill = "#D95F02")
)
create_clinical_dashboard(cptac_clean, "CPTAC", cptac_vars, "Plots/Clinical/CPTAC_clinical_dashboard.png")

message("Clinical dashboards generated successfully!")
