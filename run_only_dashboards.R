# Script to run only the clinical dashboard portion of Analysis.R
setwd("/mnt/Linux_storage/KNC")
message("Reading Analysis.R...")
analysis_lines <- readLines("Analysis.R")
message("Executing lines 1 to 532 of Analysis.R...")
# Since we modified lines, let us execute lines up to where create_clinical_dashboard(msk_clean...) is called
msk_call_idx <- grep("create_clinical_dashboard\\(msk_clean", analysis_lines)
if (length(msk_call_idx) > 0) {
  end_idx <- msk_call_idx[1]
  message(paste("Running up to line", end_idx))
  eval(parse(text = analysis_lines[1:end_idx]))
} else {
  message("Warning: MSK call not found, running up to 532")
  eval(parse(text = analysis_lines[1:532]))
}
message("Clinical dashboards generated successfully!")
