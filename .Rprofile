# Project-level R configuration for KNC
# Reads .renv and strips user home library (~/R/...) from .libPaths()

if (file.exists(".renv")) {
  readRenviron(".renv")
}

# Remove user home library directory from .libPaths()
.libPaths(.libPaths()[!grepl("^/home/[^/]+/R/", .libPaths())])
