# Define a function to get the script directory
get_script_dir <- function() {
  # First, check if RStudio's API can be used
  if ("rstudioapi" %in% installed.packages() && rstudioapi::isAvailable("0.99.467")) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  # Else, try using commandArgs()
  script_path <- commandArgs(trailingOnly = FALSE)
  script_path <- script_path[grepl("^--file=", script_path)]
  if (length(script_path) > 0) {
    return(dirname(sub("^--file=", "", script_path)))
  }
  stop("Cannot determine script's directory in this environment.")
}

# Set the working directory to the script's directory
setwd(get_script_dir())

# Check for BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Install scRNAseq using BiocManager
if (!requireNamespace("scRNAseq", quietly = TRUE)) {
  BiocManager::install("scRNAseq")
}

# Check if the intended directory exists, if not, create it
if (!dir.exists("../data/raw")) {
  dir.create("../data/raw", recursive = TRUE)
}

# Install the package containing the datasets
library(scRNAseq)

dataset_functions <- list(
  "KotliarovPBMCData" = KotliarovPBMCData,
  "ZhongPrefrontalData" = ZhongPrefrontalData,
  "XinPancreasData" = XinPancreasData,
  "LaMannoBrainData_human_ips" = function() LaMannoBrainData('human-ips'),
  "LaMannoBrainData_human_es" = function() LaMannoBrainData('human-es'),
  "BacherTCellData" = BacherTCellData,
  "DarmanisBrainData" = DarmanisBrainData,
  "LedergorMyelomaData" = LedergorMyelomaData,
  "ReprocessedFluidigmData" = ReprocessedFluidigmData
)

for (name in names(dataset_functions)) {
  cat("Processing", name, "...\n")

  # Get the dataset using the function
  sce_data <- dataset_functions[[name]]()
  
  # Save to a raw data folder as an RDS file
  save_path <- paste0("../data/raw/", name, ".rds")
  saveRDS(sce_data, file = save_path)
  
  cat(name, "saved to", save_path, "\n")
}
