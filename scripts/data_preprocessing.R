# List of packages to check/install
packages_needed <- c("scuttle", "scran", "SingleCellExperiment", "hexbin", "gridExtra", "BiocParallel", "ggplot2")

# Check which ones are not installed
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]

# If any are not installed, install them
if(length(packages_to_install) > 0) {
  BiocManager::install(packages_to_install)
}

# Load necessary libraries
library(scuttle)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(BiocParallel)


###### Setp 0: Load the raw data

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
# Set the working directory
setwd(get_script_dir())

# Construct the data path relative to the script's directory
data_path <- file.path("..", "data", "raw", "kotliarovPBMCData.rds")

# Read the RDS file
pbmc_data <- readRDS(file = data_path)


########### Step 1: Quality Control

plot_QC <- function(sce_object) {
  # Extract cell metadata for QC
  cell_metadata <- colData(sce_object)
  cell_metadata <- as.data.frame(cell_metadata)
  
  # Visualize nGene, nUMI, and pctMT metrics
  # color palette
  palette <- scales::brewer_pal(palette = "Set2")(3)
  
  # Create the plots
  plot_nGene <- ggplot(cell_metadata, aes(x = factor(0), y = nGene, fill = "nGene")) + 
    geom_violin() + 
    labs(title = "Detected Genes per Cell", x = "") +
    theme_light() +
    scale_fill_manual(values = palette[1], guide = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  plot_nUMI <- ggplot(cell_metadata, aes(x = factor(0), y = nUMI, fill = "nUMI")) + 
    geom_violin() + 
    labs(title = "Transcripts per Cell", x = "") +
    theme_light() +
    scale_fill_manual(values = palette[2], guide = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  plot_pctMT <- ggplot(cell_metadata, aes(x = factor(0), y = pctMT, fill = "pctMT")) + 
    geom_violin() + 
    labs(title = "Mitochondrial Genes Percentage", x = "") +
    theme_light() +
    scale_fill_manual(values = palette[3], guide = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  plot_scatter <- ggplot(data = cell_metadata, aes(x = nUMI, y = nGene, color = pctMT)) + 
    geom_point(alpha = 0.6) + 
    stat_smooth(method = lm, se = FALSE, color = "black") +
    scale_x_log10() + 
    scale_y_log10() + 
    labs(color = "Mitochondrial\nGene Ratio") +
    theme_light()
  
  
  # Arrange the plots in a 2x2 grid
  grid.arrange(plot_nGene, plot_nUMI, plot_pctMT, plot_scatter, ncol=2)
}

# Plot before filtering
plot_QC(pbmc_data)

# Extract cell metadata for QC
cell_metadata <- colData(pbmc_data)
cell_metadata <- as.data.frame(cell_metadata)

# nGene thresholds:
# - Lower threshold: Removes cells with an exceptionally low number of genes detected,
#   which may be indicative of dying cells or cells that didn't get captured properly.
# - Upper threshold: Excludes cells with an abnormally high number of genes detected,
#   which could be indicative of doublets (two cells mistakenly identified as one) or 
#   other technical artifacts.
upper_nGene_threshold <- 3000 
lower_nGene_threshold <- quantile(cell_metadata$nGene, 0.1)



# nUMI thresholds:
# - Lower threshold: Removes cells with very few UMIs, which might suggest poor-quality cells or 
#   issues in cell capture or library preparation.
# - Upper threshold: Filters out cells with an extremely high count of UMIs, which might suggest
#   that multiple cells have been captured together or other potential technical issues.
lower_nUMI_threshold <- quantile(cell_metadata$nUMI, 0.1)
upper_nUMI_threshold <- 10000

max_pctMT <- 0.2  # retain cells with less than 20% mitochondrial genes

# Filter based on quality metrics
good_cells <- which(
  cell_metadata$nGene > lower_nGene_threshold & 
    cell_metadata$nGene < upper_nGene_threshold & 
    cell_metadata$nUMI > lower_nUMI_threshold & 
    cell_metadata$nUMI < upper_nUMI_threshold & 
    cell_metadata$pctMT < max_pctMT
)
pbmc_data <- pbmc_data[,good_cells]

# Plot after filtering
plot_QC(pbmc_data)


########### Step 2: Normalization

# Normalize using scran pooling-based normalization
clusters <- quickCluster(pbmc_data)
pbmc_data <- computeSumFactors(pbmc_data, cluster=clusters)
pbmc_data <- logNormCounts(pbmc_data)
summary(sizeFactors(pbmc_data))


########### Step 3: Highly Variable genes

# add logcounts
pbmc_data <- logNormCounts(pbmc_data)

# model gene variance
dec <- modelGeneVar(pbmc_data)

# Plot variance
plot(dec$mean, dec$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

# Select highly variable genes

# Define marker genes   # You can modify this part to manually include any genes you ar eitnerested that might be lost during filtering
marker_genes <- c("ENSG00000198851", "ENSG00000010610", "ENSG00000153563", 
                  "ENSG00000115523", "ENSG00000105374", "ENSG00000156738", 
                  "ENSG00000170458", "ENSG00000090382", "ENSG00000101439", 
                  "ENSG00000166927", "ENSG00000203747", "ENSG00000090382", 
                  "ENSG00000101439", "ENSG00000179639", "ENSG00000101439")

# Get genes with FDR below 5%
genes_fdr <- getTopHVGs(dec, fdr.threshold=0.05)

# Exclude the genes with FDR below 5% from the original set and then get the top genes
remaining_genes <- setdiff(rownames(dec), union(genes_fdr, marker_genes))
top_remaining_genes <- head(remaining_genes, 2000 - length(genes_fdr) - length(marker_genes))

# Combine the genes: FDR genes + top genes + marker genes
final_top_hvgs <- union(union(genes_fdr, top_remaining_genes), marker_genes)

# Filtering the pbmc_data to retain only the top HVGs
# Dimensions before filtering
print(dim(pbmc_data))

# Subset the data
pbmc_data_hvg <- pbmc_data[final_top_hvgs, ]

# Dimensions after filtering
print(dim(pbmc_data_hvg))

# Save data for further analysis
save_path <- file.path("..", "data", "processed", "kotliarovPBMCData.rds")
# Save the RDS file
saveRDS(pbmc_data_hvg, file = save_path)







