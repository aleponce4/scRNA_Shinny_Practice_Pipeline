# List of packages to check/install
packages_needed <- c("uwot", "scatter", "igraph","org.Hs.eg.db")

# Check which ones are not installed
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]

# If any are not installed, install them
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if(length(packages_to_install) > 0) {
  BiocManager::install(packages_to_install)
}

# Load necessary libraries
library(scuttle)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(cluster)
library(scater)
library(uwot)
library(cowplot)
library(igraph)
library(org.Hs.eg.db)
library(RColorBrewer)
library(pheatmap)

###### Setp 0: Load the processed data

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
processed_files <- list.files(path = "../data/processed", pattern = "*.rds", full.names = TRUE)
raw_files <- list.files(path = "../data/raw", pattern = "*.rds", full.names = TRUE)

# Read the RDS file
pbmc_data <- readRDS(file = processed_files)

########## Step 1: Automatic PCA dimension reduction

# Model gene variance
dec <- modelGeneVar(pbmc_data)

# Denoise using PCA
sced <- denoisePCA(pbmc_data, dec, subset.row = NULL) # We already subseted for HVGs

# Determine the number of principal components to retain
output <- getClusteredPCs(reducedDim(sced))
npcs <- metadata(output)$chosen

# Subset the PCA result to the chosen number of PCs
reducedDim(sced, "PCAsubset") <- reducedDim(sced, "PCA")[, 1:npcs, drop = FALSE]

# Get PCA results
pca_results <- reducedDim(sced, "PCAsubset")

# Downsample the Data
set.seed(123)  # for reproducibility
n_cells_to_sample <- 500  # number of cells to sample
sced_sampled <- sced[, sample(1:ncol(sced), n_cells_to_sample)]


########## Step 2: Clustering

# Calculate the SNN graph
snn_graph <- buildSNNGraph(sced, use.dimred = "PCA")

# Detecting communities using Louvain algorithm
clustering_result <- cluster_louvain(snn_graph)

# Add the clustering result to the colData for visualization
colData(sced)$snn_cluster <- factor(clustering_result$membership)


########## Setp 3: UMAP dimensionality reduction

# Compute UMAP using the PCA results
umap_results <- umap(reducedDim(sced, "PCA"))

# Add UMAP results to colData for easy access and visualization
colData(sced)$UMAP_1 <- umap_results[, 1]
colData(sced)$UMAP_2 <- umap_results[, 2]

# Create a matrix from the UMAP data
umap_matrix <- as.matrix(colData(sced)[, c("UMAP_1", "UMAP_2")])

# Assign the UMAP matrix to the 'UMAP' slot in the reducedDims of sced object
reducedDims(sced)$UMAP <- umap_matrix


########## Step 4: Find markers

# Change ensembl ID for gene name
genes <- rownames(sced)
gene_symbols <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace NAs with original Ensemble IDs
gene_symbols[is.na(gene_symbols)] <- genes[is.na(gene_symbols)]
rownames(sced) <- gene_symbols

# Assign clusters to col labels
colLabels(sced) <- colData(sced)$snn_cluster

# Find markers function
markers <- findMarkers(sced, direction="up")

# Store markers in metadata slot of sced object
metadata(sced)$markers <- markers

# Extract and save the log-fold changes for the top 20 marker genes for each cluster and store in the sced object
list_logFCs_ordered <- list()

for(cluster in names(markers)) {
  marker.set <- markers[[cluster]]
  
  # Extract the log-fold changes for the top 30 marker genes
  logFCs <- getMarkerEffects(marker.set[1:30,])
  
  # Order clusters for heatmap
  ordered_clusters <- order(as.numeric(colnames(logFCs)))
  logFCs_ordered <- logFCs[, ordered_clusters]
  
  list_logFCs_ordered[[cluster]] <- logFCs_ordered
}

# Store ordered logFCs in metadata slot of sced object
metadata(sced)$logFCs_ordered <- list_logFCs_ordered

########## Step 5: Save results for Shinny App analysis

saveRDS(sced, file = "../results/sced_complete_analysis.rds")




