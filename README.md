# SingleCell-RNAseq--Shinny -TrainingToolkit

This pipeline is designed for the analysis and visualization of single-cell RNA sequencing (scRNA-seq) data. The pipeline was tested using publicly available scRNA-seq datasets from Bioconductor, but it should be adaptable to any scRNA-seq experiment when provided with an appropriate object.

## Directory Structure
```
scRNA_Shinny_Practice_Pipeline/
│
├── data/ # All data files
│ ├── raw/ # Raw data files
│ └── processed/ # Processed data after scripts
│
├── results/ # Results files and plots from analysis
│
├── scripts/ # Scripts for data acquisition, analysis, & processing
│
└── shiny_app/ # Shiny application for visualization
```

## Usage

1. **Data Acquisition**: 
   - Run the `data_acquisition.R` script. 
   - If you want to use a different data file, make sure to modify the file path in this script.


**Data Analysis**: 
   - Run the `data_analysis.R` script.
   - **Variance Modeling**: Carefully inspect the variance plots. Adjust settings or filters as needed based on the plot's insights.
    ![Variance Plot](/Plots/2.png)
   - **Data Quality Assessment**: Check the grid of distributions in to determine filtering thresholds.
    ![Variance Plot](/Plots/1.png)
   - This grid provides insights into:
     - Detected genes per cell
     - Mitochondrial genes per cell
     - Transcripts per cell
     - nGene vs. nUMI scatter plot colored by the mitochondrial ratio

3. **Data Processing**: 
   - Run the `data_processing.R` script. 
   - This script generates a file that will be used by the Shiny application for visualization.

4. **Shiny App Visualization**: 
   - After processing the data, launch the Shiny application to explore various visualizations.
   - The app allows you to view heatmaps, UMAP scatterplots, and other relevant outputs.
  ![Umap colored by gene markers](/Plots/4.png)
   - Besides the standard UMAP visualization, you have the option to select primary gene markers to color the cells, providing deeper insights into cell clustering based on gene expression.
  ![Umap colored by gene markers](/Plots/6.png)

## Dependencies

Before running the scripts, ensure you have the required packages installed. The project relies on packages from both CRAN and Bioconductor.

#### CRAN:
Install the necessary CRAN packages using the following command:
```
install.packages(c("shiny", "ggplot2", "gridExtra", "cluster", "scater", "uwot", "cowplot", "igraph", "RColorBrewer", "pheatmap"))
```

#### Bioconductor
First, install BiocManager if you haven't:
```
install.packages("BiocManager")
```

Then, install the required Bioconductor packages:
```
BiocManager::install(c("scRNAseq", "scuttle", "scran", "SingleCellExperiment", "org.Hs.eg.db", "BiocParallel"))
```

