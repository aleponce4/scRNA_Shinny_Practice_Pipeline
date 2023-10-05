# Load Shiny and any other libraries
library(shiny)
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


# Load results
sced <- readRDS(file.path("../results", "sced_complete_analysis.rds"))
markers <- metadata(sced)$markers
logFCs_ordered <- metadata(sced)$logFCs_ordered

# source helper R scripts
source("additional_files/utility_functions.R")



###########################       Define server     ##########################
server <- function(input, output) {
    
    # Render the selected analysis plot
    output$analysis_plot <- renderPlot({
        
        # Select PCA
        if (input$analysis_type == "PCA") {
            # PCA plotting using sced object
            plot <- plotReducedDim(sced, dimred = "PCA", ncomponents = 4, colour_by = "nGene")
            print(plot)
        
        # Select Heatmap
        } else if (input$analysis_type == "Heatmap") {
            logFC_data <- metadata(sced)$logFCs_ordered[[input$cluster]]
            # Heatmap plotting code using logFC_data
            plot <- pheatmap(logFC_data)
            print(plot)
        
        # Select UMAP
        } else if (input$analysis_type == "UMAP") {
            plot <- plotReducedDim(sced, dimred = "UMAP", 
                                   colour_by = "snn_cluster") +
                theme(legend.text = element_text(size = 12), # adjust text size
                      legend.key.size = unit(1.5, "cm")) +    # adjust key size
                guides(colour = guide_legend(override.aes = list(size = 5)))  # adjust dot size in legend
            print(plot)
            
        # Select UMAP vs Marker
        } else if (input$analysis_type == "UMAP Grid") {
            plotlist <- list()
            
            for (i in 1:9) {
                marker <- input$selected_genes[i]
                plotlist[[paste0("Gene ", i, ": ", marker)]] <- 
                    plotReducedDim(sced, dimred = "UMAP", colour_by = marker, by_exprs_values = "logcounts") +
                    scale_fill_gradientn(colours = colorRampPalette(c("grey90", "orange3", "firebrick", "red"))(10)) +
                    ggtitle(label = paste0("Gene ", i, ": ", marker)) +
                    theme(plot.title = element_text(size = 12))
            }
            
            grid <- plot_grid(plotlist = plotlist, ncol = 3)
            print(grid)
        }
    })
    
    # Add dynamic context text
    output$context_text <- renderText({
        if (input$analysis_type == "Heatmap") {
            return(paste("Heatmap of the top marker genes for cluster", input$cluster, "in the dataset, stratified by cluster."))
        } else {
            return(NULL)
        }
    })
}


##########################    Define the UI     ##########################
ui <- fluidPage(
    # Main Title
    titlePanel("ScRNAseq Analysis Explorer"),
    
    # Sidebar layout for inputs and main panel for plots
    sidebarLayout(
        
        # Sidebar Inputs
        sidebarPanel(
            # Dropdown for selecting the analysis type
            selectInput("analysis_type", 
                        "Choose Analysis Type:", 
                        choices = c("PCA", "Heatmap", "UMAP", "UMAP with Gene Markers")),
            
            # Conditional dropdown for Heatmap - to select the cluster number
            conditionalPanel(
                condition = "input.analysis_type == 'Heatmap'",
                selectInput("cluster", 
                            "Choose a cluster:", 
                            choices = names(metadata(sced)$logFCs_ordered))
            ),
            

            # Dropdown for selecting up to 9 genes for "UMAP Grid"
            conditionalPanel(
                condition = "input.analysis_type == 'UMAP with Gene Markers'",
                selectInput("selected_genes", 
                            "Select a gene:", 
                            choices = rownames(markers),
                            selected = rownames(markers)[1:9],
                            multiple = TRUE)
            )
        ),
        
        # Main Panel for displaying plots and context text
        mainPanel(
            # Text output for context
            textOutput("context_text"),
            
            # Plot output
            plotOutput("analysis_plot", height = "800px", width = "1100px")
        )
    )
)


# Run app
shinyApp(ui = ui, server = server)

