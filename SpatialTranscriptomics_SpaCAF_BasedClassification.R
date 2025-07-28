library(Seurat)
library(dplyr)
library(openxlsx)    
library(HGNChelper)  


#Loading the Seurat object (fibroblast subset after Harmony integration)
seurat_LS_fibro <- readRDS("LS_CAFs_PostHarmony.rds")


# Source scType functions (gene set preparation and scoring utilities)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")


#Running scType with custom spaCAF marker signatures
seurat_LS_fibro <- run_sctype(
  seurat_LS_fibro,
  assay = "Spatial",              
  scaled = TRUE,               
  known_tissue_type = "Pancreas",        
  custom_marker_file = "C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Toulouse/spaCAF.xlsx",    
  name = "spaCAF_class",       
  plot = FALSE                 
)


#Viewing results, cell-type by cell matrix. See the complete example below
View(seurat_LS_fibro)


#Checking the distinct spaCAF labels assigned
table(seurat_LS_fibro$spaCAF_class)


DimPlot(seurat_LS_fibro, reduction = "umap", group.by = "spaCAF_class")


#Saving updated Seurat object
saveRDS(seurat_LS_fibro, file = "LS_CAFs_SpaCAFClassified.rds")