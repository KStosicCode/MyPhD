library(Seurat)          
library(harmony)          
library(dplyr)            
library(ggplot2)         
library(patchwork)        


#Setting a reproducible seed
set.seed(1234)

#Loading the LS Seurat Object 
seurat_LS <- readRDS("LongOS_RCTD_Harmony_ISCHIA__PDAC_10Samples_SlitRobo_SeuratAnn.rds")


#Inspecting key metadata
cat("Samples present in orig.ident:\n")
print(unique(seurat_LS$orig.ident))

cat("\nCellType_Specific annotations present:\n")
print(unique(seurat_LS$CellType_Specific))

DefaultAssay(seurat_LS)
seurat_LS@reductions$pca@assay.used


#Re-run PCA before Harmony 
seurat_LS <- NormalizeData(seurat_LS, verbose = FALSE)
seurat_LS <- FindVariableFeatures(seurat_LS, selection.method = "vst", nfeatures = 2000)
seurat_LS <- ScaleData(seurat_LS, verbose = FALSE)
seurat_LS <- RunPCA(seurat_LS, npcs = 30, verbose = FALSE)


#Harmony Integration 
#Integrate the 4 LS samples using Harmony with sample-level batch correction
#Harmony aligns the four LS samples in a shared space, 
#mitigating batch or sample-specific effects
#It directly adjusts the PCA embeddings to produce a harmony embedding for each cell
seurat_LS <- RunHarmony(
  object = seurat_LS,
  group.by.vars = "orig.ident",
  plot_convergence = T,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  theta = 5,
  assay.use = "Spatial"
)

DimPlot(seurat_LS, reduction = "harmony", group.by = "orig.ident")


x <- LongOS_RCTD_Harmony_ISCHIA__PDAC_10Samples_SlitRobo_SeuratAnn

#x <- NormalizeData(x) %>% FindVariableFeatures() %>% 
  #ScaleData() %>% RunPCA(verbose = FALSE)

#x <- RunHarmony(
  #object = x,
  #group.by.vars = "orig.ident",
  #plot_convergence = T,
  #nclust = 50,
  #max_iter = 10,
  #early_stop = TRUE,
  #theta = 3,
  #assay.use = "Spatial"
#)

#x <- RunUMAP(x, reduction = "harmony", dims = 1:30)
#x <- FindNeighbors(x, reduction = "harmony", dims = 1:30)
#x <- FindClusters(x, resolution = 0.5)

#DimPlot(x, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("UMAP by Sample")
#DimPlot(x, reduction = "umap", group.by = "CellType_Specific", label = TRUE) + ggtitle("UMAP by Cell Type")


#Non-Linear Dimensionality Reduction & Clustering 
#Using the Harmony-corrected embeddings for downstream analysis
#Running UMAP and clustering on the Harmony components instead of the original PCs
seurat_LS <- FindNeighbors(seurat_LS, reduction = "harmony", dims = 1:30)
seurat_LS <- FindClusters(seurat_LS, resolution = 0.5)
seurat_LS <- RunUMAP(seurat_LS, reduction = "harmony", dims = 1:30)

#Plot to visually inspect clustering and integration
DimPlot(seurat_LS, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("UMAP by Sample")
DimPlot(seurat_LS, reduction = "umap", group.by = "CellType_Specific", label = TRUE) + ggtitle("UMAP by Cell Type")


#Save LS-only object for downstream analyses
saveRDS(seurat_LS, file = "LS_PostHarmony.rds")


#Subset Fibroblast Populations 
#Using CellType_Specific: Keep only iCAFs and myCAFs
#Original approach
Idents(seurat_LS)
table(Idents(seurat_LS))

#Cross-tabulation to show how each cluster overlaps with each cell type, to prove that
#clustering and annotation are overlapping but not identical.
table(Idents(seurat_LS), seurat_LS$CellType_Specific)


Idents(seurat_LS) <- "CellType_Specific"
seurat_LS_fibro <- subset(seurat_LS, idents = c("iCAF", "myCAF"))

#Alternative approach generating the same subsetted Seurat Object as the original one
#fibroblast_types <- c("iCAF", "myCAF")
#seurat_LS_fibro <- subset(seurat_LS, subset = CellType_Specific %in% fibroblast_types)


#Re-integrate the Fibroblast Subset 
seurat_LS_fibro <- NormalizeData(seurat_LS_fibro, verbose = FALSE)
seurat_LS_fibro <- FindVariableFeatures(seurat_LS_fibro, selection.method = "vst", nfeatures = 2000)
seurat_LS_fibro <- ScaleData(seurat_LS_fibro, verbose = FALSE)
seurat_LS_fibro <- RunPCA(seurat_LS_fibro, npcs = 30, verbose = FALSE)


#Re-run Harmony
seurat_LS_fibro <- RunHarmony(
  object = seurat_LS_fibro,
  group.by.vars = "orig.ident",
  plot_convergence = T,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  theta = 5,
  assay.use = "Spatial"
)

DimPlot(seurat_LS_fibro, reduction = "harmony", group.by = "orig.ident")


#Final Non-linear dimension reduction and clustering
#Using the Harmony-corrected embeddings for downstream analysis
#Running UMAP and clustering on the Harmony components instead of the original PCs
seurat_LS_fibro <- FindNeighbors(seurat_LS_fibro, reduction = "harmony", dims = 1:30)
seurat_LS_fibro <- FindClusters(seurat_LS_fibro, resolution = 0.4)
seurat_LS_fibro <- RunUMAP(seurat_LS_fibro, reduction = "harmony", dims = 1:30)


#UMAP plots to confirm integration and annotation
DimPlot(seurat_LS_fibro, reduction = "umap", group.by = "orig.ident", label = TRUE) + 
  ggtitle("Fibroblast UMAP by Sample")
DimPlot(seurat_LS_fibro, reduction = "umap", group.by = "CellType_Specific", label = TRUE) +
  ggtitle("Fibroblast UMAP (iCAF and myCAF only)")


#Save fibroblast-only object for downstream analyses
saveRDS(seurat_LS_fibro, file = "LS_CAFs_PostHarmony.rds")