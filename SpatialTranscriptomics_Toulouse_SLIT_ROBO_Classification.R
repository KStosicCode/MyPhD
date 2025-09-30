library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(harmony)

#Setting reproducible seed
set.seed(1234)


#Loading the newly subsetted object
seurat_LS_subset <- readRDS("LS_CAFs_Tumor_Immune_PostHarmony.rds")


#Confirming that "Spatial" is the default assay
DefaultAssay(seurat_LS_subset)
DefaultAssay(seurat_LS_subset) <- "Spatial"


#Classifying CAFs according to SLIT2 and SLIT3 expression, Tumor Epithelial cells according to
#ROBO1 expression and leaving Immune cells as they are. 
#CAFs
#Subsetting CAFs from the common Seurat object 
table(Idents(seurat_LS_subset), seurat_LS_subset$CellType_Specific)

table(Idents(seurat_LS_subset))

Idents(seurat_LS_subset) <- "CellType_Specific"

table(Idents(seurat_LS_subset))

cafs <- subset(seurat_LS_subset, idents = c("iCAF", "myCAF"))


#Re-pre-processing the Subsetted Seurat object 
cafs <- NormalizeData(cafs, verbose = FALSE)
cafs <- FindVariableFeatures(cafs, selection.method = "vst", nfeatures = 2000)
cafs <- ScaleData(cafs, verbose = FALSE)
cafs <- RunPCA(cafs, npcs = 30, verbose = FALSE)


#Re-run Harmony
cafs <- RunHarmony(
  object = cafs,
  group.by.vars = "orig.ident",
  plot_convergence = T,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  theta = 5,
  assay.use = "Spatial"
)

DimPlot(object = cafs, reduction = "harmony", group.by = "orig.ident")


#Final Non-linear dimension reduction and clustering
#Using the Harmony-corrected embeddings for downstream analysis
#Running UMAP and clustering on the Harmony components instead of the original PCs
cafs <- FindNeighbors(cafs, reduction = "harmony", dims = 1:30)
cafs <- FindClusters(cafs, resolution = 0.4)
cafs <- RunUMAP(cafs, reduction = "harmony", dims = 1:30)


#UMAP plots to confirm integration and annotation
DimPlot(cafs, reduction = "umap", group.by = "orig.ident", label = TRUE) + 
  ggtitle("LS_CAFs UMAP by Sample")

DimPlot(cafs, reduction = "umap", group.by = "CellType_Specific", label = TRUE) +
  ggtitle("LS_CAFs UMAP")

#Classifying CAFs according to SLIT2 and SLIT3 genes expression
table(cafs$CellType_Specific)

#Fetching SLIT2 and SLIT3 expression (log-normalized)
slit2_exp <- FetchData(cafs, vars = "SLIT2")[, 1]
slit3_exp <- FetchData(cafs, vars = "SLIT3")[, 1]


#Calculating non-zero means for SLIT2 and SLIT3, that is calculating means for cells that
#express the genes of interest at least once (so that there is at least one transcript), thus
#even though cells that don't express either of these two genes are kept for the analysis,
#they aren't taken for mean calculation 
slit2_nonzero_mean <- mean(slit2_exp[slit2_exp > 0])
slit3_nonzero_mean <- mean(slit3_exp[slit3_exp > 0])

#Classifying CAF cells into five categories according to SLIT2 and SLIT3 genes expression
SLIT_category <- case_when(
  slit2_exp > 0 & slit3_exp == 0 ~ "SLIT2+_SLIT3-",
  slit2_exp == 0 & slit3_exp > 0 ~ "SLIT2-_SLIT3+",
  slit2_exp > 0 & slit3_exp > 0 & slit2_exp > slit2_nonzero_mean & slit3_exp > slit3_nonzero_mean ~ "SLIT2+_SLIT3+_high",
  slit2_exp > 0 & slit3_exp > 0 & slit2_exp < slit2_nonzero_mean & slit3_exp < slit3_nonzero_mean ~ "SLIT2+_SLIT3+_low",
  slit2_exp == 0 & slit3_exp == 0 ~ "SLIT2-_SLIT3-"  
)

#Adding classification to Seurat metadata
cafs$SLIT_Class <- SLIT_category


#Visualizing on UMAP
p <- DimPlot(cafs, group.by = "SLIT_Class", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("SLIT2/SLIT3-Based Fibroblast Classification") +
  theme_minimal()

print(p)

#Showing class counts
table(cafs$SLIT_Class)

#Applying SLIT classification from cafs on the main dataset
#Initializing the Class column in the main Seurat object and setting default label for all cells,
#which can be overwritten by values like iCAF, myCAF and Tumor Epithelial cells
seurat_LS_subset@meta.data$Class <- "Immune_cells"
head(seurat_LS_subset@meta.data$Class)
tail(seurat_LS_subset@meta.data$Class)

#Assigning SLIT-based CAF (fibroblast) classifications by matching cell names
#Getting a character vector of CAF cell names (barcodes) from the subset .
fibro_cells <- rownames(cafs@meta.data) 
length(fibro_cells)

#Using a vector with individual caf names to index seurat_LS_subset@meta.data, thus ensuring
#that each matching cell gets the correct label.
#Mapping the SLIT-based classification from the subset back to the main object
seurat_LS_subset@meta.data[fibro_cells, "Class"] <- cafs$SLIT_Class

x <- seurat_LS_subset


#Tumor Epithelial cells
#Subsetting Tumor Epithelial cells from the common Seurat object 
table(Idents(seurat_LS_subset), seurat_LS_subset$CellType_Specific)

table(Idents(seurat_LS_subset))

Idents(seurat_LS_subset) <- "CellType_Specific"

table(Idents(seurat_LS_subset))

tecs <- subset(seurat_LS_subset, idents = "Tumor Epithelial cells")


#Re-pre-processing the Subsetted Seurat object 
tecs <- NormalizeData(tecs, verbose = FALSE)
tecs <- FindVariableFeatures(tecs, selection.method = "vst", nfeatures = 2000)
tecs <- ScaleData(tecs, verbose = FALSE)
tecs <- RunPCA(tecs, npcs = 30, verbose = FALSE)


#Re-run Harmony
tecs <- RunHarmony(
  object = tecs,
  group.by.vars = "orig.ident",
  plot_convergence = T,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  theta = 5,
  assay.use = "Spatial"
)

DimPlot(object = tecs, reduction = "harmony", group.by = "orig.ident")


#Final Non-linear dimension reduction and clustering
#Using the Harmony-corrected embeddings for downstream analysis
#Running UMAP and clustering on the Harmony components instead of the original PCs
tecs <- FindNeighbors(tecs, reduction = "harmony", dims = 1:30)
tecs <- FindClusters(tecs, resolution = 0.4)
tecs <- RunUMAP(tecs, reduction = "harmony", dims = 1:30)


#UMAP plots to confirm integration and annotation
DimPlot(tecs, reduction = "umap", group.by = "orig.ident", label = TRUE) + 
  ggtitle("LS_Tumor Epithelial cells UMAP by Sample")

DimPlot(tecs, reduction = "umap", group.by = "CellType_Specific", label = TRUE) +
  ggtitle("LS_Tumor Epithelial cells UMAP")

#Classifying Tumor Epithelial cells according to ROBO1 genes expression
table(tecs$CellType_Specific)

#Fetching ROBO1 expression (log-normalized)
robo1_exp <- FetchData(tecs, vars = "ROBO1")[, 1]

#Classifying Tumor Epithelial cells into two categories according to ROBO1 gene expression
ROBO_category <- case_when(
  robo1_exp > 0 ~ "ROBO1+",
  robo1_exp == 0 ~ "ROBO1-"
)

#Adding classification to Seurat metadata
tecs$ROBO_Class <- ROBO_category


#Visualizing on UMAP
p <- DimPlot(tecs, group.by = "ROBO_Class", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("ROBO-Based Tumor Epithelial cells Classification") +
  theme_minimal()

print(p)

#Showing class counts
table(tecs$ROBO_Class)

#Applying ROBO classification from tecs on the main dataset
head(seurat_LS_subset@meta.data$Class)
tail(seurat_LS_subset@meta.data$Class)
table(seurat_LS_subset$Class)

#Assigning ROBO-based Tumor Epithelial cells classifications by matching cell names
#Getting a character vector of Tumor Epithelial cells names (barcodes) from the subset.
tumor_cells <- rownames(tecs@meta.data)
length(tumor_cells)

#Using a vector with individual tecs names to index seurat_LS_subset@meta.data, thus ensuring
#that each matching cell gets the correct label.
#Mapping the ROBO-based classification from the subset back to the main object
seurat_LS_subset@meta.data[tumor_cells, "Class"] <- tecs$ROBO_Class

x <- seurat_LS_subset
y <- x@meta.data


#Checking the new classification
table(seurat_LS_subset$Class)
table(seurat_LS_subset$Class, seurat_LS_subset$CellType_Specific)


#Saving the updated Seurat object
saveRDS(seurat_LS_subset, file = "LS_Classified_CAFs_Tumor_Immune.rds")