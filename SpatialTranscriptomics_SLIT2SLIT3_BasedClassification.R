library(Seurat)
library(ggplot2)
library(dplyr)


#Setting reproducible seed
set.seed(1234)


#Loading the fibroblast-only object
seurat_LS_fibro <- readRDS("LS_CAFs_PostHarmony.rds")


#Confirming that "Spatial" is the default assay
DefaultAssay(seurat_LS_fibro) <- "Spatial"


#Fetching SLIT2 and SLIT3 expression (log-normalized)
slit2_exp <- FetchData(seurat_LS_fibro, vars = "SLIT2")[, 1]
slit3_exp <- FetchData(seurat_LS_fibro, vars = "SLIT3")[, 1]


#Calculating non-zero means for SLIT2 and SLIT3, that is calculating means for cells that
#express the genes of interest at least once (so that there is at least one transcript), thus
#even though cells that don't express either of these two genes are kept for the analysis,
#they aren't taken for mean calculation 
slit2_nonzero_mean <- mean(slit2_exp[slit2_exp > 0])
slit3_nonzero_mean <- mean(slit3_exp[slit3_exp > 0])


cat("Mean (SLIT2 > 0):", slit2_nonzero_mean, "\n")
cat("Mean (SLIT3 > 0):", slit3_nonzero_mean, "\n")


#Classifying cells into five categories
SLIT_category <- case_when(
  slit2_exp > 0 & slit3_exp == 0 ~ "SLIT2+_SLIT3-",
  slit2_exp == 0 & slit3_exp > 0 ~ "SLIT2-_SLIT3+",
  slit2_exp > 0 & slit3_exp > 0 & slit2_exp > slit2_nonzero_mean & slit3_exp > slit3_nonzero_mean ~ "SLIT2+_SLIT3+_high",
  slit2_exp > 0 & slit3_exp > 0 & slit2_exp < slit2_nonzero_mean & slit3_exp < slit3_nonzero_mean ~ "SLIT2+_SLIT3+_low",
  slit2_exp == 0 & slit3_exp == 0 ~ "SLIT2-_SLIT3-",
  TRUE ~ "Unclassified"  
)


#Adding classification to Seurat metadata
seurat_LS_fibro$SLIT_Class <- SLIT_category


#Visualizing on UMAP
p <- DimPlot(seurat_LS_fibro, group.by = "SLIT_Class", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("SLIT2/SLIT3-Based Fibroblast Classification") +
  theme_minimal()
print(p)


#Showing class counts
table(seurat_LS_fibro$SLIT_Class)


#Saving updated Seurat object
saveRDS(seurat_LS_fibro, file = "LS_CAFs_SLITClassified.rds")