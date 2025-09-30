library(Seurat)
library(tidyverse)
library(pheatmap)
library(readxl)
library(reshape2)


#Loading the Seurat object
seurat_LS_fibro <- readRDS("LS_CAFs_PostHarmony.rds")
DefaultAssay(seurat_LS_fibro) <- "Spatial"


#SLIT genes of interest
slit_genes <- c("SLIT2", "SLIT3")


#Preparing to store signature scores
signature_scores <- list()


#Loading a path to the Excel file 
signatures_table_path <- "C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Toulouse/gene_sets_Zeynep2.xlsx"


#genes_to_check <- c("NCRNA00008", "LINC00008", "D11S813E", "MIR675HG")
#genes_to_check %in% rownames(seurat_LS_fibro)
#"DHHC-11" %in% rownames(seurat_LS_fibro) 


#Looping over each sheet (excel_sheets() function loads in all the sheets, whereas read_excel()
#function loads in only the first sheet and it loads it in as a data frame)
for (sheet in readxl::excel_sheets(signatures_table_path)) {
  
  #Loading the sheet as a data frame
  marker_table <- readxl::read_excel(signatures_table_path, sheet = sheet)
  
  #Looping over each column (CAF signature) in that sheet
  for (sig_col in colnames(marker_table)) {
    
    genes <- marker_table[[sig_col]]
    genes <- toupper(genes)
    genes <- genes[!is.na(genes)]
    
    if (length(genes) > 0) {
      
      #Creating unique name: "sheet_column"
      sig_name <- paste(sheet, sig_col, sep = "_")
      
      #Adding module score
      seurat_LS_fibro <- AddModuleScore(
        object = seurat_LS_fibro,
        features = list(genes),
        name = sig_name
      )
      
      #Storing generated column name (e.g. "sheet_col1")
      signature_scores[[sig_name]] <- paste0(sig_name, "1")
    }
  }
}


#Extracting expression data containing SLIT genes + Module scores for each CAF signature
expr_data <- FetchData(seurat_LS_fibro, vars = c(slit_genes, unlist(signature_scores)))


#Computing Pearson correlations
#Function expand.grid() makes all combinations of CAF signature × SLIT gene
#Adding empty columns cor and pval to later store Pearson correlation coefficient in cor and
#p-value for the correlation in pval
results <- expand.grid(Signature = unlist(signature_scores), Gene = slit_genes)
results$cor <- NA
results$pval <- NA

#Looping to compute Pearson correlations
for (i in 1:nrow(results)) {
  #Variable g represents current SLIT gene
  g <- results$Gene[i]
  #Variable s represents current CAF signature
  s <- results$Signature[i]
  #Argument expr_data[[g]] represents numeric vector of that SLIT gene expression per-cell
  #Argument expr_data[[s]] represents numeric vector of CAF signature score per-cell
  test <- cor.test(expr_data[[g]], expr_data[[s]], method = "pearson")
  #Assigning $estimate (r value) to the previously introduced cor column
  results$cor[i] <- test$estimate
  results$pval[i] <- test$p.value
}


#Creating correlation matrix (CAF Signatures as rows, SLIT genes as columns)
#Converting the long-format table (results) into a matrix, because that is the format required
#by the pheatmap()  
heatmap_matrix <- reshape2::acast(results, Signature ~ Gene, value.var = "cor")


#Adding asterisks for p-values
#Creating an empty matrix of the same dimensions as the correlation matrix
asterisks <- matrix("", nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))

#Assigning matching row/column names so they align perfectly when plotted
rownames(asterisks) <- rownames(heatmap_matrix)
colnames(asterisks) <- colnames(heatmap_matrix)

for (i in 1:nrow(results)) {
  s <- results$Signature[i]
  g <- results$Gene[i]
  p <- results$pval[i]
  if (p < 0.001) asterisks[s, g] <- "***"
  else if (p < 0.01) asterisks[s, g] <- "**"
  else if (p < 0.05) asterisks[s, g] <- "*"
}


#Plotting heatmap with correct axes: rows = signatures, columns = SLIT genes
pheatmap(heatmap_matrix,
         display_numbers = asterisks,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Fibroblasts LS: SLIT2/3 genes and CAF Signatures")


#Checking how SLIT2 and SLIT3 genes expressions are distributed in space
FeaturePlot(seurat_LS_fibro, features = c("SLIT2", "SLIT3"))

VlnPlot(seurat_LS_fibro, features = c("SLIT2", "SLIT3"))

SpatialFeaturePlot(seurat_LS_fibro, features = c("SLIT2", "SLIT3"))