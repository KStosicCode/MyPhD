library(Seurat)
library(readxl)
library(reshape2)
library(pheatmap)
library(dplyr)


#Loading the path to the Excel table with CAF signatures
signatures_table_path <- "C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Toulouse/gene_sets_Zeynep2.xlsx"


#Loading Seurat object
seurat <- readRDS("LS_CAFs_SLITClassified.rds")
DefaultAssay(seurat) 


#Creating a function for adding all module scores for each gene signature from Excel 
add_all_signatures <- function(seu, xlsx_path, min_genes = 2) {
  #Creating a blank list to store the names of the new columns in Seurat object
  #that will contain each module score
  score_cols_names <- list()
  #Looping over all the sheets in the Excel table
  for (sheet in readxl::excel_sheets(xlsx_path)) {
    #Assigning that particular (current) sheet to a variable "signatures_table" as a df
    signatures_table <- readxl::read_excel(xlsx_path, sheet = sheet)
    #Looping over each column in that data frame where sig_col is the name of the current column
    for (sig_col in colnames(signatures_table)) {
      
      genes <- signatures_table[[sig_col]] |> toupper() |> na.omit()
      genes <- genes[genes %in% rownames(seu)]
      
      #If fewer than min_genes remain after filtering, skip this signature entirely
      if (length(genes) < min_genes) next
      #Creating a unique name for each signature (sheet name_column name)
      sig_name <- paste(sheet, sig_col, sep = "_")
      #Assessing the expression level of genes in that signature for each cell individually
      #It will create a new column in meta.data named <sig_name>1
      seu <- AddModuleScore(seu, features = list(genes), name = sig_name)
      #Storing the new column name of this module score
      #Module score columns always have a “1” suffix (from AddModuleScore behavior).
      score_cols_names[[sig_name]] <- paste0(sig_name, "1")
    }
  }
  #In a function, variables created inside the function don’t exist outside unless explicitly 
  #returned. To use the updated Seurat object and the list of column names later, 
  #they must be returned.
  list(seu = seu, score_cols_names = unlist(score_cols_names, use.names = TRUE))
}


#Adding scores using the previously created function "add_all_signatures"
added <- add_all_signatures(seu = seurat, xlsx_path = signatures_table_path, min_genes = 2)

added$seu
seurat
#Now, variable "added" contains an object "seu" and a character vector "score_cols_names"
#Object "seu" is the same as "seurat", only upgraded with module scores for each signature
added$score_cols_names

seurat <- added$seu
signature_scores <- added$score_cols_names


#Fetching module scores and SLIT class labels
stopifnot("SLIT_Class" %in% colnames(seurat@meta.data))
vars <- c(signature_scores, "SLIT_Class")
df <- FetchData(seurat, vars = vars)


#Building 0/1 membership columns for each SLIT class because the next step is point-biserial  
#correlation (to assess the strength and direction of the relationship between one continuous  
#variable and one dichotomous (binary) variable).
#A categorical variable (SLIT class) can't directly correlate with a continuous variable (score)
#The standard trick is point-biserial correlation (just Pearson r where one variable is 0/1).
classes <- na.omit(unique(df$SLIT_Class))

for (cl in classes) {
  #Creating a binary column for each class (, i. e. generating new column for each class in
  #such a way that "1" means that it's that class and "0" means that it's one of the other 5)
  df[[paste0("is_", cl)]] <- as.integer(df$SLIT_Class == cl)
}


#Computing point-biserial correlations: r between each signature score and each is_<class>
#Function expand.grid() builds all pairs: every signature score × every class
results <- expand.grid(Signature = signature_scores, Class = paste0("is_", classes))
results$cor <- NA_real_
results$pval <- NA_real_

for (i in seq_len(nrow(results))) {
  #Variable s_col represents current CAF signature
  s_col <- results$Signature[i]
  #Variable c_col is a binary indicator representing current class from SLIT classification 
  c_col <- results$Class[i]
  #Extracting the two numeric vectors of interest to compare in the correlation test
  x <- df[[s_col]]
  y <- df[[c_col]]
  #Running cor.test(score, is_class) for each pair
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  results$cor[i]  <- unname(ct$estimate)
  results$pval[i] <- ct$p.value
}


#Making matrix: rows = signatures, cols = SLIT classes
mat <- acast(results, Signature ~ Class, value.var = "cor")


#Asterisks
aster <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                #Using dimnames() to copy both rownames and colnames from mat immediately,
                #otherwise they should be assigned after this step, but before the for loop as
                #rownames(aster) <- rownames(mat) and colnames(aster) <- colnames(mat)
                dimnames = dimnames(mat))

for (i in seq_len(nrow(results))) {
  s <- results$Signature[i]; c <- results$Class[i]; p <- results$pval[i]
  if (is.na(p)) next
  aster[s, c] <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}


#Dropping the "is_" prefix so columns read like actual class names.
colnames(mat)   <- sub("^is_", "", colnames(mat))
colnames(aster) <- colnames(mat)


#Heatmap showing correlation
pheatmap(mat,
         display_numbers = aster,
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue","white","red"))(100),
         main = "SLIT classes vs CAF signatures (point-biserial r)")