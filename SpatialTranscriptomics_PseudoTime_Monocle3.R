install.packages("devtools")
remotes::install_github("bnprks/BPCells", subdir = "r")

devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")

install.packages("R.utils")

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)        
library(readxl)


set.seed(1234)
#Loading in the subsetted data according to survival (LS) and SLIT expression
seu <- readRDS("LS_CAFs_SLITClassified.rds")

DefaultAssay(seu) <- "Spatial"

seu

unique(seu@meta.data$seurat_clusters)
unique(seu@meta.data$orig.ident)
unique(seu@meta.data$CompositionCluster_CC)
unique(seu@meta.data$Spatial_snn_res.0.4)
unique(seu@meta.data$Spatial_snn_res.0.8)
unique(seu@meta.data$Spatial_snn_res.0.6)
unique(seu@active.ident)

unique(seu@meta.data$SLIT_Class)
table(seu$SLIT_Class)

a1 <- DimPlot(seu, reduction = 'umap', group.by = 'SLIT_Class', label = T)
a2 <- DimPlot(seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2


#Monocle3 requires cell_data_set (cds) object
#Converting Seurat Object to cell_data_set object for monocle3
cds <- as.cell_data_set(seu)
cds

#Getting cell metadata
colData(cds)
#Getting gene metdata
fData(cds)
rownames(fData(cds))[1:10]
#Since it misses the gene_short_name column, adding it here
#Option 1
rowData(cds)$gene_short_name <- rownames(cds)

#Option 2 (fData() is featureData() and it's the equivalent to rowData())
fData(cds)$gene_short_name <- rownames(fData(cds))

#Counts
counts(cds)


#Clustering cells (using clustering info from Seurat's UMAP)
#Monocle3 needs clustering information (either with it's own function or to retrieve
#the clustering information previously obtained with Seurat pipeline). 
#Monocle3 also has functions to cluster cells, where it not only determines the clusters,
#but also partitions and these partitions are nothing but super clusters.
#Three things are needed here: adding partitions information, adding clustering information
#and UMAP cell embeddings which are UMAP coordinates.

#Assigning partitions (Adding partitions information)
#And asigning all the cells to the one partition by crating a named vector or a named list
#For the named list, the names should be cell IDs and the partition values will be 1, because
#all those cells will be assigned to the same partition
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition


#Assigning the cluster info (Adding clustering information from Seurat)
list_cluster <- seu@active.ident
cds@clusters$UMAP$clusters <- list_cluster


#Assigning UMAP coordinates - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seu@reductions$umap@cell.embeddings


#Plot before learning the trajectory
#In Monocle3 plot() cells is the equivalent of DimPlot(), thus allowing
#UMAP generation and visualization of cells
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")


cluster.names <- plot_cells(cds,
                            color_cells_by = "SLIT_Class",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")


cluster.before.trajectory | cluster.names


#Learning trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'SLIT_Class',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


#Ordering the cells in pseudotime
#Defining root cells from SLIT_Class 
root_label <- "SLIT2+_SLIT3+_high"

if (!root_label %in% unique(colData(cds)$SLIT_Class)) {
  stop(paste0("Requested root_label '", root_label, "' is not present in SLIT_Class."))
}

root_cells <- rownames(colData(cds))[colData(cds)$SLIT_Class == root_label]


cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)


#Plots
p_pt <- plot_cells(
  cds, color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE, label_branch_points = TRUE, label_leaves = TRUE
) + ggtitle(paste0("Monocle3 pseudotime (root = ", root_label, ")"))

p_slit <- plot_cells(
  cds, color_cells_by = "SLIT_Class",
  label_groups_by_cluster = FALSE, label_branch_points = TRUE
) + ggtitle("SLIT classes on trajectory")

print(p_pt); print(p_slit)


#Cells ordered by Monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
#Checking if monocle3_pseudotime column is created and present there
colData(cds)

data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, SLIT_Class, fill = SLIT_Class)) +
  geom_boxplot()
#Re-ordering the y-axis according to the median values of the pseudotime 
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(SLIT_Class, monocle3_pseudotime, median), fill = SLIT_Class)) +
  geom_boxplot()


#Finding genes that change expression as cells progress along a trajectory, that is
#finding genes that change as a function of pseudotime, by testing genes for
#differential expression based on the low dimension embeddings.
#This statistics tells whether the cells at nearby positions on a trajectory will have
#similiar or disimiliar expression levels for the genes that are being tested
deg_cds <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_cds %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

#Understanding and looking at these differentially expressed genes can give us an idea on the
#pathways or the regulatory mechanisms that are differentially regulated when certain group of
#cells take a shorter or longer trajectory.
#These DEGs can be leveraged to perform a pathway analysis or GSEA to understand what pathways
#are differentially expressed, especially if studying the mechanism of a response to certain
#therapy or understanding the mechanism of development of resistance to a therapy
FeaturePlot(seu, features = c('ISG15', 'AGRN', 'MXRA8', 'ERRFI1', 'CTRC', 'CELA2A'))


#Visualizing pseudotime in Seurat
seu$pseudotime <- pseudotime(cds)
FeaturePlot(seu, features = "pseudotime")
Idents(seu) <- seu$SLIT_Class
FeaturePlot(seu, features = "pseudotime", label = T)


#Plots
p_pt <- plot_cells(
  cds, color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE, label_branch_points = TRUE, label_leaves = TRUE
) + ggtitle(paste0("Monocle3 pseudotime (root = ", root_label, ")"))

p_slit <- plot_cells(
  cds, color_cells_by = "SLIT_Class",
  label_groups_by_cluster = FALSE, label_branch_points = TRUE
) + ggtitle("SLIT classes on trajectory")

print(p_pt); print(p_slit)


#Exporting pseudotime back to Seurat 
pt <- monocle3::pseudotime(cds)
seu <- AddMetaData(seu, metadata = pt, col.name = "pseudotime_SLIT")


#Summaries: median pseudotime per SLIT class to reveal an ordering
slit_order <- seu@meta.data %>%
  dplyr::group_by(SLIT_Class) %>%
  dplyr::summarize(n = dplyr::n(),
                   median_pt = median(pseudotime_SLIT, na.rm = TRUE),
                   mean_pt   = mean(pseudotime_SLIT, na.rm = TRUE)) %>%
  dplyr::arrange(median_pt)
print(slit_order)


#Boxplot of pseudotime by SLIT class
bp <- ggplot(seu@meta.data,
             aes(x = factor(SLIT_Class, levels = slit_order$SLIT_Class),
                 y = pseudotime_SLIT)) +
  geom_boxplot(outlier.size = 0.8) +
  xlab("SLIT class (ordered by median pseudotime)") +
  ylab("Pseudotime") +
  ggtitle("Pseudotime by SLIT class")
print(bp)


#Visualizing SLIT2 and SLIT3 along pseudotime 
#(expecting monotone trends if trajectory aligns with SLIT axis)
genes_to_plot <- intersect(c("SLIT2","SLIT3"), rownames(cds))
if (length(genes_to_plot) > 0) {
  for (g in genes_to_plot) {
    print(plot_genes_in_pseudotime(cds[genes_to_plot[genes_to_plot==g],], color_cells_by = "SLIT_Class") +
            ggtitle(paste("Expression of", g, "along pseudotime")))
  }
} else {
  message("SLIT2/SLIT3 not found in features for plot_genes_in_pseudotime.")
}


#Saving the outputs 
saveRDS(cds, file = "LS_CAFs_Monocle3_cds_SLIT.rds")
saveRDS(seu, file = "LS_CAFs_SLITClassified_withPseudotime.rds")
saveRDS(deg_cds, file = "LS_CAFs_SLIT_Pseudotime_DGE.rds")