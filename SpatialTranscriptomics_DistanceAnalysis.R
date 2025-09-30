install.packages("STDistance")
library(devtools)
devtools::install_github("PrinceWang2018/ST_Distance")
library(STDistance)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(Hmisc)
library(scales)
library(stats)
library(RColorBrewer)
library(tidyr)
library(stringr)


#Loading Seurat object with Class labels
seurat_LS <- readRDS("LS_Classified_CAFs_Tumor_Immune.rds")
x <- seurat_LS


#Replacing "-" with "neg" and "+" with "pos" in the Class column
seurat_LS$Class <- gsub("-", "neg", seurat_LS$Class)
seurat_LS$Class <- gsub("\\+", "pos", seurat_LS$Class)

#Checking the result
table(seurat_LS$Class)


#Spatial coordinates files
#Loading tissue_positions.csv for LS samples
setwd("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Toulouse Fibroblasts/Distance Analysis/wetransfer_tissue_positions_pdac3-csv_2025-09-29_0747")

pdac4 <- read.csv("tissue_positions_PDAC4.csv", header = TRUE)
pdac4$Sample <- "PDAC4"
pdac4$Newbarcode <- paste(pdac4$Sample, pdac4$barcode, sep = "_")

pdac5 <- read.csv("tissue_positions_PDAC5.csv", header = TRUE)
pdac5$Sample <- "PDAC5"
pdac5$Newbarcode <- paste(pdac5$Sample, pdac5$barcode, sep = "_")

pdac7 <- read.csv("tissue_positions_PDAC7.csv", header = TRUE)
pdac7$Sample <- "PDAC7"
pdac7$Newbarcode <- paste(pdac7$Sample, pdac7$barcode, sep = "_")

pdac8 <- read.csv("tissue_positions_PDAC8.csv", header = TRUE)
pdac8$Sample <- "PDAC8"
pdac8$Newbarcode <- paste(pdac8$Sample, pdac8$barcode, sep = "_")


setwd("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Seurat Objects For Spatial Transcriptomics")
#Merging into one LS positions file
ls_positions <- rbind(pdac4, pdac5, pdac7, pdac8)


#Building a metadata data.frame with rownames as the join key
metadata <- seurat_LS@meta.data


#Normalizing spatial coordinates 
ls_positions_normalized <- normalize_spatial(ls_positions)


#Merging metadata and spatial information
posi <- merge(
  x = ls_positions_normalized,
  y = metadata,
  by.x = "Newbarcode",
  by.y = "row.names",
  all.y = TRUE
)


#Saving the merged posi object
saveRDS(posi, file = "Distance_Analysis_LS_posi.rds")


table(posi$Class)
names(table(posi$Class))

#Making sure NA values don't interfere with calculating distances 
posi$Class[is.na(posi$Class)]='unclassified'

#Calculating nearest distances between cell types
#Distances from ROBO1+
distance_results_plus <- calculate_nearest_distances(
  spatial_data   = posi,
  reference_type = "ROBO1pos",
  target_types   = c("Immune_cells", "SLIT2neg_SLIT3neg", "SLIT2neg_SLIT3pos", 
                     "SLIT2pos_SLIT3neg", "SLIT2pos_SLIT3pos_high", "SLIT2pos_SLIT3pos_low"),
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  id_col   = "Newbarcode",
  type_col = "Class"
)

#Saving the merged posi object
saveRDS(distance_results_plus, file = "Distance_Analysis_LS_ROBO1+.rds")

#Comparing distance among subgroups
#With log10
plot_distance_boxplot(
  distance_results_plus,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "log10",
  palette = "Dark2"
)
#Without log10
plot_distance_boxplot(
  distance_results_plus,
  id_col = "Newbarcode",
  show_points = TRUE,
  palette = "Dark2"
)

#Creating radial network visualization 
plot_radial_distance(
  distance_results_plus,
  id_col = "Newbarcode",
  reference_type = "ROBO1pos",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
)

#Density plot
require(reshape2)
m = melt(distance_results_plus[,2:7])

library(dplyr)
ma <- m %>%
  group_by(variable) %>%
  dplyr::summarise(mean=mean(value))

More1_Den <- ggplot(data = m, mapping = aes(x = value, fill = variable)) + 
  geom_line(stat = "density")+geom_density(alpha = 0.5, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#008B00", "#CDC9C9","#A020F0", 
                             "#FF0000", "#00FFFF", "#FFD700"))+
  labs(x = "Distance", y = "Density", title = "Distance from ROBO1+ in LS") +
  theme_minimal()+geom_vline(data = ma, aes(xintercept = mean, 
                                            color = variable), linetype="dashed", size=1)+ 
  theme(legend.position="bottom")

More1_Den

#Spatially visualizing interactions
#Between two cell types
visualize_spatial_network(
  posi,
  sample = "PDAC4",
  reference_type = "ROBO1pos",
  target_type = "SLIT2pos_SLIT3neg",
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  type_col = "Class",
  color_palette = c("ROBO1pos" = "#008B00", "SLIT2pos_SLIT3neg" = "#00FFFF"),
  alpha = 0.7
)  

#Between reference and multiple target types
visualize_spatial_multinetwork(
  posi,
  sample = "PDAC4",
  reference_type = "ROBO1pos",
  target_type = c("Immune_cells", "SLIT2neg_SLIT3neg", "SLIT2neg_SLIT3pos", 
                  "SLIT2pos_SLIT3neg", "SLIT2pos_SLIT3pos_high", "SLIT2pos_SLIT3pos_low"),
  type_col = "Class",
  color_palette = c("ROBO1pos" = "#008B00", 
                    "Immune_cells" = "#CDC9C9",
                    "SLIT2neg_SLIT3neg" = "#A020F0",
                    "SLIT2neg_SLIT3pos" = "#E41A1C",
                    "SLIT2pos_SLIT3neg" = "#00FFFF",
                    "SLIT2pos_SLIT3pos_high" = "#FFD700",
                    "SLIT2pos_SLIT3pos_low" = "#FF1493"),
  point_alpha = 0.7
)


#Distances from ROBO1-
distance_results_minus <- calculate_nearest_distances(
  spatial_data   = posi,
  reference_type = "ROBO1neg",
  target_types   = c("Immune_cells", "SLIT2neg_SLIT3neg", "SLIT2neg_SLIT3pos", 
                     "SLIT2pos_SLIT3neg", "SLIT2pos_SLIT3pos_high", "SLIT2pos_SLIT3pos_low"),
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  id_col   = "Newbarcode",
  type_col = "Class"
)

#Saving the merged posi object
saveRDS(distance_results_minus, file = "Distance_Analysis_LS_ROBO1-.rds")

#Comparing distance among subgroups
#With log10
plot_distance_boxplot(
  distance_results_minus,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "log10",
  palette = "Dark2"
)
#Without log10
plot_distance_boxplot(
  distance_results_minus,
  id_col = "Newbarcode",
  show_points = TRUE,
  palette = "Dark2"
)

#Creating radial network visualization 
plot_radial_distance(
  distance_results_minus,
  id_col = "Newbarcode",
  reference_type = "ROBO1neg",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
)

#Density plot
require(reshape2)
m2 = melt(distance_results_minus[,2:7])

library(dplyr)
ma2 <- m2 %>%
  group_by(variable) %>%
  dplyr::summarise(mean=mean(value))

More2_Den <- ggplot(data = m2, mapping = aes(x = value, fill = variable)) + 
  geom_line(stat = "density")+geom_density(alpha = 0.5, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#008B00", "#CDC9C9","#A020F0", 
                             "#FF0000", "#00FFFF", "#FFD700"))+
  labs(x = "Distance", y = "Density", title = "Distance from ROBO1- in LS") +
  theme_minimal()+geom_vline(data = ma2, aes(xintercept = mean, 
                                            color = variable), linetype="dashed", size=1)+ 
  theme(legend.position="bottom")

More2_Den

#Spatially visualizing interactions
#Between two cell types
visualize_spatial_network(
  posi,
  sample = "PDAC4",
  reference_type = "ROBO1neg",
  target_type = "SLIT2pos_SLIT3neg",
  x_col = "pxl_row_in_fullres",
  y_col = "pxl_col_in_fullres",
  type_col = "Class",
  color_palette = c("ROBO1neg" = "#008B00", "SLIT2pos_SLIT3neg" = "#00FFFF"),
  alpha = 0.7
)  

#Between reference and multiple target types
visualize_spatial_multinetwork(
  posi,
  sample = "PDAC4",
  reference_type = "ROBO1neg",
  target_type = c("Immune_cells", "SLIT2neg_SLIT3neg", "SLIT2neg_SLIT3pos", 
                  "SLIT2pos_SLIT3neg", "SLIT2pos_SLIT3pos_high", "SLIT2pos_SLIT3pos_low"),
  type_col = "Class",
  color_palette = c("ROBO1neg" = "#008B00", 
                    "Immune_cells" = "#CDC9C9",
                    "SLIT2neg_SLIT3neg" = "#A020F0",
                    "SLIT2neg_SLIT3pos" = "#E41A1C",
                    "SLIT2pos_SLIT3neg" = "#00FFFF",
                    "SLIT2pos_SLIT3pos_high" = "#FFD700",
                    "SLIT2pos_SLIT3pos_low" = "#FF1493"),
  point_alpha = 0.7
)