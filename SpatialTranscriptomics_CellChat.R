devtools::install_github("jinworks/CellChat")
library(Seurat)    
library(CellChat)
library(patchwork)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)


set.seed(1234)
#Loading the Seurat Object of interest
seurat_obj <- readRDS("LS_CAFs_PostHarmony.rds")
#seurat_obj <- readRDS("LS_PostHarmony.rds")


#Subsetting the main Seurat Object according to the sample origin
#groups <- seurat_obj$sample
seu1 <- subset(x = seurat_obj, subset = seurat_obj$sample == "PDAC4")
seu2 <- subset(x = seurat_obj, subset = seurat_obj$sample == "PDAC5")
seu3 <- subset(x = seurat_obj, subset = seurat_obj$sample == "PDAC7")
seu4 <- subset(x = seurat_obj, subset = seurat_obj$sample == "PDAC8")
#Checking if the subsetting was done well
dim(seurat_obj) # 18085  9390
dim(seu1) # 18085  2651
dim(seu2) # 18085  3491
dim(seu3) # 18085  2387
dim(seu4) # 18085   861

#Retrieving info needed for four key input arguments for createCellChat()
#The normalized expression matrix
data.input1 <- GetAssayData(seu1, assay = "Spatial", slot = "data") 
data.input2 <- GetAssayData(seu2, assay = "Spatial", slot = "data")
data.input3 <- GetAssayData(seu3, assay = "Spatial", slot = "data")
data.input4 <- GetAssayData(seu4, assay = "Spatial", slot = "data")

length(rownames(data.input1))
length(rownames(data.input2))
length(rownames(data.input3))
length(rownames(data.input4))

length(colnames(data.input1))
length(colnames(data.input2))
length(colnames(data.input3))
length(colnames(data.input4))

#No need to run this because colnames of each individual data.input are already preceded by
#the ordinal number of the samples. If this wasn't the case, PDACn_ and colnames(data.input_n)
#should get pasted
#colnames(data.input1) <- paste0("PDAC4_", colnames(data.input1))
#colnames(data.input2) <- paste0("PDAC5_", colnames(data.input2))
#colnames(data.input3) <- paste0("PDAC7_", colnames(data.input3))
#colnames(data.input4) <- paste0("PDAC8_", colnames(data.input4))

data.input_total <- cbind(data.input1, data.input2, data.input3, data.input4)

#Cell type labels from the Seurat object
#metadata1  <- data.frame(labels = seu1$CellType_Specific, 
                        #samples = seu1$sample, 
                        #row.names = colnames(seu1))

#metadata2  <- data.frame(labels = seu2$CellType_Specific, 
                         #samples = seu2$sample, 
                         #row.names = colnames(seu2))

#metadata3  <- data.frame(labels = seu3$CellType_Specific, 
                         #samples = seu3$sample, 
                         #row.names = colnames(seu3))

#metadata4  <- data.frame(labels = seu4$CellType_Specific, 
                         #samples = seu4$sample, 
                         #row.names = colnames(seu4))

metadata1  <- data.frame(labels = seu1$CellType_Specific, samples = seu1$sample)
metadata2  <- data.frame(labels = seu2$CellType_Specific, samples = seu2$sample)    
metadata3  <- data.frame(labels = seu3$CellType_Specific, samples = seu3$sample)
metadata4  <- data.frame(labels = seu4$CellType_Specific, samples = seu4$sample)

head(colnames(seu1))
head(colnames(data.input1))
head(colnames(seu1)) == head(colnames(data.input1)) # TRUE TRUE TRUE TRUE TRUE TRUE
tail(colnames(seu1)) == tail(colnames(data.input1)) # TRUE TRUE TRUE TRUE TRUE TRUE

metadata_total <- rbind(metadata1, metadata2, metadata3, metadata4)

head(rownames(metadata_total)) == head(colnames(data.input_total))
tail(rownames(metadata_total)) == tail(colnames(data.input_total))

#No need for rownames(metadata_total) <- colnames(data.input_total), because it's already 
#like that

#Checking if the columns in metadata_total are factored and factoring them if they aren't
is.factor(metadata_total$labels) # TRUE
is.factor(metadata_total$samples) #FALSE

metadata_total$samples <- factor(metadata_total$samples, 
                                 levels = c("PDAC4", "PDAC5", "PDAC7", "PDAC8"))
class(metadata_total$samples)

unique(metadata_total$labels)
unique(metadata_total$samples)

#Spatial coordinates/ locations of each cell/ spot centroid
head(seurat_obj@images$PDAC4$centroids@coords)
seu1@images$PDAC4$centroids@coords == seurat_obj@images$PDAC4$centroids@coords

head(seurat_obj@images$PDAC5$centroids@coords)
seu2@images$PDAC5$centroids@coords == seurat_obj@images$PDAC5$centroids@coords

head(seurat_obj@images$PDAC7$centroids@coords)
seu3@images$PDAC7$centroids@coords == seurat_obj@images$PDAC7$centroids@coords

head(seurat_obj@images$PDAC8$centroids@coords)
seu4@images$PDAC8$centroids@coords == seurat_obj@images$PDAC8$centroids@coords

#It's important to take each sample individually as an input for the coordinates
spatial.locs1 <- Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("x", "y"))
spatial.locs2 <- Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("x", "y"))
spatial.locs3 <- Seurat::GetTissueCoordinates(seu3, scale = NULL, cols = c("x", "y"))
spatial.locs4 <- Seurat::GetTissueCoordinates(seu4, scale = NULL, cols = c("x", "y"))

spatial.locs_total <- rbind(spatial.locs1, spatial.locs2, spatial.locs3, spatial.locs4)

spatial.locs_total <- spatial.locs_total[,1:2]

spatial.locs_total <- as.matrix(spatial.locs_total)

head(rownames(spatial.locs_total)) == head(colnames(data.input_total))
tail(rownames(spatial.locs_total)) == tail(colnames(data.input_total))

#No need for rownames(spatial.locs_total) <- colnames(data.input_total), because it's already
#set that way

#Spatial factors of spatial distance
#Doing it separately for each samples (PDAC4, PDAC5, PDAC7 and PDAC8)
#PDAC4
scalefactors1  = jsonlite::fromJSON(txt = file.path("C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Json files/wetransfer_scale-factors_2025-08-18_0954", 'PDAC4scalefactors_json.json'))
spot.size = 55 # the theoretical spot size (um) in 10X Visium v2.0
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)

#x <- spatial.locs1[,1:2]
#d.spatial <- computeCellDistance(coordinates = x, 
                                 #ratio = spatial.factors1$ratio, 
                                 #tol = spatial.factors1$tol)
#y <- spatial.locs2[,1:2]

#min(d.spatial[d.spatial!=0]) 

#Quick check of cell type labels 
head(metadata_total$labels)

#PDAC5
scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Json files/wetransfer_scale-factors_2025-08-18_0954", 'pdac5_scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)

#PDAC7
scalefactors3 = jsonlite::fromJSON(txt = file.path("C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Json files/wetransfer_scale-factors_2025-08-18_0954", 'pdac7_scalefactors_json.json'))
conversion.factor3 = spot.size/scalefactors3$spot_diameter_fullres
spatial.factors3 = data.frame(ratio = conversion.factor3, tol = spot.size/2)

#PDAC8
scalefactors4 = jsonlite::fromJSON(txt = file.path("C:/Users/Admin/Desktop/Spatial Transcriptomics Data/Json files/wetransfer_scale-factors_2025-08-18_0954", 'pdac8_scalefactors_json.json'))
conversion.factor4 = spot.size/scalefactors4$spot_diameter_fullres
spatial.factors4 = data.frame(ratio = conversion.factor4, tol = spot.size/2)
  
#Binding them altogether 
spatial.factors_total <- rbind(spatial.factors1, 
                         spatial.factors2,
                         spatial.factors3, 
                         spatial.factors4)

rownames(spatial.factors_total) <- c("PDAC4", "PDAC5","PDAC7","PDAC8")


#Creating a CellChat object from the data, grouping cells by the specified cell type label
cellchat <- createCellChat(object = data.input_total, 
                           meta = metadata_total,
                           group.by = "labels",
                           datatype = "spatial",
                           coordinates = spatial.locs_total,
                           spatial.factors = spatial.factors_total)  


#Setting the CellChat ligand-receptor database to the human dataset (for Homo sapiens)
#Loading the default human CellChat database of interactions
CellChatDB <- CellChatDB.human      

#Using a subset of CellChatDB for cell-cell communication analysis
#Using Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

#Assigning the database to the CellChat object
cellchat@DB <- CellChatDB.use

#(Optional) Subsetting the database to a category if needed, for example only secreted signaling
#cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")


#Pre-processing: filtering data and identifying over-expressed ligands or receptors in 
#one cell group and then identifying over-expressed ligand-receptor interactions if either 
#ligand or receptor is over-expressed.
#If Subsetting was done in the previous step, expression data here should also get subsetted  
#to signaling genes for efficiency. But even if no subsetting was done in the previous step,
#this step is necessary even if using the whole database
cellchat <- subsetData(cellchat)   
future::plan("multisession", workers = 4) 
#Finding genes in each group that are highly expressed
cellchat <- identifyOverExpressedGenes(cellchat)  
#Predicting over-expressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)  


#Inferring cell-cell communication network
#Calculating communication probability between all cell pairs
#cellchat <- computeCommunProb(cellchat, 
                              #type = "triMean",
                              #trim = 0.1,
                              #distance.use = FALSE, 
                              #interaction.range = 250, 
                              #contact.knn = TRUE, 
                              #contact.knn.k = 5)

library(future)
plan(multisession, workers = 4)             
options(future.globals.maxSize = 8 * 1024^3)

cellchat <- computeCommunProb(cellchat,
                              type = "truncatedMean",
                              trim = 0.1,
                              distance.use = TRUE,
                              interaction.range = 250,
                              scale.distance = 0.01,
                              contact.dependent = TRUE,
                              contact.range = 100,
                              contact.knn.k = 5)

#Saving the output of computeCommunProb() 
saveRDS(cellchat, file = "Cellchat_LS_CAF_CommunProb.rds")


#Filtering out weak links from groups with few cells
#Argument min.cells=10 ensures having enough cells per group to consider 
#communications significant
View(Cellchat_LS_CAF_CommunProb)
cellchat <- Cellchat_LS_CAF_CommunProb
cellchat <- filterCommunication(cellchat, min.cells = 10)  


#Aggregating signals at the pathway level and create an overall network
#Summing communication probabilities by signaling pathway
cellchat <- computeCommunProbPathway(cellchat) 
head(cellchat@net)
head(cellchat@netP)
#Aggregating the network for all interactions
#cellchat@net will contain the inferred number and strength of interactions between 
#each pair of cell groups
cellchat <- aggregateNet(cellchat)                    


#Visualizing the global communication network among all cell types
#Number of cells (spots) in each cell type group
groupSize <- as.numeric(table(cellchat@idents))   
#Plotting side by side
par(mfrow=c(1,2), xpd=TRUE)       
#Circle plot of interaction counts
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions") 
#Circle plot of interaction weights
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction strength")  
#Heatmap 
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

#Upon inferring the cell-cell communication network, CellChat provides various functionality for
#further data exploration, analysis, and visualization. For example, the `circle plot` 
#and the new `spatial plot`.
#Visualization of cell-cell communication at different levels: One can visualize the inferred
#communication network of signaling pathways using netVisual_aggregate` and visualize the 
#inferred communication networks of individual L-R pairs associated with that signaling pathway 
#using netVisual_individual.
#All the signaling pathways showing significant communications can be accessed by
#cellchat@netP$pathways.
#Circle plot for particular pathways of interest
cellchat@netP$pathways
#Single pathway 
pathways.show <- c("TGFb")
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#Three pathways simultaneously 
#First three
sigs <- c("TGFb","CXCL","CCL")

op <- par(mfrow = c(1, length(sigs)), xpd = TRUE, mar = c(1,1,2,1))

for (sg in sigs) {
  netVisual_aggregate(cellchat, signaling = sg, layout = "circle")
}
par(op)

#Second three
sigs <- c("PDGF","FGF","CSF")

op <- par(mfrow = c(1, length(sigs)), xpd = TRUE, mar = c(1,1,2,1))

for (sg in sigs) {
  netVisual_aggregate(cellchat, signaling = sg, layout = "circle")
}
par(op)

#Third three
sigs <- c("MIF","PERIOSTIN","GALECTIN")

op <- par(mfrow = c(1, length(sigs)), xpd = TRUE, mar = c(1,1,2,1))

for (sg in sigs) {
  netVisual_aggregate(cellchat, signaling = sg, layout = "circle")
}
par(op)

#Fourth three
sigs <- c("IGF","IGFBP","GAS")

op <- par(mfrow = c(1, length(sigs)), xpd = TRUE, mar = c(1,1,2,1))

for (sg in sigs) {
  netVisual_aggregate(cellchat, signaling = sg, layout = "circle")
}
par(op)

#All 12 pathways at once
sigs <- c("TGFb","CXCL","CCL","PDGF","FGF","CSF","MIF",
          "PERIOSTIN","GALECTIN","IGF","IGFBP","GAS")

#Deciding grid size: e.g. 3 rows × 4 columns (since there are 12 pathways)
op <- par(mfrow = c(3, 4), xpd = TRUE, mar = c(1,1,2,1))

for (sg in sigs) {
  netVisual_aggregate(cellchat, signaling = sg, layout = "circle")
}

par(op) 

#Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "PDAC4", 
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 0)

netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "PDAC5", 
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 0)

netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "PDAC7", 
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 0)

netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "PDAC8", 
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 0)



#Checking out what L-R interactions are the significant ones in case of the listed pathways
comm_tgfb <- subsetCommunication(cellchat, signaling = "TGFb")
head(comm_tgfb)
#Picking the top pair (based on communication probability)
top_pair <- comm_tgfb[which.max(comm_tgfb$prob), ]
#Plotting that top pair (netVisual_individual)
netVisual_individual(cellchat,
                     signaling = "TGFb",
                     pairLR.use = top_pair$interaction_name,
                     layout = "circle")

top3 <- head(comm_tgfb[order(-comm_tgfb$prob), ], 3)

netVisual_individual(cellchat,
                     signaling = "TGFb",
                     pairLR.use = top3$interaction_name,
                     layout = "circle")


#Computing and visualizing the network centrality scores
#Computing centrality to identify dominant senders (sources) and receivers (targets) 
#in the network
#Analyzing network centrality for pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")   
#Scatter plot: outgoing vs incoming communication role for each cell type
netAnalysis_signalingRole_scatter(cellchat)           

#Visualizing the computed centrality scores using heatmap, allowing ready identification of 
#major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, 
                                  height = 2.5, font.size = 10)
#Spatial plots
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "PDAC4", 
                    layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                    vertex.weight = "incoming", vertex.size.max = 6, vertex.label.cex = 0)


#Focusing analysis on fibroblasts (CAFs)
#Examining which interactions fibroblasts send and receive.
#Index of the fibroblast group
#fibro_ind <- which(levels(cellchat@idents) == "iCAF" & levels(cellchat@idents) == "myCAF")  
#Bubble plot of all significant ligand-receptor interactions from fibroblasts 
#to other cell types
#netVisual_bubble(cellchat, sources.use = fibro_ind, 
                 #targets.use = setdiff(1:length(levels(cellchat@idents)), fibro_ind), 
                 #remove.isolate = FALSE, title.name = "Fibroblasts as signal senders")
#Bubble plot of all significant interactions from all other cell types targeting fibroblasts 
#netVisual_bubble(cellchat, sources.use = setdiff(1:length(levels(cellchat@idents)), fibro_ind), 
                 #targets.use = fibro_ind, 
                 #remove.isolate = FALSE, title.name = "Fibroblasts as signal receivers")


#Identifying key signaling pathways involving fibroblasts
#Listing the top pathways that fibroblasts send or receive
#All communications sent by fibroblasts
#fibro_out <- subsetCommunication(cellchat, sources.use = fibro_ind)   
#All communications received by fibroblasts
#fibro_in  <- subsetCommunication(cellchat, targets.use = fibro_ind)   
#Printing a few top ligand-receptor pairs sent by fibroblasts
#head(fibro_out$interaction_name) 
#Printing a few top pairs received by fibroblasts
#head(fibro_in$interaction_name)   


#Computing the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)


#Spatial Feature Plots
#TGFB3_ACVR1_TGFBR1
#PDAC4
#Taking an input of a few genes
spatialFeaturePlot(cellchat, features = c("TGFB3","ACVR1","TGFBR1"), sample.use = "PDAC4", 
                   point.size = 0.8, color.heatmap = "Reds", direction = 1)

#Taking an input of a ligand-receptor pair
spatialFeaturePlot(cellchat, pairLR.use = "TGFB3_ACVR1_TGFBR1", sample.use = "PDAC4", 
                   point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
                   
#Taking an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "TGFB3_ACVR1_TGFBR1", sample.use = "PDAC4", 
                   point.size = 1.5, do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1)

#PDAC5, PDAC7 and PDAC8


#TGFB1_ACVR1_TGFBR1
spatialFeaturePlot(cellchat, features = c("TGFB1","ACVR1","TGFBR1"), sample.use = "PDAC4", 
                   point.size = 0.8, color.heatmap = "Reds", direction = 1)



#Saving the CellChat object 
saveRDS(cellchat, file = "Cellchat_LS_CAF.rds")

x <- cellchat
seurat <- LS_CAFs_SLITClassified
all(colnames(seurat) %in% rownames(cellchat@meta))
head(colnames(seurat) %in% rownames(cellchat@meta))
tail(colnames(seurat) %in% rownames(cellchat@meta))

#Adding SLIT class to CellChat meta
all(colnames(seurat) %in% rownames(x@meta))
x@meta$SLIT_Class <- seurat$SLIT_Class[match(rownames(cellchat@meta), colnames(seurat))]
x@meta$SLIT_Class

#All TGFb pairs that passed the filters
comm_tgfb <- subsetCommunication(cellchat, signaling = "TGFb")

#Ligands are already gene symbols
lig <- unique(comm_tgfb$ligand)

#Receptors can be multi-subunit, like "ACVR1_TGFBR1" or "TGFBR1_R2"
#Splitting them on "_" and collecting unique parts
rec_parts <- unique(unlist(str_split(comm_tgfb$receptor, "_")))

#Full gene list for plotting
tgfb_genes <- unique(c(lig, rec_parts))
tgfb_genes <- c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "R2", "ACVR1B", "TGFBR2", "ACVR1", "TGFBR")

c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "R2", "ACVR1B", "TGFBR2", "ACVR1", "TGFBR") %in% rownames(seurat)
#TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE

#Keeping only genes present in Seurat assay
tgfb_genes <- tgfb_genes[tgfb_genes %in% rownames(seurat)]
length(tgfb_genes); tgfb_genes


VlnPlot(seurat, features = tgfb_genes, group.by = "SLIT_Class", pt.size = 0.2, stack = FALSE)