library(Seurat) 
library(CellChat)
library(patchwork)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)


set.seed(1234)


#Loading CellChat object of each dataset and then merge together
#What is needed is to run CellChat on each dataset seperately and then merge different 
#CellChat objects together.

#Loading in each CellChat output 
cellchat.SS <- Cellchat_SS_CAF
cellchat.LS <- Cellchat_LS_CAF

#cellchat.SS <- updateCellChat(cellchat.SS)
#cellchat.LS <- updateCellChat(cellchat.LS)

object.list <- list(SS = cellchat.SS, LS = cellchat.LS)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

levels(cellchat@meta$datasets)


#Predicting general principles of cell-cell communication
#CellChat starts with the big picture to predict general principles of cell-cell communication.
#When comparing cell-cell communication among multiple biological conditions, it can answer the
#following biological questions: 
#* Whether the cell-cell communication is enhanced or not
#* The interaction between which cell types is significantly changed
#* How the major sources and targets change from one condition to another
#Comparing the total number of interactions and interaction strength
#To answer whether the cell-cell communication is enhanced or not, CellChat compares the total
#number of interactions and interaction strength of the inferred cell-cell communication networks
#from different biological conditions. 

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


#Comparing the number of interactions and interaction strength among different cell populations.
#To identify the interaction between which cell populations showing significant changes, 
#CellChat compares the number of interactions and interaction strength among different 
#cell populations.
#Differential number of interactions or interaction strength among different cell populations;
#The differential number of interactions or interaction strength in the cell-cell communication
#network between two datasets can be visualized using circle plot, 
#where $\color{red}{\text{red}}$ (or $\color{blue}{\text{blue}}$) colored edges represent 
#$\color{red}{\text{increased}}$ (or $\color{blue}{\text{decreased}}$) signaling in the 
#second dataset compared to the first one (the second dataset is the second element of the 
#object list (LS), thus Red edges = interactions/signaling that are increased in LS relative 
#to SS, whereas Blue edges = interactions/signaling that are decreased in LS relative to SS)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#We can also show differential number of interactions or interaction strength 
#in a greater details using a heatmap. The top colored bar plot represents the sum of column
#of values displayed in the heatmap (incoming signaling). The right colored bar plot represents
#the sum of row of values (outgoing signaling). In the colorbar, 
#$\color{red}{\text{red}}$ (or $\color{blue}{\text{blue}}$) represents 
#$\color{red}{\text{increased}}$ (or $\color{blue}{\text{decreased}}$) signaling in the 
#second dataset compared to the first one (the second dataset is the second element of the 
#object list (LS), thus Red edges = interactions/signaling that are increased in LS relative 
#to SS, whereas Blue edges = interactions/signaling that are decreased in LS relative to SS)

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


#The differential network analysis only works for pairwise datasets. If there are more datasets
#for comparison, we can directly show the number of interactions or interaction strength between
#any two cell populations in each dataset. 
#To better control the node size and edge weights of the inferred networks across 
#different datasets, we compute the maximum number of cells per cell group and the maximum number
#of interactions (or interaction weights) across all datasets. 

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#Differential number of interactions or interaction strength among different cell types
#To simplify the complicated network and gain insights into the cell-cell communication at the 
#cell type level, we can aggregate the cell-cell communication based on the defined cell groups.
#Categorization of the cell populations into two cell types, and then re-merge the list of 
#CellChat object. 
levels(object.list$LS@idents)
levels(object.list$SS@idents)
#Object list already exists as object.list <- list(SS = cellchat.SS, LS = cellchat.LS)

#This is the part where in case of existence of multiple subclusters (for example, there are 4
#subcluster of fibroblasts, 4 of dendritic cells and 4 of T cells), they should all be categorized
#into major clusters (so 4 FIB subclusters become 1 FIB cluster, 4 DC become DC and 4 TC become
#TC). But here there are only iCAFs and myCAFs and all cell already belong to one of those two.
#Building the grouping vector for each object (order must match levels(x@idents))
#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
#group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
make_grouping <- function(x) {
  lev <- levels(x@idents)           # should be c("iCAF","myCAF")
  grp <- factor(lev, levels = c("iCAF","myCAF"))  # map each level to itself
  grp
}

#object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
object.list <- lapply(object.list, function(x) mergeInteractions(x, make_grouping(x)))

#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

#We then can show the number of interactions or interaction strength between any two cell types 
#in each dataset. 
#Comparing number/strength at the iCAF–myCAF level for each dataset
weight.max <- getMaxWeight(object.list, slot.name = c("idents","net","net"),
                           attribute = c("idents","count","count.merged"))

#Circle Plots for SS and LS showing numbers of interactions explicitly labeling the count of
#significant ligand–receptor pairs
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, 
                   weight.scale = T, label.edge= T, 
                   edge.weight.max = weight.max[3], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Circle Plots for SS and LS showing interaction strength and labeling the aggregated 
#communication probability
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, 
                   weight.scale = T, label.edge= T, 
                   edge.weight.max = weight.max[3], edge.width.max = 12, 
                   title.name = paste0("Interaction strength - ", names(object.list)[i]))
}


#Similarly, we can also show the differential number of interactions or interaction strength 
#between any two cell types using circle plot. Red (or blue) colored edges represent increased 
#(or decreased) signaling in the second dataset compared to the first one. 
par(mfrow = c(1,2), xpd=TRUE)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


#Comparing the major sources and targets in 2D space
#Comparing the outgoing and incoming interaction strength in 2D space allows ready identification
#of the cell populations with significant changes in sending or receiving signals between 
#different datasets. 
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

gg <- list()

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)


#From the scatter plot, we can see that Inflam.DC and cDC1 emerge as one of the major source and 
#targets in LS compared to NL. Fibroblast populations also become the major sources in LS. 
#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and
#LS. 
#Identifying signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))


#Identifying the conserved and context-specific signaling pathways
#CellChat then can identify signaling networks with larger (or less) difference, signaling groups,
#and the conserved and context-specific signaling pathways based on their cell-cell communication
#networks among multiple biological conditions. 

#Identifying signaling networks with larger (or less) difference as well as signaling groups based
#on their functional/structure similarity

#CellChat performs joint manifold learning and classification of the inferred communication
#networks based on their functional and topological similarity. 
#NB: Such analysis is applicable to more than two datasets. 

#**Functional similarity**: High degree of functional similarity indicates major senders and 
#*receivers are similar, and it can be interpreted as the two signaling pathways or two 
#*ligand-receptor pairs exhibit similar and/or redundant roles. 
#***NB**: Functional similarity analysis is not applicable to multiple datsets with 
#*different cell type composition. 

#**Structural similarity**: A structural similarity was used to compare their signaling network 
#*structure, without considering the similarity of senders and receivers. 
#***NB**: Structural similarity analysis is applicable to multiple datsets with the same cell type composition or the vastly different cell type composition. 

#Here we can run the manifold and classification learning analysis based on the functional 
#similarity because the two datasets have the the same cell type composition. 


#Identifying signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")

#Visualizing in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


#Identifying signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")

cellchat <- netEmbedding(cellchat, type = "structural")

cellchat <- netClustering(cellchat, type = "structural")

#Visualizing in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)


#Computing and visualizing the pathway distance in the learned joint manifold
#We can identify the signaling networks with larger (or less) difference based on their Euclidean
#distance in the shared two-dimensions space. Larger distance implies larger difference of the 
#communication networks between two datasets in terms of either functional or structure 
#similarity. **NB**: We only compute the distance of overlapped signaling pathways between two 
#datasets. Those signaling pathways that are only identified in one dataset are not considered 
#here. If there are more than three datasets, one can do pairwise comparisons by defining 
#`comparison` in the function `rankSimilarity`. 
rankSimilarity(cellchat, type = "functional")


#Identifying and visualize the conserved and context-specific signaling pathways
#By comparing the information flow/interaction strengh of each signaling pathway, we can identify
#signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their
#information flow at one condition as compared to another condition. 

#Comparing the overall information flow of each signaling pathway
#We can identify the conserved and context-specific signaling pathways by simply comparing the 
#information flow for each signaling pathway, which is defined by the sum of communication 
#probability among all pairs of cell groups in the inferred network (i.e., the total weights in 
#the network). 

#This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were 
#ranked based on differences in the overall information flow within the inferred networks between
#NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin. 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

gg1 + gg2


#Comparing outgoing (or incoming) signaling associated with each cell population
#The above analysis summarize the information from the outgoing and incoming signaling together. 
#We can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing
#to identify signaling pathways/ligand-receptors that exhibit different signaling patterns. 

#We can combine all the identified signaling pathways from different datasets and thus compare 
#them side by side, including outgoing signaling, incoming signaling and overall signaling by 
#aggregating outgoing and incoming signaling together. NB: `rankNet` also shows the comparison of
#overall signaling, but it does not show the signaling strength in specific cell populations. 
library(ComplexHeatmap)

i = 1

#Combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i],
                                        width = 5, height = 6)

ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 5, height = 6, color.heatmap = "GnBu")

ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6, color.heatmap = "GnBu")

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 5, height = 6, color.heatmap = "OrRd")

ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 5, height = 6, color.heatmap = "OrRd")

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Identifying the upgulated and down-regulated signaling ligand-receptor pairs

#Identifying dysfunctional signaling by comparing the communication probabities  
#We can compare the communication probabilities mediated by ligand-receptor pairs from some cell 
#groups to other cell groups. This can be done by setting `comparison` in the function 
#`netVisual_bubble`. 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), 
                 angle.x = 45)


#Moreover, we can identify the upgulated (increased) and down-regulated (decreased) signaling 
#ligand-receptor pairs in one dataset compared to the other dataset. This can be done by 
#specifying `max.dataset` and `min.dataset` in the function `netVisual_bubble`. 
#The increased signaling means these signaling have higher communication probability (strength) 
#in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, 
                        remove.isolate = T)

gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, 
                        remove.isolate = T)

gg1 + gg2

#NB: The ligand-receptor pairs shown in the bubble plot can be accessed via `signaling.
#LSIncreased = gg1$data`. 


#Identifying dysfunctional signaling by using differential expression analysis
#The above method for identifying the upgulated and down-regulated signaling is performed by 
#comparing the communication probability between two datasets for each L-R pair and each pair of
#cell groups. Alternative, we can identify the upgulated and down-regulated signaling 
#ligand-receptor pairs based on the differential gene expression analysis. Specifically, 
#we perform differential expression analysis between two biological conditions 
#(i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated 
#signaling based on the fold change of ligands in the sender cells and receptors in the receiver 
#cells. Such analysis can be done as follows.

#Defining a positive dataset, i.e., the dataset with positive fold change against the other 
#dataset
pos.dataset = "LS"

#Defining a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#Performing differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, 
                                       thresh.p = 1)

#Mapping the results of differential expression analysis onto the inferred cell-cell 
#communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extracting the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, 
                              receptor.logFC = NULL)

#Extracting the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, 
#i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, 
                                receptor.logFC = -0.1)

#Since the signaling genes in the `net.up` and `net.down` might be complex with multi-subunits, 
#we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)

gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


#We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble 
#plot or chord diagram. 
pairLR.use.up = net.up[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, 
                        targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, 
                        remove.isolate = T, 
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, 
                        targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

gg1 + gg2


#Visualizing the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
#Chord diagram
par(mfrow = c(1,2), xpd=TRUE)

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


#Visualizing the enriched ligands, signaling,or ligand-receptor pairs in one condition compared to
#another condition using wordcloud
#Visualizing the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human')

#Visualizing the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human')


#Visually comparing cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
#Similar to the CellChat analysis of individual dataset, we can visualize the cell-cell 
#communication network using Hierarchy plot, Circle plot or Chord diagram. 

#**Edge color/weight, node color/size/shape**: In all visualization plots, edge colors are 
#*consistent with the sources as sender, and edge weights are proportional to the interaction 
#*strength. Thicker edge line indicates a stronger signal. In the **Hierarchy plot and 
#*Circle plot**, circle sizes are proportional to the number of cells in each cell group. In the
#* hierarchy plot, solid and open circles represent source and target, respectively. In the
#* **Chord diagram**, the inner thinner bar colors represent the targets that receive signal from
#*  the corresponding outer bar. The inner bar size is proportional to the signal strength 
#*  received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. 
#*  Note that there exist some inner bars without any chord for some cell groups, please just 
#*  ignore it because this is an issue that has not been addressed by 
#*  [circlize](https://github.com/jokergoo/circlize) package. 
pathways.show <- c("CXCL") 

#Controlling the edge weights across different datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("CXCL") 

par(mfrow = c(1,2), xpd=TRUE)

ht <- list()

for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


#Chord diagram
pathways.show <- c("CXCL") 

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


#For the chord diagram, CellChat has an independent function `netVisual_chord_cell` to flexibly 
#visualize the signaling network by adjusting different parameters in the 
#[circlize](https://github.com/jokergoo/circlize) package. For example, we can define a named 
#char vector `group` to create multiple-group chord diagram, e.g., grouping cell clusters into 
#different cell types. 
#Chord diagram
#Grouping cell clusters into fibroblast, DC and TC cells
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) 

names(group.cellType) <- levels(object.list[[1]]@idents)

pathways.show <- c("CXCL") 

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", 
                                           names(object.list)[i]))
}

#Using chord diagram, CellChat provides two functions `netVisual_chord_cell` and 
#`netVisual_chord_gene` for visualizing cell-cell communication with different purposes and 
#different levels. `netVisual_chord_cell` is used for visualizing the cell-cell communication 
#between different cell groups (where each sector in the chord diagram is a cell group), and 
#`netVisual_chord_gene` is used for visualizing the cell-cell communication mediated by mutiple 
#ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, 
#receptor or signaling pathway.)
par(mfrow = c(1, 2), xpd=TRUE)
#Comparing all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, 
                       title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

#Comparing all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  
                       title.name = paste0("Signaling received by Inflam.DC and .TC - ", 
                                           names(object.list)[i]), legend.pos.x = 10)
}

#Showing all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),
                       slot.name = "netP", 
                       title.name = paste0("Signaling pathways sending from fibroblast - ", 
                                           names(object.list)[i]), legend.pos.x = 10)
}

#NB: Please ignore the note when generating the plot such as "Note: The first link end is drawn 
#out of sector 'MIF'.". If the gene names are overlapped, you can adjust the argument `small.gap`
#by decreasing the value. 



#Comparing the signaling gene expression distribution between different datasets 
#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling
#pathway using a Seurat wrapper function `plotGeneExpression`.  
#Setting factor level
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) 

plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)



#Saving the merged CellChat object
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_LS_vs_SS.rds")