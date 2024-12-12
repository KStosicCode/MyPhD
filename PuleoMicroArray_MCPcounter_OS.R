install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(hgu219.db)
library(edgeR)
library(tidyverse)
library(readr)
library(DGEobj.utils)
library(ggpubr)
library(ggsignif)



expression_data <- Puleo_AllData$exp
metadata <- Puleo_AllData$samannot
gene_names <- Puleo_AllData$probeannot
survival <- Puleo_AllData$survdf
ensembl_export <- read.csv("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/MCP_Counter/mart_export.txt")



#Adding the gene symbols column to the ensembl export
#ordb <- org.Hs.eg.db
#keytypes(org.Hs.eg.db)
#columns(org.Hs.eg.db)

hgu <- hgu219.db
keytypes(hgu219.db)

ensembl_export$SYMBOL <- mapIds(hgu,
                                keys=ensembl_export$Gene.stable.ID,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")

#Removing NA values
length(ensembl_export$SYMBOL)
index <- which(is.na(ensembl_export$SYMBOL))
ensembl_export <- ensembl_export[-index,]
length(ensembl_export$SYMBOL)
#Removing the duplicated names
index = which(duplicated(ensembl_export$SYMBOL))
ensembl_export <- ensembl_export[-index,]
length((ensembl_export$SYMBOL))



#Converting probe IDs in gene names table into gene symbols 
#ordb <- org.Hs.eg.db
#keytypes(org.Hs.eg.db)
#columns(org.Hs.eg.db)

hgu <- hgu219.db
keytypes(hgu219.db)

gene_names$SYMBOL <- mapIds(hgu,
                            keys=gene_names$Probe.Set.ID,
                            column="SYMBOL",
                            keytype="PROBEID",
                            multiVals="first")

#Removing NA values
length(gene_names$SYMBOL)
index <- which(is.na(gene_names$SYMBOL))
gene_names <- gene_names[-index,]
length(gene_names$SYMBOL)
#Removing the duplicated names
index = which(duplicated(gene_names$SYMBOL))
gene_names <- gene_names[-index,]
length((gene_names$SYMBOL))

head(rownames(gene_names))
head(rownames(expression_data))
tail(rownames(gene_names))
tail(rownames(expression_data))



#Merging counting matrix with gene names table in order to get gene symbols
expression_data$Probe.Set.ID <- rownames(expression_data)
length(expression_data$Probe.Set.ID)
length(gene_names$Probe.Set.ID)
expression_data <- merge(expression_data, gene_names, by="Probe.Set.ID")
#Removing the duplicated names
index = which(duplicated(expression_data$Gene.Symbol))
expression_data <- expression_data[-index,]
length((expression_data$Gene.Symbol))
length(unique(expression_data$Gene.Symbol))
#Removing NA values
length(expression_data$Gene.Symbol)
index <- which(is.na(expression_data$Gene.Symbol))
expression_data <- expression_data[-index,]
length(expression_data$Gene.Symbol)



#Naming rows of the counting matrix according to gene symbols
rownames(expression_data) <- expression_data$SYMBOL
#Getting rid of the new columns (keeping only columns that represent samples)
expression_data <- expression_data[,2:310]



#Saving the counting matrix with gene symbols as rownames
save(expression_data, file="PuleoMicroArray_CountingMatrix_MCPcounter.Rdata")



#Getting the gene lengths
gene_lengths_and_symbols <- ensembl_export[,c(5,6)]
#Keeping only those present in the counting matrix
gene_lengths_and_symbols <- gene_lengths_and_symbols[gene_lengths_and_symbols$SYMBOL %in% rownames(expression_data), ]
#Equalization of the number of rows in expression_data and gene_lengths_and_symbols
expression_data <- expression_data[rownames(expression_data) %in% gene_lengths_and_symbols$SYMBOL, ]
#Putting SYMBOL column from gene_lengths_and_symbols in the same order as rownames(expression_data)
identical(gene_lengths_and_symbols$SYMBOL, rownames(expression_data))
gene_lengths_and_symbols <- gene_lengths_and_symbols[match(rownames(expression_data), gene_lengths_and_symbols$SYMBOL), ]
identical(gene_lengths_and_symbols$SYMBOL, rownames(expression_data))

df <- gene_lengths_and_symbols
df2 <- expression_data
#Checking the class for expression_data, because the input for convertCounts() should be a matrix
class(expression_data)
expression_data <- as.matrix(expression_data)
class(expression_data)



#TPM normalisation
tpm <- convertCounts(
  countsMatrix=expression_data,
  unit="TPM",
  geneLength=gene_lengths_and_symbols$Transcript.length..including.UTRs.and.CDS.,
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)



#MCP-counter
counts <- tpm

deconvolution_methods
res_mcp_counter <- deconvolute(tpm, "mcp_counter")



#Merging MCP-counter output data frame and metadata
#Transposition of the MCP-counter output
#First remember the cell_type
n <- res_mcp_counter$cell_type
#Transpose all, but the first column (cell_type)
transposed <- as.data.frame(t(res_mcp_counter[,-1]))
colnames(transposed) <- n
#Merging with the metadata
m <- metadata
class(transposed)
class(metadata)
identical(rownames(transposed), rownames(metadata))
metadata <- merge(transposed, metadata, by="row.names")
rownames(metadata) <- metadata$Row.names
class(metadata)



#Saving the metadata upgraded with MCP-counter output info
save(metadata, file="PuleoMicroArray_MetadataPlusMCPcounterOutput.Rdata")



#Filtering metadata according to Long-term survivors vs Short-term survivors
#OS 36 vs 18 months 
#Creating a new column, which will eventually become a new OS column with categorical values 
#over and under a threshold
length(metadata$update.OS.time)
metadata$OS <- rep(10, times=309)

#Checking for the NAs and "xx" values and removing them
#NAs
which(is.na(metadata$update.OS.time))
length(which(is.na(metadata$update.OS.time)))

index <- which(is.na(metadata$update.OS.time))
#Even though you get the ordinal number of the indexed rows from the OS column,
#when it comes to the assignment, it should be done on the whole metadata
metadata <- metadata[-index,]

which(is.na(metadata$update.OS.time))

#"xx"
which(metadata$update.OS.time=="xx")
length(which((metadata$update.OS.time=="xx")))

index <- which(metadata$Last.news.RFS.m.1=="xx")
#Even though you get the ordinal number of the indexed rows from the RFS column,
#when it comes to the assignment, it should be done on the whole metadata
metadata <- metadata[-index,]

which(metadata$update.OS.time=="xx")

class(metadata$update.OS.time)
#Converting the character into numeric values
as.numeric(metadata$update.OS.time)
#It's mandatory to assign the change to the variable of interest. Otherwise, it'll stay the same
metadata$update.OS.time <- as.numeric(metadata$update.OS.time)
class(metadata$update.OS.time)



#Setting the threshold
index=which(metadata$update.OS.time > 36)
metadata$OS[index]="MoreThan36Months"

index=which(metadata$update.OS.time < 18)
metadata$OS[index]="LessThan18Months"

metadata$OS
index=which(metadata$OS=="10")
metadata <- metadata[-index,]
length(metadata$OS)
metadata$OS



#Plotting
colnames(metadata)[2] <- "T.cell"
colnames(metadata)[3] <- "T.cell.CD8"
colnames(metadata)[4] <- "Cytotoxicity.score"
colnames(metadata)[5] <- "NK.cell"
colnames(metadata)[6] <- "B.Cell"
colnames(metadata)[8] <- "Macrophage.Monocyte"
colnames(metadata)[9] <- "Myeloid.dendritic.cell"
colnames(metadata)[11] <- "Endothelial.cell"
colnames(metadata)[12] <- "Cancer.associated.fibroblast"



#T cells

t <- ggplot(metadata,aes(OS,T.cell)) + 
  geom_boxplot() +
  ggtitle("T cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

t <- ggplot(metadata,aes(OS,T.cell)) + 
  geom_boxplot() +
  ggtitle("T cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t <- ggplot(metadata,aes(factor(OS),T.cell)) + 
  geom_boxplot() +
  ggtitle("T cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

t


#T cells CD8+

t_cd8 <- ggplot(metadata,aes(OS,T.cell.CD8)) + 
  geom_boxplot() +
  ggtitle("T cells CD8+") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

t_cd8 <- ggplot(metadata,aes(OS,T.cell.CD8)) + 
  geom_boxplot() +
  ggtitle("T cells CD8+") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd8 <- ggplot(metadata,aes(factor(OS),T.cell.CD8)) + 
  geom_boxplot() +
  ggtitle("T cells CD8+") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

t_cd8


#Cytotoxicity score

cytotoxicity <- ggplot(metadata,aes(OS,Cytotoxicity.score)) + 
  geom_boxplot() +
  ggtitle("Cytotoxicity score") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

cytotoxicity <- ggplot(metadata,aes(OS,Cytotoxicity.score)) + 
  geom_boxplot() +
  ggtitle("Cytotoxicity score") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

cytotoxicity <- ggplot(metadata,aes(factor(OS),Cytotoxicity.score)) + 
  geom_boxplot() +
  ggtitle("Cytotoxicity score") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

cytotoxicity


#NK cells

nk <- ggplot(metadata,aes(OS, NK.cell)) + 
  geom_boxplot() +
  ggtitle("NK cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

nk <- ggplot(metadata,aes(OS, NK.cell)) + 
  geom_boxplot() +
  ggtitle("NK cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

nk <- ggplot(metadata,aes(factor(OS), NK.cell)) + 
  geom_boxplot() +
  ggtitle("NK cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

nk

#T cells, T cells that are CD8+, Cytotoxicity score and NK cells
ggarrange(t, t_cd8, cytotoxicity, nk)



#B cells

b <- ggplot(metadata,aes(OS, B.Cell)) + 
  geom_boxplot() +
  ggtitle("B Cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

b <- ggplot(metadata,aes(OS, B.Cell)) + 
  geom_boxplot() +
  ggtitle("B Cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

b <- ggplot(metadata,aes(factor(OS), B.Cell)) + 
  geom_boxplot() +
  ggtitle("B Cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

b


#Monocytes

monocytes <- ggplot(metadata,aes(OS,Monocyte)) + 
  geom_boxplot() +
  ggtitle("Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

monocytes <- ggplot(metadata,aes(OS,Monocyte)) + 
  geom_boxplot() +
  ggtitle("Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

monocytes <- ggplot(metadata,aes(factor(OS),Monocyte)) + 
  geom_boxplot() +
  ggtitle("Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

monocytes


#Macrophage/Monocyte

macrophages_monocytes <- ggplot(metadata,aes(OS,Macrophage.Monocyte)) + 
  geom_boxplot() +
  ggtitle("Macrophages/Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

macrophages_monocytes <- ggplot(metadata,aes(OS,Macrophage.Monocyte)) + 
  geom_boxplot() +
  ggtitle("Macrophages/Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

macrophages_monocytes <- ggplot(metadata,aes(factor(OS),Macrophage.Monocyte)) + 
  geom_boxplot() +
  ggtitle("Macrophages/Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

macrophages_monocytes


#Myeloid dendritic cell

myeloid_dendritic_cell <- ggplot(metadata,aes(OS,Myeloid.dendritic.cell)) + 
  geom_boxplot() +
  ggtitle("Myeloid.dendritic.cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

myeloid_dendritic_cell <- ggplot(metadata,aes(OS,Myeloid.dendritic.cell)) + 
  geom_boxplot() +
  ggtitle("Myeloid.dendritic.cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

myeloid_dendritic_cell <- ggplot(metadata,aes(factor(OS),Myeloid.dendritic.cell)) + 
  geom_boxplot() +
  ggtitle("Myeloid.dendritic.cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 


myeloid_dendritic_cell

#B cells, Monocytes, Macrophage/Monocyte and Myeloide dendritic cell

ggarrange(b, monocytes, macrophages_monocytes, myeloid_dendritic_cell)



#Neutrophils

neutrophils <- ggplot(metadata,aes(OS,Neutrophil)) + 
  geom_boxplot() +
  ggtitle("Neutrophil") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

neutrophils <- ggplot(metadata,aes(OS,Neutrophil)) + 
  geom_boxplot() +
  ggtitle("Neutrophil") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

neutrophils <- ggplot(metadata,aes(factor(OS),Neutrophil)) + 
  geom_boxplot() +
  ggtitle("Neutrophil") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

neutrophils


#Endothelial cell

endothelial_cells <- ggplot(metadata,aes(OS,Endothelial.cell)) + 
  geom_boxplot() +
  ggtitle("Endothelial cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

endothelial_cells <- ggplot(metadata,aes(OS,Endothelial.cell)) + 
  geom_boxplot() +
  ggtitle("Endothelial cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

endothelial_cells <- ggplot(metadata,aes(factor(OS),Endothelial.cell)) + 
  geom_boxplot() +
  ggtitle("Endothelial cell") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

endothelial_cells


#Cancer associated fibroblast 

caf <- ggplot(metadata,aes(OS,Cancer.associated.fibroblast)) + 
  geom_boxplot() +
  ggtitle("Cancer associated fibroblast") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means()

caf <- ggplot(metadata,aes(OS,Cancer.associated.fibroblast)) + 
  geom_boxplot() +
  ggtitle("Cancer associated fibroblast") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan18Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

caf <- ggplot(metadata,aes(factor(OS),Cancer.associated.fibroblast)) + 
  geom_boxplot() +
  ggtitle("Cancer associated fibroblast") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

caf

#Neutrophils, Endothelial cells and Cancer associated fibroblasts
ggarrange(neutrophils, endothelial_cells, caf)



#All
puleo_microarray_outputOS1 <- ggarrange(t, t_cd8, cytotoxicity, nk)
puleo_microarray_outputOS2 <- ggarrange(b, monocytes, macrophages_monocytes, myeloid_dendritic_cell)
puleo_microarray_outputOS3 <- ggarrange(neutrophils, endothelial_cells, caf)
