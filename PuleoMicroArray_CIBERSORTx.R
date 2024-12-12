install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")
install.packages("remotes")
remotes::install_github("IOBR/IOBR")
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
library(preprocessCore)
library(parallel)
library(e1071)



expression_data <- Puleo_AllData$exp
metadata <- Puleo_AllData$samannot
gene_names <- Puleo_AllData$probeannot
survival <- Puleo_AllData$survdf
sig_matrix <- read.csv("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/CIBERSORTx/LM22.txt", header=TRUE, sep="\t")
set_cibersort_binary("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/CIBERSORTx/CIBERSORT.R")
set_cibersort_mat("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/CIBERSORTx/LM22_2.txt")



#Putting gene names as row names in signature matrix
head(rownames(sig_matrix))
rownames(sig_matrix) <- sig_matrix$Gene.symbol
head(rownames(sig_matrix))
sig_matrix <- sig_matrix[,-1]



#Saving the signature matrix with gene symbols as rownames
write.table(sig_matrix, file = "LM22_2.txt", row.names = F, sep = "\t", quote = F)



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
expression_data <- expression_data[,c(315,2:314,1)]
expression_data <- expression_data[,1:310]
expression_data <- expression_data[,2:310]



#Saving the counting matrix with gene symbols as rownames
save(expression_data, file="PuleoMicroArray_CountingMatrix_CIBERSORT.Rdata")
write.table(expression_data, file = "PuleoMicroArray_CountingMatrix_CIBERSORT.txt", row.names = F, sep = "\t", quote = F)
write.table(expression_data, file = "PuleoMicroArray_CountingMatrix_CIBERSORT_2.txt", row.names = F, sep = "\t", quote = F)



#Run CIBERSORT 
source("CIBERSORT.R")

res_cibersort <- CIBERSORT(
  sig_matrix,
  expression_data,
  perm=100,
  QN = TRUE,
  absolute=FALSE,
)


cibersort <- CIBERSORT(
  sig_matrix,
  expression_data,
  perm = 100,
  QN = TRUE,
  absolute = TRUE,
  abs_method = "sig.score"
)


res_ciber <- deconvolute(expression_data, "cibersort")



#CIBERSORTx output
res_cibersortx <- read.csv("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/CIBERSORTx/CIBERSORTx_Job6_Results.txt", header=TRUE, sep="\t")



#Merging CIBERSORTx output data frame and metadata
#Checking if CIBERSORTx output is a data frame
class(res_cibersortx)
#Merging with the metadata
m <- metadata
class(res_cibersortx)
class(metadata)

identical(rownames(res_cibersortx), rownames(metadata))
rownames(res_cibersortx) <- res_cibersortx$Mixture
identical(rownames(res_cibersortx), rownames(metadata))

metadata <- merge(res_cibersortx, metadata, by="row.names")
rownames(metadata) <- metadata$Row.names
class(metadata)



#Saving the metadata upgraded with CIBERSORTx output info
save(metadata, file="PuleoMicroArray_MetadataPlusCIBERSORTxOutput.Rdata")



#Filtering metadata according to Long-term survivors vs Short-term survivors
#RFS 36 vs 6 months 
#Creating a new column, which will eventually become a new RFS column with categorical values 
#over and under a threshold
length(metadata$Last.news.RFS.m.1)
metadata$RFS <- rep(10, times=309)

#Checking for the NAs and "xx" values and removing them
#NAs
which(is.na(metadata$Last.news.RFS.m.1))
length(which(is.na(metadata$Last.news.RFS.m.1)))

index <- which(is.na(metadata$Last.news.RFS.m.1))
#Even though you get the ordinal number of the indexed rows from the RFS column,
#when it comes to the assignment, it should be done on the whole metadata
metadata <- metadata[-index,]

which(is.na(metadata$Last.news.RFS.m.1))

#"xx"
which(metadata$Last.news.RFS.m.1=="xx")
length(which((metadata$Last.news.RFS.m.1=="xx")))

index <- which(metadata$Last.news.RFS.m.1=="xx")
#Even though you get the ordinal number of the indexed rows from the RFS column,
#when it comes to the assignment, it should be done on the whole metadata
metadata <- metadata[-index,]

which(metadata$Last.news.RFS.m.1=="xx")

class(metadata$Last.news.RFS.m.1)
#Converting the character into numeric values
as.numeric(metadata$Last.news.RFS.m.1)
#It's mandatory to assign the change to the variable of interest. Otherwise, it'll stay the same
metadata$Last.news.RFS.m.1 <- as.numeric(metadata$Last.news.RFS.m.1)
class(metadata$Last.news.RFS.m.1)



#Setting the threshold
index=which(metadata$Last.news.RFS.m.1 > 36)
metadata$RFS[index]="MoreThan36Months"

index=which(metadata$Last.news.RFS.m.1 < 6)
metadata$RFS[index]="LessThan6Months"

metadata$RFS
index=which(metadata$RFS=="10")
metadata <- metadata[-index,]
length(metadata$RFS)
metadata$RFS



#Keeping only those samples that have p-value below 0.05
index=which(metadata$P.value>0.05)
metadata <- metadata[-index,]

index=which(metadata$RFS=="MoreThan36Months")
index=which(metadata$RFS=="LessThan6Months")



#Plotting 
#B cells naive
ggplot(metadata,aes(RFS,B.cells.naive)) + 
  geom_boxplot() +
  ggtitle("B cells naive") +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_compare_means() 

bn <- ggplot(metadata,aes(RFS,B.cells.naive)) + 
  geom_boxplot() +
  ggtitle("B cells naive") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

bn

ggplot(metadata,aes(RFS,B.cells.naive)) + 
  geom_boxplot() +
  ggtitle("B cells naive") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

#B cells memory 
bm <- ggplot(metadata,aes(RFS,B.cells.memory)) + 
  geom_boxplot() +
  ggtitle("B cells memory") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

bm


#Plasma cells
pc <- ggplot(metadata,aes(RFS,Plasma.cells)) + 
  geom_boxplot() +
  ggtitle("Plasma cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

pc


#T cells that are CD8+
t_cd8 <- ggplot(metadata,aes(RFS, T.cells.CD8)) + 
  geom_boxplot() +
  ggtitle("T cells CD8+") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd8

#Naive B cells, Memory B cells, Plasma cells and T cells that are CD8+
ggarrange(bn, bm, pc, t_cd8)



#T cells that are CD4+ and naive
t_cd4_n <- ggplot(metadata,aes(factor(RFS), T.cells.CD4.naive)) + 
  geom_boxplot() +
  ggtitle("T cells CD4+ Naive") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd4_n


#T cells that are CD4+ and memory resting

ggplot(metadata,aes(factor(RFS), T.cells.CD4.memory.resting)) + 
  geom_boxplot() +
  ggtitle("T cells CD4+ Memory Resting") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd4_m_r <- ggplot(metadata,aes(factor(RFS), T.cells.CD4.memory.resting)) + 
  geom_boxplot() +
  ggtitle("T cells CD4+ Memory Resting") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd4_m_r


#T cells that are CD4+ and memory activated
t_cd4_m_a <- ggplot(metadata,aes(factor(RFS), T.cells.CD4.memory.activated)) + 
  geom_boxplot() +
  ggtitle("T cells CD4+ Memory Acivated") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_cd4_m_a


#T cells follicular helper
t_f_h <- ggplot(metadata,aes(factor(RFS),T.cells.follicular.helper)) + 
  geom_boxplot() +
  ggtitle("T Cells Follicular Helper") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)


t_f_h

#T CD4+ naive cells, T CD4+ memory resting cells , T CD4+ activated cells and T cells follicular helper

ggarrange(t_cd4_n, t_cd4_m_r, t_cd4_m_a, t_f_h)



#T cells regulatory - Treg
t_reg <- ggplot(metadata,aes(RFS,T.cells.regulatory..Tregs.)) + 
  geom_boxplot() +
  ggtitle("Regulatory T cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_reg


#T cells gamma delta
t_g_d <- ggplot(metadata,aes(factor(RFS),T.cells.gamma.delta)) + 
  geom_boxplot() +
  ggtitle("T cells Gamma Delta") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

t_g_d


#NK cells resting
nk_r <- ggplot(metadata,aes(factor(RFS),NK.cells.resting)) + 
  geom_boxplot() +
  ggtitle("NK Cells REsting") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

nk_r


#NK cells activated
nk_a <- ggplot(metadata,aes(factor(RFS),NK.cells.activated)) + 
  geom_boxplot() +
  ggtitle("NK cells Activated") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

nk_a

#Regulatory T cells, T cells gamma delta, NK cells resting and NK cells activated
ggarrange(t_reg, t_g_d, nk_r, nk_a)


#Monocytes
monocytes <- ggplot(metadata,aes(factor(RFS),Monocytes)) + 
  geom_boxplot() +
  ggtitle("Monocytes") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

monocytes


#Macrophages M0
m_m0 <- ggplot(metadata,aes(RFS, Macrophages.M0)) + 
  geom_boxplot() +
  ggtitle("Macrophages M0") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

m_m0


#Macrophages M1
m_m1 <- ggplot(metadata,aes(factor(RFS), Macrophages.M1)) + 
  geom_boxplot() +
  ggtitle("Macrophages M1") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

m_m1


#Macrophages M2
m_m2 <- ggplot(metadata,aes(factor(RFS), Macrophages.M2)) + 
  geom_boxplot() +
  ggtitle("Macrophages M2") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)


m_m2

#Monocytes, Macrophages M0, Macrophages M1 and Macrophages M2
ggarrange(monocytes, m_m0, m_m1, m_m2)


#Resting Dendritic cells
dendritic_r <- ggplot(metadata,aes(factor(RFS),Dendritic.cells.resting)) + 
  geom_boxplot() +
  ggtitle("Resting Dendritic Cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

dendritic_r


#Activated Dendritic cells
dendritic_a <- ggplot(metadata,aes(factor(RFS),Dendritic.cells.activated)) + 
  geom_boxplot() +
  ggtitle("Activated Dendritic cells") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

dendritic_a


#Mast cells resting
mast_r <- ggplot(metadata,aes(factor(RFS), Mast.cells.resting)) + 
  geom_boxplot() +
  ggtitle("Mast Cells Resting") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

mast_r


#Mast cells activated
mast_a <- ggplot(metadata,aes(factor(RFS), Mast.cells.activated)) + 
  geom_boxplot() +
  ggtitle("Mast Cells Activated") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

mast_a

#Resting Dendritic Cells, Activated Dendritic Cells, Mast Cells resting, Mast cells activated
ggarrange(dendritic_r, dendritic_a, mast_r, mast_a)


#Eosinophils
eosinophils <- ggplot(metadata,aes(factor(RFS), Eosinophils)) + 
  geom_boxplot() +
  ggtitle("Eosinophils") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

eosinophils


#Neutrophils
neutrophils <- ggplot(metadata,aes(factor(RFS), Neutrophils)) + 
  geom_boxplot() +
  ggtitle("Neutrophils") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_signif(comparisons = list(c("LessThan6Months", "MoreThan36Months")), 
              map_signif_level=TRUE)

neutrophils

#Eosinophils and Neutrophils
ggarrange(eosinophils, neutrophils)

#All
puleo_microarray_RFS_cibersortx_output <- ggarrange(t_cd8, t_reg, m_m0)
