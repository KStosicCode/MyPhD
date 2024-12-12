library(edgeR)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(oligo)
library(affy)

expression_data <- Puleo_AllData$exp
metadata <- Puleo_AllData$samannot
gene_names <- Puleo_AllData$probeannot
survival <- Puleo_AllData$survdf

library(readxl)
library("xlsx")
library("openxlsx")
library(readxl)
library(writexl)
table <- read_xlsx("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/SUBCLIN-24-3-20(1).xlsx")


#Table with counts
counts <- expression_data

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

index <- which(metadata$update.OS.time=="xx")
#Even though you get the ordinal number of the indexed rows from the OS column,
#when it comes to the assignment, it should be done on the whole metadata
metadata <- metadata[-index,]

which(metadata$update.OS.time=="xx")



#Checking if the column of interest has numeric or character values
class(metadata$update.OS.time)



#Setting the threshold
index=which(metadata$update.OS.time > 36)
metadata$OS[index]="MoreThan36Months"

index=which(metadata$update.OS.time < 12)
metadata$OS[index]="LessThan12Months"

metadata$OS
index=which(metadata$OS=="10")
metadata <- metadata[-index,]
length(metadata$OS)
metadata$OS


#DGE
#Performing RMA for Normalization, Variance stabilizing transformation and 
#Background noise correction
expression_data <- affy::rma(expression_data)
#Getting rid of those columns in the expression_data that were already removed
#from the rows in metadata (because columns in the expression_data correspond to
#the rows in metadata), by keeping only those columns in the expression_data
#that have column names corresponding to the row names of the metadata
rownames(metadata) <- metadata$idcit
expression_data <- expression_data[colnames(expression_data) %in% rownames(metadata)]
#Setting those two vectors in the same order, just in case that they got 
#accidentally reordered during the previous step
identical(colnames(expression_data), rownames(metadata))
expression_data <- expression_data[match(rownames(metadata), colnames(expression_data))]
colnames(expression_data)
rownames(metadata)
identical(colnames(expression_data), rownames(metadata))

#Create a design (model) matrix
group <- factor(metadata$OS)
design <- model.matrix(~0 + group)
design
colnames(design) <- levels(group)
design

#Creating a contrast matrix
cont.matrix <- makeContrasts(LongerOS_vs_ShorterOS=MoreThan36Months-LessThan12Months,
                             levels=design)
cont.matrix
dim(cont.matrix)

#Fitting the data into lm model
fit <- lmFit(expression_data,design)

#Creating a contrast matrix
fit.cont <- contrasts.fit(fit,cont.matrix)

#Empirical Bayes Statistics
fit.cont <- eBayes(fit.cont)
plotSA(fit.cont, main="Final model: Mean-variance trend")
dim(fit.cont)
summa.fit <- decideTests(fit.cont, p.value = 0.05, lfc=1)
summary(summa.fit)

#Extracting a general table of the top-ranked genes from a linear model fit and sorting them by p-values
top.table <- topTable(fit.cont, sort.by = "P", n = Inf) 
top.table<- topTable(fit.cont,coef=NULL,number=Inf,sort.by = "B",resort.by="logFC")
head(top.table, 5) 
#Getting rid of NAs in top.table
length(rownames(top.table))
index <- which(is.na(rownames(top.table)))
top.table <- top.table[-index,]
length(rownames(top.table))
head(top.table, 5) 
class(top.table)

#Merging this toptable with the gene_names one
top.table$Probe.Set.ID <- rownames(top.table)
length(top.table$Probe.Set.ID)
length(gene_names$Probe.Set.ID)
top.table2 <- merge(top.table, gene_names, by="Probe.Set.ID")
rownames(top.table2) <- top.table2$Gene.Symbol

#Saving the general table
save(top.table2, file="OS_36vs12_toptable.Rdata")
write.csv(x=top.table2,"OS_36vs12_toptable.csv")
write.table(top.table2, file = "OS_36vs12_toptable.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(top.table2,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_DGE\\OS_36vs12_toptable.xlsx")
#Importing the saved table
w <- read.delim("OS_36vs12_toptable.txt")

#Up and down regulated genes in LongerOS vs ShorterOS for OS 36 vs 12
library(dplyr)
LongerOS_vs_ShorterOS <- dplyr::filter(top.table2, (logFC>1|logFC< -1)&P.Value<0.05)
#Take only those that have no missing values
LongerOS_vs_ShorterOS <-  LongerOS_vs_ShorterOS[complete.cases(LongerOS_vs_ShorterOS), ]
#Up regulated
LongerOS_vs_ShorterOS_up <- LongerOS_vs_ShorterOS[LongerOS_vs_ShorterOS$logFC>0,]
#Down regulated
LongerOS_vs_ShorterOS_down <- LongerOS_vs_ShorterOS[LongerOS_vs_ShorterOS$logFC<0,]

#Saving the top table for Up regulated genes for OS 36 vs 12
save(LongerOS_vs_ShorterOS_up,file="OS_36vs12_lfc1_up.Rdata")
write.csv(x=LongerOS_vs_ShorterOS_up,"OS_36vs12_lfc1_up.csv")
write.table(LongerOS_vs_ShorterOS_up, file = "OS_36vs12_lfc1_up.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(LongerOS_vs_ShorterOS_up,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_DGE\\OS_36vs12_lfc1_up.xlsx")

#Saving the top table for Down regulated genes for OS 36 vs 12
save(LongerOS_vs_ShorterOS_down,file="OS_36vs12_lfc1_down.Rdata")
write.csv(x=LongerOS_vs_ShorterOS_down,"OS_36vs12_lfc1_down.csv")
write.table(LongerOS_vs_ShorterOS_down, file = "OS_36vs12_lfc1_down.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(LongerOS_vs_ShorterOS_down,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_DGE\\OS_36vs12_lfc1_down.xlsx")

#Volcano Plot
head((fit.cont$coefficients))
head(rownames(fit.cont$coefficients))
volcanoplot(fit.cont, coef = 1, highlight = 10, names = rownames(fit.cont$coefficients),
            main = "LongerOS_vs_ShorterOS")