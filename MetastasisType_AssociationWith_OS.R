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

count_matrix <- Puleo_AllData$exp
metadata <- Puleo_AllData$samannot
gene_names <- Puleo_AllData$probeannot
survival <- Puleo_AllData$survdf

library(readxl)
table <- read_xlsx("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/SUBCLIN-24-3-20(1).xlsx")



#Modifying the Liver Metastasis column
index=which(table$foie==1)
table$foie[index]="LiverMetastasis"

index=which(table$foie==0 | table$foie==2)
table$foie[index]="NoLiverMetastasis"

index=which(metadata$foie==1)
metadata$foie[index]="LiverMetastasis"

index=which(metadata$foie==0 | metadata$foie==2)
metadata$foie[index]="NoLiverMetastasis"


#Modifying the Lung Metastasis column
index=which(table$poumon==1)
table$poumon[index]="LungMetastasis"

index=which(table$poumon==0 | table$poumon==2)
table$poumon[index]="NoLungMetastasis"

index=which(metadata$poumon==1)
metadata$poumon[index]="LungMetastasis"

index=which(metadata$poumon==0 | metadata$poumon==2)
metadata$poumon[index]="NoLungMetastasis"



#Saving the modified version of the table
library("xlsx")
library("openxlsx")
library(readxl)
library(writexl)
write.csv(x=table,"PuleoTable.csv")
write_xlsx(as.data.frame(table),"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray\\PuleoTable.xlsx")



#Merging metadata table with the table containing Puleo subtypes and PurIST subtypes 
colnames(metadata)[7]
colnames(metadata)[7] <- "PatientID"
colnames(metadata)[7]

colnames(table)[1]
colnames(table)[1] <- "PatientID"
colnames(table)[1]

metadata <- merge(table, metadata, by="PatientID")

length(which(metadata$foie.x=="LiverMetastasis"))
length(which(metadata$foie.x=="NoLiverMetastasis"))

length(which(metadata$poumon.x=="LungMetastasis"))
length(which(metadata$poumon.x=="NoLungMetastasis"))

#Modifying the OS column
index=which(metadata$update.OS.event==1)
metadata$update.OS.event[index]="WorseOS"

index=which(metadata$update.OS.event==0 | metadata$update.OS.event==2)
metadata$update.OS.event[index]="BetterOS"



#Association of different metastasis with OS
#Liver and OS
#Contingency table
table(metadata$foie.x, metadata$update.OS.event)
metadata$foie.x <- factor(metadata$foie.x)
metadata$foie.x <- relevel(metadata$foie.x, ref="NoLiverMetastasis")

#Chi-squared test
chisq.test(metadata$foie.x, metadata$update.OS.event)
#Fisher's test to get odds ratio
fisher.test(metadata$foie.x, metadata$update.OS.event)


#Lung and OS
#Contingency table
table(metadata$poumon.x, metadata$update.OS.event)
metadata$poumon.x <- factor(metadata$poumon.x)
metadata$poumon.x <- relevel(metadata$poumon.x, ref="NoLungMetastasis")

#Chi-squared test
chisq.test(metadata$poumon.x, metadata$update.OS.event)
#Fisher's test to get odds ratio
fisher.test(metadata$poumon.x, metadata$update.OS.event)
