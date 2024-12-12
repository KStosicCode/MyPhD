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


#Checking if column of interest has numeric or character values
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

#Merging metadata table with the table containing Puleo subtypes and PurIST subtypes 
colnames(metadata)[7]
colnames(metadata)[7] <- "PatientID"
colnames(metadata)[7]

colnames(table)[1]
colnames(table)[1] <- "PatientID"
colnames(table)[1]

metadata <- merge(table, metadata, by="PatientID")

#Due to merging, some samples (those that were not present in both tables)
#disappeared, so their number should be checked after merging
index=which(metadata$OS=="MoreThan36Months")
index=which(metadata$OS=="LessThan18Months")



#Creating new columns with only Y and N for each Puleo subtype
#Pure classical column
length(metadata$soustype)
metadata$PureClassical=rep(10, times=201)

index=which(metadata$soustype=="PurClassic") 
metadata$PureClassical[index]="Yes"

index=which(metadata$PureClassical==10)
metadata$PureClassical[index]="No"

#Immune classical column
metadata$ImmuneClassical=rep(10, times=201)

index=which(metadata$soustype=="ImmuClassic")
metadata$ImmuneClassical[index]="Yes"

index=which(metadata$ImmuneClassical==10)
metadata$ImmuneClassical[index]="No"

#Desmoplastic column
metadata$Desmoplastic=rep(10, times=201)

index=which(metadata$soustype=="Desmo")
metadata$Desmoplastic[index]="Yes"

index=which(metadata$Desmoplastic==10)
metadata$Desmoplastic[index]="No"

#Stroma activated column
metadata$StromaActivated=rep(10, times=201)

index=which(metadata$soustype=="Activ")
metadata$StromaActivated[index]="Yes"

index=which(metadata$StromaActivated==10)
metadata$StromaActivated[index]="No"

#Pure basal-like column
metadata$PureBasalLike=rep(10, times=201)

index=which(metadata$soustype=="PurBasal")
metadata$PureBasalLike[index]="Yes"

index=which(metadata$PureBasalLike==10)
metadata$PureBasalLike[index]="No"


#Saving the modified version of the table
library("xlsx")
library("openxlsx")
library(readxl)
library(writexl)
write.csv(x=metadata,"OSTable.csv")
write_xlsx(as.data.frame(metadata),"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray\\OSTable.xlsx")



#Association of an overall survival with Puleo subtypes in primary tumour
#OS and Pure Classical Puleo subtype 
#Contingency table, Factorization and Relevelling 
table(metadata$OS, metadata$PureClassical)
metadata$OS <- factor(metadata$OS)
metadata$OS <- relevel(metadata$OS, ref="LessThan18Months")
#Chi-squared test
chisq.test(metadata$OS, metadata$PureClassical)
#Fisher's test to get odds ratio
fisher.test(metadata$OS, metadata$PureClassical)


#OS and Immune Classical Puleo subtype
#Contingency table
table(metadata$OS, metadata$ImmuneClassical)
#Chi-squared test
chisq.test(metadata$OS, metadata$ImmuneClassical)
#Fisher's test to get odds ratio
fisher.test(metadata$OS, metadata$ImmuneClassical)


#OS and Desmoplastic Puleo subtype
#Contingency table
table(metadata$OS, metadata$Desmoplastic)
#Chi-squared test
chisq.test(metadata$OS, metadata$Desmoplastic)
#Fishers's test to get odds ratio
fisher.test(metadata$OS, metadata$Desmoplastic)


#OS and Stroma Activated Puleo subtype
#Contingency table
table(metadata$OS, metadata$StromaActivated)
#Chi-squared test
chisq.test(metadata$OS, metadata$StromaActivated)
#Fisher's test to get odds ratio
fisher.test(metadata$OS, metadata$StromaActivated)


#OS and Pure basal-like Puleo subtype
#Contingency table
table(metadata$OS, metadata$PureBasalLike)
#Fisher's test because one category has less than 5 observations and also to get odds ratio
fisher.test(metadata$OS, metadata$PureBasalLike)



#Association of an overall survival with PurIST subtypes
#OS and Classical/Basal-like PurIST subtype
#Contingency table
table(metadata$OS, metadata$purist)
#Chi-squared test
chisq.test(metadata$OS, metadata$purist)
#Fisher's test to get odds ratio
fisher.test(metadata$OS, metadata$purist)
