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

Data@assayData$exprs
Data@annotation

library(readxl)
table <- read_xlsx("C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/SUBCLIN-24-3-20(1).xlsx")

#Creating new columns with only Y and N for each Puleo subtype
#Pure classical column
length(table$soustype)
table$PureClassical=rep(10, times=471)

index=which(table$soustype=="PurClassic") 
table$PureClassical[index]="Yes"

index=which(table$PureClassical==10)
table$PureClassical[index]="No"

#Immune classical column
table$ImmuneClassical=rep(10, times=471)

index=which(table$soustype=="ImmuClassic")
table$ImmuneClassical[index]="Yes"

index=which(table$ImmuneClassical==10)
table$ImmuneClassical[index]="No"

#Desmoplastic column
table$Desmoplastic=rep(10, times=471)

index=which(table$soustype=="Desmo")
table$Desmoplastic[index]="Yes"

index=which(table$Desmoplastic==10)
table$Desmoplastic[index]="No"

#Stroma activated column
table$StromaActivated=rep(10, times=471)

index=which(table$soustype=="Activ")
table$StromaActivated[index]="Yes"

index=which(table$StromaActivated==10)
table$StromaActivated[index]="No"

#Pure basal-like column
table$PureBasalLike=rep(10, times=471)

index=which(table$soustype=="PurBasal")
table$PureBasalLike[index]="Yes"

index=which(table$PureBasalLike==10)
table$PureBasalLike[index]="No"


#Modifying the Liver Metastasis column
index=which(table$foie==1)
table$foie[index]="LiverMetastasis"

index=which(table$foie==0 | table$foie==2)
table$foie[index]="NoLiverMetastasis"


#Modifying the Lung Metastasis column
index=which(table$poumon==1)
table$poumon[index]="LungMetastasis"

index=which(table$poumon==0 | table$poumon==2)
table$poumon[index]="NoLungMetastasis"

#Saving the modified version of the table
library("xlsx")
library("openxlsx")
library(readxl)
library(writexl)
write.csv(x=table,"PuleoTable.csv")
write_xlsx(as.data.frame(table),"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray\\PuleoTable.xlsx")



#Association of different metastasis with Puleo subtypes in primary tumour
#Liver metastasis
#Liver metastasis and Pure classical Puleo subtype
#Contingency table, Factorization and Relevelling  
table(table$foie, table$PureClassical)
table$foie <- factor(table$foie)
table$foie <- relevel(table$foie, ref="NoLiverMetastasis")

#Chi-squared test
chisq.test(table$foie, table$PureClassical)
#Fisher's test to get odds ratio (given that level "NoLiverMetastasis" was put 
#as a reference, if odds ratios is more than one, then LiverMetastasis is more
#likely to be associated with Pure classical Puleo subtype in a statistically
#significant way, whereas if it's smaller than one, then NoLiverMetastasis is 
#more likely to be associated with this Puleo subtype)
fisher.test(table$foie, table$PureClassical)

#Liver metastasis and Immune classical Puleo subtype
#Contingency table
table(table$foie, table$ImmuneClassical)

#Chi-squared test
chisq.test(table$foie, table$ImmuneClassical)

#Liver metastasis and Desmoplastic Puleo subtype
#Contingency table
table(table$foie, table$Desmoplastic)

#Chi-squared test
chisq.test(table$foie, table$Desmoplastic)

#Liver metastasis and Stroma activated Puleo subtype
#Contingency table
table(table$foie, table$StromaActivated)

#Chi-squared test
chisq.test(table$foie, table$StromaActivated)

#Liver metastasis and Pure basal-like Puleo subtype
#Contingency table
table(table$foie, table$PureBasalLike)

#Chi-squared test
chisq.test(table$foie, table$PureBasalLike)
#Fisher's test to get the odds ratio
fisher.test(table$foie, table$PureBasalLike)


#Lung metastasis
#Lung metastasis and Pure classical Puleo subtype
#Contingency table
table(table$poumon, table$PureClassical)
table$poumon <- factor(table$poumon)
table$poumon <- relevel(table$poumon, ref="NoLungMetastasis")

#Chi-squared test
chisq.test(table$poumon, table$PureClassical)

#Lung metastasis and Immune classical Puleo subtype
#Contingency table
table(table$poumon, table$ImmuneClassical)

#Fisher's test
fisher.test(table$poumon, table$ImmuneClassical)

#Lung metastasis and Desmoplastic Puleo subtype
#Contingency table
table(table$poumon, table$Desmoplastic)

#Chi-squared test
chisq.test(table$poumon, table$Desmoplastic)

#Lung metastasis and Stroma activated Puleo subtype
#Contingency table
table(table$poumon, table$StromaActivated)

#Chi-squared test
chisq.test(table$poumon, table$StromaActivated)

#Lung metastasis and Pure basal-like Puleo subtype
#Contingency table
table(table$poumon, table$PureBasalLike)

#Fisher test
fisher.test(table$poumon, table$PureBasalLike)



#Association of different metastasis with PurIST subtypes
#Liver
#Liver and Classical/Basal-like
#Contingency table
table(table$foie, table$purist)

#Chi-squared test
chisq.test(table$foie, table$purist)
#Fisher's test to get odds ratio
fisher.test(table$foie, table$purist)


#Lung
#Lung and classical/Basal-like
#Contingency table
table(table$poumon, table$purist)

#Chi-squared test
chisq.test(table$poumon, table$purist)