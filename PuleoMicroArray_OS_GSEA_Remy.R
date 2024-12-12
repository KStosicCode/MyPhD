if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msigdb")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("hgu219.db")
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(tidyverse)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(hgu219.db)

devtools::install_github('GeNeHetX/CancerRNASig')
library(CancerRNASig)
CancerRNASig <- CancerRNASig::signatures$geneset

#That CancerRNASig object as input for fgsea, equivalent to "msigdb_ids"


library(fgsea)
library(data.table)
library(ggplot2)


#Adding gene names to the top table
top_table <- top.table2
top_table <- top_table[,1:7]
#ordb <- org.Hs.eg.db
#keytypes(org.Hs.eg.db)
#columns(org.Hs.eg.db)

hgu <- hgu219.db
keytypes(hgu219.db)

top_table$SYMBOL <- mapIds(hgu,
                           keys=top_table$Probe.Set.ID,
                           column="SYMBOL",
                           keytype="PROBEID",
                           multiVals="first")

#Removing NA values
length(top_table$SYMBOL)
index <- which(is.na(top_table$SYMBOL))
top_table <- top_table[-index,]
length(top_table$SYMBOL)

#Removing the duplicated names
index = which(duplicated(top_table$SYMBOL))
top_table <- top_table[-index,]
length((top_table$SYMBOL))

#Adding the corrected logFC and gene symbols
original_gene_list <- top_table$t
names(original_gene_list) <- top_table$SYMBOL
head(original_gene_list)


#Running fgsea for Pathways Analysis
fgseaRes <- fgsea(pathways = CancerRNASig, 
                  stats    = original_gene_list,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes)
head(fgseaRes[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(CancerRNASig[["PDAC_Bailey16_Immunogenic"]],
               original_gene_list) + labs(title="PDAC_Bailey16_Immunogenic")


#A table plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(CancerRNASig[topPathways], original_gene_list, fgseaRes,
              gseaParam = 0.5)

library(xlsx)
output1 <- fgseaRes
output1 <- apply(output1,2,as.character)
#Saving the fgsea results for the Pathways analysis
save(fgseaRes, file="OS_FGSEA_Remy_36vs18.Rdata")
write.csv(x=output1,"OS_FGSEA_Remy_36vs18.csv")
write.table(output1, file = "OS_FGSEA_Remy_36vs18.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_18_FGSEA\\OS_FGSEA_Remy_36vs18.xlsx")
#Importing the saved table
w <- read.delim("OS_FGSEA_Remy_36vs18.txt")