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

install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

#Package msigdb provides the gene sets in the MSigDB in the form of GeneSet
#objects. Accessory functions implemented in the GSEABase package provide a neat
#interface to interact with GeneSet objects. The ExperimentHub package processes
#the latest version of the MSigDB database into R objects that can be queried 
#using the GSEABase R/Bioconductor package. The entire database is stored in a 
#GeneSetCollection  object which in turn stores each signature as a GeneSet 
#object. Data in this package can be downloaded using the ExperimentHub interface.
#To download the data, what is needed first is to get a list of of the data
#available in the msigdb package. The query function assists in getting that list. 
eh = ExperimentHub()
query(eh , 'msigdb')

#Downloading the data using the custom accessor msigdb::getMsigdb()
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
msigdb.hs

#KEGG pathways (gene sets) cannot be included (integrated) within this
#ExperimentHub package due to licensing limitations. However, the data can be
#downloaded, processed and integrated directly from MSigDB by using appendKEGG() 
msigdb.hs = appendKEGG(msigdb.hs)
msigdb.hs

length(msigdb.hs)

listCollections(msigdb.hs)

listSubCollections(msigdb.hs)

#Gene Set Enrichment Analysis - Setting GeneSets for GO Analysis
subsetCollection(msigdb.hs, "c5")
#Retrieve collections
gene_set_collection_c5 <- subsetCollection(msigdb.hs, "c5")
gene_set_collection_c5 <- geneIds(gene_set_collection_c5)

#Gene Set Enrichment Analysis - Setting GeneSets for Pathways Analysis
subsetCollection(msigdb.hs, "c2", "CP")
#Retrieve collections
gene_set_collection_c2 <- subsetCollection(msigdb.hs, subcollection = "CP:REACTOME")
gene_set_collection_c2 <- geneIds(gene_set_collection_c2)


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
fgseaRes <- fgsea(pathways = gene_set_collection_c2, 
                  stats    = original_gene_list,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes)
head(fgseaRes[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(gene_set_collection_c2[["REACTOME_M_PHASE"]],
               original_gene_list) + labs(title="REACTOME_M_PHASE")


#A table plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gene_set_collection_c2[topPathways], original_gene_list, fgseaRes,
              gseaParam = 0.5)

library(xlsx)
output1 <- fgseaRes
output1 <- apply(output1,2,as.character)
#Saving the fgsea results for the Pathways analysis
save(fgseaRes, file="OS_36vs18_pathways.Rdata")
write.csv(x=output1,"OS_36vs18_pathways.csv")
write.table(output1, file = "OS_36vs18_pathways.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_18_FGSEA\\OS_36vs18_pathways.xlsx")
#Importing the saved table
w <- read.delim("OS_36vs18_pathways.txt")



#Running fgsea for GO Analysis
fgseaRes2 <- fgsea(pathways = gene_set_collection_c5, 
                   stats    = original_gene_list,
                   minSize  = 15,
                   maxSize  = 500)

head(fgseaRes2)
head(fgseaRes2[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(gene_set_collection_c5[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]],
               original_gene_list) + labs(title="GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS")


#A table plot
topPathwaysUp2 <- fgseaRes2[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown2 <- fgseaRes2[ES < 0][head(order(pval), n=10), pathway]
topPathways2 <- c(topPathwaysUp2, rev(topPathwaysDown2))
plotGseaTable(gene_set_collection_c5[topPathways2], original_gene_list, fgseaRes2,
              gseaParam = 0.5)


#collapsedPathways <- collapsePathways(fgseaRes2[order(pval)][padj < 0.01],
#gene_set_collection_c5, original_gene_list)

#mainPathways <- fgseaRes2[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

#plotGseaTable(gene_set_collection_c5[mainPathways], original_gene_list, fgseaRes2,
#gseaParam = 0.5)


output2 <- fgseaRes2
output2 <- apply(output2,2,as.character)
#Saving the fgsea results for the Gene Ontology analysis
save(fgseaRes2, file="OS_36vs18_geneontology.Rdata")
write.csv(x=output2,"OS_36vs18_geneontology.csv")
write.table(output2, file = "OS_36vs18_geneontology.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes2,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_18_FGSEA\\OS_36vs18_geneontology.xlsx")
fwrite(fgseaRes2, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))
#Importing the saved table
w <- read.delim("OS_36vs18_geneontology.txt")