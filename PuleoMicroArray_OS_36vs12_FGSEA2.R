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

devtools::install_github('GeNeHetX/CancerRNASig')
library(CancerRNASig)
CancerRNASig <- CancerRNASig::signatures$geneset

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

#Gene Set Enrichment Analysis - Setting GeneSets for TFT Analysis
subsetCollection(msigdb.hs, "c3", "TFT")
#Retrieve the "C3" collection for transcription factor targets
gene_set_collection_c3 <- subsetCollection(msigdb.hs, "c3", subcollection = "TFT:GTRD")
gene_set_collection_c3 <- geneIds(gene_set_collection_c3)

#Gene Set Enrichment Analysis - Setting GeneSets for Immunologic Analysis
subsetCollection(msigdb.hs, "c7")
#Retrieve the "C7" collection for immunologic gene sets
gene_set_collection_c7 <- subsetCollection(msigdb.hs, "c7")
gene_set_collection_c7 <- geneIds(gene_set_collection_c7)



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


#Defining a function to color the pathways based on their enrichment score
qGseaTable_erraz=function (pathways, stats, fgseaRes, gseaParam = 1,
                           colwidths = c(5,3, 0.8, 1.2, 1.2),rename=NULL) {
  require(fgsea)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  n <- length(statsAdj)
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  if(!is.null(rename)&length(rename)==length(pathways)){
    names(rename)=names(pathways)
    
  }else{
    rename=names(pathways)
    names(rename)=names(pathways)
    
  }
  
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(rename[pn], just = "right", x = unit(0.95, "npc"),gp=gpar(col=if (sprintf("%.2f", annotation$NES)<0) {"blue"} else {"red"})),
         ggplot() + geom_segment(aes(x = p, xend = p, y = 0,
                                     yend = statsAdj[p],col=p), size = 0.2) + scale_x_continuous(limits = c(0,
                                                                                                            length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1,
                                                                                                                                                                                 1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
           scale_colour_gradient2(midpoint=n/2,low = '#ca0020', mid="#f7f7f7",  high = '#0571b0') +
           theme(panel.background = element_blank(), legend.position="none", axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0,
                                                                                 "null"), 4)), textGrob(sprintf("%.2f", annotation$NES),gp=gpar(col=if (sprintf("%.2f", annotation$NES)<0) {"blue"} else {"red"})),
         textGrob(sprintf("%.1e", annotation$pval),gp=gpar(col=if (sprintf("%.2f", annotation$NES)<0) {"blue"} else {"red"})), textGrob(sprintf("%.1e",annotation$padj),gp=gpar(col=if (sprintf("%.2f", annotation$NES)<0) {"blue"} else {"red"})))
  })
  rankPlot <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0,
                                                                      length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1,
                                                                                                                                           1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) + theme(panel.background = element_blank(),
                                                                                                                                                                                                   axis.line = element_blank(), axis.text.y = element_blank(),
                                                                                                                                                                                                   axis.ticks.y = element_blank(), panel.grid = element_blank(),
                                                                                                                                                                                                   axis.title = element_blank(), plot.margin = unit(c(0,
                                                                                                                                                                                                                                                      0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0,
                                                                                                                                                                                                                                                                                                 0, 0), "npc"))
  grobs <- c(lapply(c("Signature", "Gene ranks", "NES", "pval",
                      "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(),
                                                                              rankPlot, nullGrob(), nullGrob(), nullGrob()))
  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))
  gridExtra::grid.arrange(grobs = grobs[grobsToDraw], ncol = sum(colwidths !=
                                                                   0), widths = colwidths[colwidths != 0])
  
  
}



#An enrichment plot for one of the most prominent output pathways
plotEnrichment(gene_set_collection_c2[["REACTOME_M_PHASE"]],
               original_gene_list) + labs(title="REACTOME_M_PHASE")


#A table plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=17), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gene_set_collection_c2[topPathways], original_gene_list, fgseaRes,
              gseaParam = 0.5)

#Coloring the pathways at the table plot according to their enrichment score
library(fgsea)
library(gridExtra)
library(ggplot2)
library(grid)

qGseaTable_erraz(
  pathways = gene_set_collection_c2[topPathways],
  stats = original_gene_list,
  fgseaRes = fgseaRes
)

#Specify pathways of interest
pathways_of_interest <- c("REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES", 
                          "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS", 
                          "REACTOME_PD_1_SIGNALING",
                          "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
                          "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
                          "REACTOME_INTERFERON_GAMMA_SIGNALING",
                          "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY",
                          "REACTOME_SELENOAMINO_ACID_METABOLISM",
                          "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
                          "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS")

#Filter the FGSEA results
fgseaRes_filtered <- fgseaRes[fgseaRes$pathway %in% pathways_of_interest, ]

#Split the filtered FGSEA results into upregulated and downregulated pathways
fgseaRes_up <- fgseaRes_filtered[fgseaRes_filtered$NES > 0, ]
fgseaRes_down <- fgseaRes_filtered[fgseaRes_filtered$NES < 0, ]

#Sort each group by adjusted p-value (padj)
fgseaRes_up <- fgseaRes_up[order(fgseaRes_up$padj), ]
fgseaRes_down <- fgseaRes_down[order(fgseaRes_down$padj, decreasing = TRUE), ]

#Combine the results with upregulated pathways first, followed by downregulated pathways
fgseaRes_sorted <- rbind(fgseaRes_up, fgseaRes_down)

#Use the sorted pathways for the GSEA table plot
qGseaTable_erraz(
  pathways = gene_set_collection_c2[fgseaRes_sorted$pathway],
  stats = original_gene_list,
  fgseaRes = fgseaRes_sorted
)

library(xlsx)
library(writexl)
output1 <- fgseaRes_sorted
output1 <- apply(output1,2,as.character)
#Saving the fgsea results for the Pathways analysis
save(fgseaRes_sorted, file="OS_36vs12_pathways.Rdata")
write.csv(x=output1,"OS_36vs12_pathways.csv")
write.table(output1, file = "OS_36vs12_pathways.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_FGSEA_2\\OS_36vs12_pathways.xlsx")
#Importing the saved table
w <- read.delim("OS_36vs12_pathways.txt")



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

#Coloring the gene ontology at the table plot according to their enrichment score

qGseaTable_erraz(
  pathways = gene_set_collection_c5[topPathways2],
  stats = original_gene_list,
  fgseaRes = fgseaRes2
)


#Specify gene ontologies of interest
gene_ontologies_of_interest <- c("GOBP_ADAPTIVE_IMMUNE_RESPONSE", 
                                 "GOBP_T_CELL_ACTIVATION", 
                                 "GOMF_IMMUNE_RECEPTOR_ACTIVITY",
                                 "GOBP_LYMPHOCYTE_DIFFERENTIATION",
                                 "GOBP_T_CELL_DIFFERENTIATION",
                                 "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION",
                                 "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_ACTIVATION",
                                 "GOBP_ALPHA-BETA_T_CELL_ACTIVATION",
                                 "GOBP_REGULATION_OF_LEUKOCYTE_DIFFERENTIATION",
                                 "GOBP_TRANSLATIONAL_INITIATION",
                                 "GOBP_PROTEIN_TARGETING_TO_ER",
                                 "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
                                 "GOMF_CADHERIN_BINDING",
                                 "GOCC_FOCAL_ADHESION",
                                 "GOCC_CELL-SUBSTRATE_JUNCTION",
                                 "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
                                 "GOCC_RIBOSOMAL_SUBUNIT",
                                 "GOCC_RIBOSOME")

#Filter the FGSEA results
fgseaRes2_filtered <- fgseaRes2[fgseaRes2$pathway %in% gene_ontologies_of_interest, ]

#Split the filtered FGSEA results into upregulated and downregulated pathways
fgseaRes2_up <- fgseaRes2_filtered[fgseaRes2_filtered$NES > 0, ]
fgseaRes2_down <- fgseaRes2_filtered[fgseaRes2_filtered$NES < 0, ]

#Sort each group by adjusted p-value (padj)
fgseaRes2_up <- fgseaRes2_up[order(fgseaRes2_up$padj), ]
fgseaRes2_down <- fgseaRes2_down[order(fgseaRes2_down$padj, decreasing = TRUE), ]

#Combine the results with upregulated pathways first, followed by downregulated pathways
fgseaRes2_sorted <- rbind(fgseaRes2_up, fgseaRes2_down)

#Use the sorted pathways for the GSEA table plot
qGseaTable_erraz(
  pathways = gene_set_collection_c5[fgseaRes2_sorted$pathway],
  stats = original_gene_list,
  fgseaRes = fgseaRes2_sorted
)


output2 <- fgseaRes2_sorted
output2 <- apply(output2,2,as.character)
#Saving the fgsea results for the Gene Ontology analysis
save(fgseaRes2_sorted, file="OS_36vs12_geneontology.Rdata")
write.csv(x=output2,"OS_36vs12_geneontology.csv")
write.table(output2, file = "OS_36vs12_geneontology.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes2,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_FGSEA_2\\OS_36vs12_geneontology.xlsx")
fwrite(fgseaRes2_sorted, file="OS_36vs12_geneontology.txt", sep="\t", sep2=c("", " ", ""))
#Importing the saved table
w <- read.delim("OS_36vs12_geneontology.txt")



#That CancerRNASig object as input for fgsea, equivalent to "msigdb_ids"
#Running fgsea for CancerRNASig Pathways Analysis
fgseaRes3 <- fgsea(pathways = CancerRNASig, 
                   stats    = original_gene_list,
                   minSize  = 15,
                   maxSize  = 500)

head(fgseaRes3)
head(fgseaRes3[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(CancerRNASig[["PDAC_Bailey16_Immunogenic"]],
               original_gene_list) + labs(title="PDAC_Bailey16_Immunogenic")

renamed_pathways <- c(
  "PDAC_CSY20_ClassicalA" = "Chan-Seng-Yue et al. Classical A",
  "PDAC_Hwang22_Malignlineage.Classicallike" = "Hwang et al. Lineage Classical like",
  "PDAC_Bailey16_Immunogenic" = "Bailey et al. Immunogenic",
  "CAF_Turley20_hCAF1early" = "Dominguez et al. Early CAF1 (iCAFs)",
  "PDAC_CSY20_ClassicalB" = "Chan-Seng-Yue et al. Classical B",
  "PDAC_Hwang22_Malignlineage.Neurallikeprogenitor" = "Hwang et al. Lineage Neural like progenitor",
  "PDAC_Bailey16_Progenitor" = "Bailey et al. Progenitor",
  "CAF_FMG20_Normal.Fibroblast" = "Kieffer et al. Normal Fibroblasts",
  "CAF_FMG20_detox.iCAF" = "Kieffer et al. Detoxification Pathway iCAFs",
  "CAF_FMG20_IFNγ.iCAF" = "Kieffer et al. IFNγ iCAFs",
  "PDAC_Hwang22_Malignstate.Interferonsignaling" = "Hwang et al. State IFN signaling",
  "CAF_Neuzillet22_PDPN" = "Neuzillet et al. CAFs expressing immune-related pathways",
  "PDAC_Hwang22_Malignstate.Adhesive" = "Hwang et al. State Adhesive",
  "PDAC_Hwang22_Malignstate.CyclingG2M" = "Hwang et al. State Cycling G2/M",
  "ECM_Helms22_PSCcaf" = "Helms et al. PSC-derived CAFs",
  "CAF_Turley20_hCAF0TGFB" = "Dominguez et al. TGFβ CAFs (myCAFs)",
  "PDAC_Hwang22_Malignstate.TNFNFkBsignaling" = "Hwang et al. State TNF NFkB signaling",
  "PDAC_CSY20_BasallikeA" = "Chan-Seng-Yue et al. Basal like A",
  "PDAC_Hwang22_Malignlineage.Squamoid" = "Hwang et al. Lineage Squamoid",
  "PDAC_Hwang22_Malignstate.Ribosomal" = "Hwang et al. State Ribosomal",
  "PDAC_Hwang22_Malignlineage.Basaloid" = "Hwang et al. Lineage Basaloid"
)

# Rename pathways in fgseaRes3
fgseaRes3$pathway <- ifelse(fgseaRes3$pathway %in% names(renamed_pathways),
                            renamed_pathways[fgseaRes3$pathway],
                            fgseaRes3$pathway)

# Rename pathways in CancerRNASig
names(CancerRNASig) <- ifelse(names(CancerRNASig) %in% names(renamed_pathways),
                              renamed_pathways[names(CancerRNASig)],
                              names(CancerRNASig))

#A table plot
topPathwaysUp3 <- fgseaRes3[ES > 0][head(order(pval), n=25), pathway]
topPathwaysDown3 <- fgseaRes3[ES < 0][head(order(pval), n=25), pathway]
topPathways3 <- c(topPathwaysUp3, rev(topPathwaysDown3))
plotGseaTable(CancerRNASig[topPathways3], original_gene_list, fgseaRes3,
              gseaParam = 0.5)

#Coloring the CancerRNASig at the table plot according to their enrichment score

qGseaTable_erraz(
  pathways = CancerRNASig[topPathways3],
  stats = original_gene_list,
  fgseaRes = fgseaRes3
)

#Specify pathways of interest (taking the newly named signatures from renamed_pathways)
pathways_of_interest3 <- c("Chan-Seng-Yue et al. Classical A", 
                          "Hwang et al. Lineage Classical like", 
                          "Bailey et al. Immunogenic",
                          "Dominguez et al. Early CAF1 (iCAFs)",
                          "Chan-Seng-Yue et al. Classical B",
                          "Hwang et al. Lineage Neural like progenitor",
                          "Bailey et al. Progenitor",
                          "Kieffer et al. Normal Fibroblasts",
                          "Kieffer et al. Detoxification Pathway iCAFs",
                          "Kieffer et al. IFNγ iCAFs",
                          "Hwang et al. State IFN signaling",
                          "Neuzillet et al. CAFs expressing immune-related pathways",
                          "Hwang et al. State Adhesive",
                          "Hwang et al. State Cycling G2/M",
                          "Helms et al. PSC-derived CAFs",
                          "Dominguez et al. TGFβ CAFs (myCAFs)",
                          "Hwang et al. State TNF NFkB signaling",
                          "Chan-Seng-Yue et al. Basal like A",
                          "Hwang et al. Lineage Squamoid",
                          "Hwang et al. State Ribosomal",
                          "Hwang et al. Lineage Basaloid")

#Filter the FGSEA results
fgseaRes3_filtered <- fgseaRes3[fgseaRes3$pathway %in% pathways_of_interest3, ]

#Split the filtered FGSEA results into upregulated and downregulated pathways
fgseaRes3_up <- fgseaRes3_filtered[fgseaRes3_filtered$NES > 0, ]
fgseaRes3_down <- fgseaRes3_filtered[fgseaRes3_filtered$NES < 0, ]

#Sort each group by adjusted p-value (padj)
fgseaRes3_up <- fgseaRes3_up[order(fgseaRes3_up$padj), ]
fgseaRes3_down <- fgseaRes3_down[order(fgseaRes3_down$padj, decreasing = TRUE), ]

#Combine the results with upregulated pathways first, followed by downregulated pathways
fgseaRes3_sorted <- rbind(fgseaRes3_up, fgseaRes3_down)

#Use the sorted pathways for the GSEA table plot
qGseaTable_erraz(
  pathways = CancerRNASig[fgseaRes3_sorted$pathway],
  stats = original_gene_list,
  fgseaRes = fgseaRes3_sorted
)

#The same, just without CAF signatures
#Specify pathways of interest (taking the newly named signatures from renamed_pathways)
pathways_of_interest32 <- c("Chan-Seng-Yue et al. Classical A", 
                           "Hwang et al. Lineage Classical like", 
                           "Bailey et al. Immunogenic",
                           "Chan-Seng-Yue et al. Classical B",
                           "Hwang et al. Lineage Neural like progenitor",
                           "Bailey et al. Progenitor",
                           "Hwang et al. State IFN signaling",
                           "Hwang et al. State Adhesive",
                           "Hwang et al. State Cycling G2/M",
                           "Hwang et al. State TNF NFkB signaling",
                           "Chan-Seng-Yue et al. Basal like A",
                           "Hwang et al. Lineage Squamoid",
                           "Hwang et al. State Ribosomal",
                           "Hwang et al. Lineage Basaloid")

#Filter the FGSEA results
fgseaRes3_filtered <- fgseaRes3[fgseaRes3$pathway %in% pathways_of_interest32, ]

#Split the filtered FGSEA results into upregulated and downregulated pathways
fgseaRes3_up <- fgseaRes3_filtered[fgseaRes3_filtered$NES > 0, ]
fgseaRes3_down <- fgseaRes3_filtered[fgseaRes3_filtered$NES < 0, ]

#Sort each group by adjusted p-value (padj)
fgseaRes3_up <- fgseaRes3_up[order(fgseaRes3_up$padj), ]
fgseaRes3_down <- fgseaRes3_down[order(fgseaRes3_down$padj, decreasing = TRUE), ]

#Combine the results with upregulated pathways first, followed by downregulated pathways
fgseaRes3_sorted <- rbind(fgseaRes3_up, fgseaRes3_down)

#Use the sorted pathways for the GSEA table plot
qGseaTable_erraz(
  pathways = CancerRNASig[fgseaRes3_sorted$pathway],
  stats = original_gene_list,
  fgseaRes = fgseaRes3_sorted
)

#The same, but only CAFs signatures
#Specify pathways of interest (taking the newly named signatures from renamed_pathways)
pathways_of_interest33 <- c("Dominguez et al. Early CAF1 (iCAFs)",
                           "Kieffer et al. Normal Fibroblasts",
                           "Kieffer et al. Detoxification Pathway iCAFs",
                           "Kieffer et al. IFNγ iCAFs",
                           "Neuzillet et al. CAFs expressing immune-related pathways",
                           "Helms et al. PSC-derived CAFs",
                           "Dominguez et al. TGFβ CAFs (myCAFs)")

#Filter the FGSEA results
fgseaRes3_filtered <- fgseaRes3[fgseaRes3$pathway %in% pathways_of_interest33, ]

#Split the filtered FGSEA results into upregulated and downregulated pathways
fgseaRes3_up <- fgseaRes3_filtered[fgseaRes3_filtered$NES > 0, ]
fgseaRes3_down <- fgseaRes3_filtered[fgseaRes3_filtered$NES < 0, ]

#Sort each group by adjusted p-value (padj)
fgseaRes3_up <- fgseaRes3_up[order(fgseaRes3_up$padj), ]
fgseaRes3_down <- fgseaRes3_down[order(fgseaRes3_down$padj, decreasing = TRUE), ]

#Combine the results with upregulated pathways first, followed by downregulated pathways
fgseaRes3_sorted <- rbind(fgseaRes3_up, fgseaRes3_down)

#Use the sorted pathways for the GSEA table plot
qGseaTable_erraz(
  pathways = CancerRNASig[fgseaRes3_sorted$pathway],
  stats = original_gene_list,
  fgseaRes = fgseaRes3_sorted
)

library(xlsx)
output3 <- fgseaRes3_sorted
output3 <- apply(output3,2,as.character)
#Saving the fgsea results for the Pathways analysis
save(fgseaRes3_sorted, file="OS_36vs12_Remy.Rdata")
write.csv(x=output3,"OS_36vs12_Remy.csv")
write.table(output3, file = "OS_36vs12_Remy.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes3_sorted,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_FGSEA_2\\OS_36vs12_Remy.xlsx")
#Importing the saved table
w <- read.delim("OS_36vs12_Remy.txt")



#Running fgsea for TFT:GTRD Analysis
fgseaRes4 <- fgsea(pathways = gene_set_collection_c3, 
                   stats    = original_gene_list,
                   minSize  = 15,
                   maxSize  = 500)

head(fgseaRes4)
head(fgseaRes4[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(gene_set_collection_c3[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]],
               original_gene_list) + labs(title="GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS")


#A table plot
topPathwaysUp4 <- fgseaRes4[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown4 <- fgseaRes4[ES < 0][head(order(pval), n=10), pathway]
topPathways4 <- c(topPathwaysUp4, rev(topPathwaysDown4))
plotGseaTable(gene_set_collection_c3[topPathways4], original_gene_list, fgseaRes4,
              gseaParam = 0.5)

#Coloring the gene ontology at the table plot according to their enrichment score

qGseaTable_erraz(
  pathways = gene_set_collection_c3[topPathways4],
  stats = original_gene_list,
  fgseaRes = fgseaRes4
)


output4 <- fgseaRes4
output4 <- apply(output4,2,as.character)
#Saving the fgsea results for the Transcription Factor Target analysis
save(fgseaRes4, file="OS_36vs12_transcriptionfactortarget.Rdata")
write.csv(x=output4,"OS_36vs12_transcriptionfactortarget.csv")
write.table(output4, file = "OS_36vs12_transcriptionfactortarget.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes4,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_FGSEA_2\\OS_36vs12_transcriptionfactortarget.xlsx")
fwrite(fgseaRes4, file="OS_36vs12_transcriptionfactortarget.txt", sep="\t", sep2=c("", " ", ""))
#Importing the saved table
w <- read.delim("OS_36vs12_transcriptionfactortarget.txt")



#Running fgsea for Immunologic Gene Set Analysis
fgseaRes5 <- fgsea(pathways = gene_set_collection_c7, 
                   stats    = original_gene_list,
                   minSize  = 15,
                   maxSize  = 500)

head(fgseaRes5)
head(fgseaRes5[order(pval), ])


#An enrichment plot for one of the most prominent output pathways
plotEnrichment(gene_set_collection_c3_2[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]],
               original_gene_list) + labs(title="GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS")


#A table plot
topPathwaysUp5 <- fgseaRes5[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown5 <- fgseaRes5[ES < 0][head(order(pval), n=20), pathway]
topPathways5 <- c(topPathwaysUp5, rev(topPathwaysDown5))
plotGseaTable(gene_set_collection_c7[topPathways5], original_gene_list, fgseaRes5,
              gseaParam = 0.5)

#Coloring the gene ontology at the table plot according to their enrichment score

qGseaTable_erraz(
  pathways = gene_set_collection_c7[topPathways5],
  stats = original_gene_list,
  fgseaRes = fgseaRes5
)


output5 <- fgseaRes5
output5 <- apply(output5,2,as.character)
#Saving the fgsea results for the Immunologic analysis
save(fgseaRes5, file="OS_36vs12_immunologic.Rdata")
write.csv(x=output5,"OS_36vs12_immunologic.csv")
write.table(output5, file = "OS_36vs12_immunologic.txt", row.names = F, sep = "\t", quote = F)
write_xlsx(fgseaRes5,"C:/Users/Admin/Desktop/R Studio for Bulk and Single Cell Analysis/Puleo Microarray/OS_36_vs_12_FGSEA_2\\OS_36vs12_immunologic.xlsx")
fwrite(fgseaRes5, file="OS_36vs12_immunologic.txt", sep="\t", sep2=c("", " ", ""))
#Importing the saved table
w <- read.delim("OS_36vs12_immunologic.txt")