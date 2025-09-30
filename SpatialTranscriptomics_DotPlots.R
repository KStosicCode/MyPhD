library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)


set.seed(1234)
#Loading the Seurat Object of interest
seurat_obj <- readRDS("LS_CAFs_SLITClassified.rds")


#Introducing signatures for different types of CAFs
adipoCAFs <- c("PTGDS", "LUM", "CFD", "FBLN1", "APOD", "DCN", "SFRP2", "MMP2")
c("PTGDS", "LUM", "CFD", "FBLN1", "APOD", "DCN", "SFRP2", "MMP2") %in% rownames(seurat_obj)

progenitorCAFs <- c("C7", "IGF1", "OGN", "PI16", "VEGFA")
c("C7", "IGF1", "OGN", "PI16", "VEGFA") %in% rownames(seurat_obj)

myCAFs <- c("LRRC15", "COL1A1", "CTHRC1", "POSTN", "ACTA2", "RGS5", "THBS1", "FN1", "TGFB1")
c("LRRC15", "COL1A1", "CTHRC1", "POSTN", "ACTA2", "RGS5", "THBS1", "FN1", "TGFB1") %in% rownames(seurat_obj)

engCAFs <- c("ENG")
"ENG" %in% rownames(seurat_obj)

adipoCAFs <- c("CCL2", "CXCL14", "CXCL12", "PDPN")
c("CCL2", "CXCL14", "CXCL12", "PDPN") %in% rownames(seurat_obj)


#SLIT classes that will be on the y-axis 
slit_levels <- c("SLIT2-_SLIT3-",
                 "SLIT2-_SLIT3+",
                 "SLIT2+_SLIT3-",
                 "SLIT2+_SLIT3+_low",
                 "SLIT2+_SLIT3+_high",
                 "Unclassified")    


#Defining the features (genes) to show on the x-axis of the first Dot Plot
features = c("PTGDS", "LUM", "CFD", "FBLN1", "APOD", "DCN", "SFRP2", "MMP2", 
             "C7", "IGF1", "OGN", "PI16", "VEGFA", "LRRC15", "COL1A1", "CTHRC1",
             "POSTN", "ACTA2", "RGS5", "THBS1", "FN1", "TGFB1", "ENG", "CCL2", 
             "CXCL14", "CXCL12", "PDPN")


#Setting identities 
#Making sure identities are set to SLIT classes
class(seurat_obj$SLIT_Class)

seurat_obj$SLIT_Class <- factor(seurat_obj$SLIT_Class, levels = slit_levels)

Idents(seurat_obj) <- seurat_obj$SLIT_Class
head(Idents(seurat_obj))


#Activated CAFs 
activatedCAFs_signature <- c("MMP11", "INHA", "COL10A1", "CTHRC1", "POSTN", "HOPX", "COL1A1",
                             "ACTA2", "COL8A1", "THY1", "COL11A1", "NTM", "COL5A1", "COL3A1",
                             "PDLIM3", "CST1", "PDGFC", "ITGBL1", "CDH11")
                             

c("MMP11", "INHA", "COL10A1", "CTHRC1", "POSTN", "HOPX", "COL1A1",
  "ACTA2", "COL8A1", "THY1", "COL11A1", "NTM", "COL5A1", "COL3A1",
  "PDLIM3", "CST1", "PDGFC", "ITGBL1", "CDH11") %in% rownames(seurat_obj)

#Basal CAFs
basalCAFs_signature <- c("GPX3", "CFD", "C7", "APOE", "LAMA2", "NAMPT", "APOD", "NR4A1", 
                         "STEAP4", "IGF1", "ADAMTS1", "CCL2", "SAT1", "SRPX", "NOVA1", "TGFBR3",
                         "S100A4", "SPARCL1", "LAMC1")
                     

c("GPX3", "CFD", "C7", "APOE", "LAMA2", "NAMPT", "APOD", "NR4A1", "STEAP4",
  "IGF1", "ADAMTS1", "CCL2", "SAT1", "SRPX", "NOVA1", "TGFBR3", "S100A4",
  "SPARCL1", "LAMC1") %in% rownames(seurat_obj)
                         
                         
#Dot plots - the size of the dot corresponds to the percentage of cells expressing the
#feature in each cluster. The color represents the average expression level
p <- DotPlot(seurat_obj, features = features, assay = "Spatial", dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.y = element_text(face = "bold"))

print(p)


q <- DotPlot(seurat_obj, features = activatedCAFs_signature, assay = "Spatial", dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.y = element_text(face = "bold"))

print(q)


r <- DotPlot(seurat_obj, features = basalCAFs_signature, assay = "Spatial", dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.y = element_text(face = "bold"))

print(r)