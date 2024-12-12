library(edgeR)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(openxlsx)
library(dplyr)
library(stringr)
library(cowplot)
library("ggpubr")

setwd("/Users/Admin/Desktop/PhD_2023/ChristelleData/")

#Count Matrix Christelle ####
#Load Count Table
counts<- read.table("LCM_count_table.txt")
pattern <- "X"
str_detect(colnames(counts), pattern)
colnames(counts) <- str_replace(colnames(counts),pattern,"")

samples_keep <- c("210","58","67","211","59","213","102")

samplesdata <-read.csv2("LCM_samplesheet.csv",header=T)
samplesdata
index <- which(samplesdata$sample%in%samples_keep)
samplesdata <- samplesdata[index,]

#Introduce samples:
index <- which(colnames(counts)%in%samples_keep)
counts <- counts[,index]


final_matrix <- DGEList(counts,samples=colnames(counts),group=samplesdata$group)


#Filtering
keep <- filterByExpr(final_matrix)
final_matrix <- final_matrix[keep,]

#See quantification values and filter
myCPM <- cpm(final_matrix)
lcpm <- cpm(final_matrix,log=T)

#counts por millon reads bigger than 15
treshold <- myCPM>=20
head(treshold)


#True summary
table(rowSums(treshold))
keep <- rowSums(treshold)>=2
summary(keep)
#final_matrix  <- final_matrix [keep,keep.lib.sizes=T]

#Normalization of expression distributions
final_matrix  <- calcNormFactors(final_matrix)

#See quantification values 
myCPM <- cpm(final_matrix )
lcpm <- cpm(final_matrix ,log=T)

#Graphs
plot(myCPM[,1],final_matrix$counts[,1])

#QC
final_matrix$samples$lib.size
barplot(final_matrix$samples$lib.size/1e06,names=colnames(final_matrix),las=2,ann=F,cex.names=0.75)
mtext(side=1,text="Samples",line=4)
mtext(side=2,text="Library size (millions)",line=3)
title("Barplot library size")


ordb <- org.Hs.eg.db

rownames(final_matrix) <- mapIds(ordb,
                              keys=rownames(final_matrix),
                             column="SYMBOL",
                            keytype="ENSEMBL",
                           multiVals="first")

#log2 counts por million - 
logcounts <- cpm(final_matrix,log=T)
boxplot(logcounts,xlab="samples",ylab="log2 counts per million",las=2,outline=F)
abline(h=median(logcounts),col="blue")
title("Logcounts boxplot")

# Scale data to Z-scores
logCPM_z <- t(scale(t(logcounts)))
#Gene 1: Robo1
logcpmROBO1 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="ROBO1")

logcpmROBO1 <- as.data.frame(t(logcpmROBO1))
logcpmROBO1$samples <- rownames(logcpmROBO1)
#Add the group
logcpmROBO1 <- logcpmROBO1 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmROBO1$SampleType <- "Normal"
tumor <- which(logcpmROBO1$group=="t")
logcpmROBO1$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmROBO1 %>%
  ggplot(aes(x = SampleType, y = ROBO1, fill = group)) +
  geom_boxplot() +geom_jitter()
  theme(aspect.ratio = 1) +
  ggtitle("Expression of ROBO1 gene")
logcpmROBO1_nostat <- logcpmROBO1 %>%
  ggplot(aes(x = SampleType, y = ROBO1, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of ROBO1 gene")+labs(y= "ROBO1 log CPM", x = "Tissue")+theme_bw()
logcpmROBO1_nostat
#Add statistics 
library(ggpubr)
logcpmROBO1_stat <- logcpmROBO1_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmROBO1_stat

#classical: ####

##Gata6, ####
logcpmGATA6 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="GATA6")

logcpmGATA6 <- as.data.frame(t(logcpmGATA6))
logcpmGATA6$samples <- rownames(logcpmGATA6)
#Add the group
logcpmGATA6 <- logcpmGATA6 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmGATA6$SampleType <- "Normal"
tumor <- which(logcpmGATA6$group=="t")
logcpmGATA6$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmGATA6 %>%
  ggplot(aes(x = SampleType, y = GATA6, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of GATA6 gene")
logcpmGATA6_nostat <- logcpmGATA6 %>%
  ggplot(aes(x = SampleType, y = GATA6, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of GATA6 gene")+labs(y= "GATA6 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmGATA6_stat <- logcpmGATA6_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmGATA6_stat

##Tff1, ####
logcpmTFF1 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="TFF1")

logcpmTFF1 <- as.data.frame(t(logcpmTFF1))
logcpmTFF1$samples <- rownames(logcpmTFF1)
#Add the group
logcpmTFF1 <- logcpmTFF1 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmTFF1$SampleType <- "Normal"
tumor <- which(logcpmTFF1$group=="t")
logcpmTFF1$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmTFF1 %>%
  ggplot(aes(x = SampleType, y = TFF1, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of TFF1 gene")
logcpmTFF1_nostat <- logcpmTFF1 %>%
  ggplot(aes(x = SampleType, y = TFF1, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of TFF1 gene")+labs(y= "TFF1 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmTFF1_stat <- logcpmTFF1_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmTFF1_stat

##Eps8l3,####
logcpmEPS8L3 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="EPS8L3")

logcpmEPS8L3 <- as.data.frame(t(logcpmEPS8L3))
logcpmEPS8L3$samples <- rownames(logcpmEPS8L3)
#Add the group
logcpmEPS8L3 <- logcpmEPS8L3 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmEPS8L3$SampleType <- "Normal"
tumor <- which(logcpmEPS8L3$group=="t")
logcpmEPS8L3$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmEPS8L3 %>%
  ggplot(aes(x = SampleType, y = EPS8L3, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of GATA6 gene")
logcpmEPS8L3_nostat <- logcpmEPS8L3 %>%
  ggplot(aes(x = SampleType, y = EPS8L3, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of EPS8L3 gene")+labs(y= "EPS8L3 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmEPS8L3_stat <- logcpmEPS8L3_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmEPS8L3_stat

#Ctse,
logcpmCTSE <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="CTSE")

logcpmCTSE <- as.data.frame(t(logcpmCTSE))
logcpmCTSE$samples <- rownames(logcpmCTSE)
#Add the group
logcpmCTSE <- logcpmCTSE %>%
  left_join(final_matrix$samples, by = "samples")

logcpmCTSE$SampleType <- "Normal"
tumor <- which(logcpmCTSE$group=="t")
logcpmCTSE$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmCTSE %>%
  ggplot(aes(x = SampleType, y = CTSE, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of CTSE gene")
logcpmCTSE_nostat <- logcpmCTSE %>%
  ggplot(aes(x = SampleType, y = CTSE, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of CTSE gene")+labs(y= "CTSE log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmCTSE_stat <- logcpmCTSE_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmCTSE_stat

#Nri2, 
logcpmNRI2 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="NRIP2")

logcpmNRI2 <- as.data.frame(t(logcpmNRI2))
logcpmNRI2$samples <- rownames(logcpmNRI2)
#Add the group
logcpmNRI2 <- logcpmNRI2 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmNRI2$SampleType <- "Normal"
tumor <- which(logcpmNRI2$group=="t")
logcpmNRI2$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmNRI2 %>%
  ggplot(aes(x = SampleType, y = NRIP2, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of NRI2 gene")
logcpmNRI2_nostat <- logcpmNRI2 %>%
  ggplot(aes(x = SampleType, y = NRIP2, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of NRIP2 gene")+labs(y= "NRIP2 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmNRI2_stat <- logcpmNRI2_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmNRI2_stat

#Cldn18
logcpmCLDN18 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="CLDN18")

logcpmCLDN18 <- as.data.frame(t(logcpmCLDN18))
logcpmCLDN18$samples <- rownames(logcpmCLDN18)
#Add the group
logcpmCLDN18 <- logcpmCLDN18 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmCLDN18$SampleType <- "Normal"
tumor <- which(logcpmCLDN18$group=="t")
logcpmCLDN18$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmCLDN18 %>%
  ggplot(aes(x = SampleType, y = CLDN18, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of CLDN18 gene")
logcpmCLDN18_nostat <- logcpmCLDN18 %>%
  ggplot(aes(x = SampleType, y = CLDN18, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of CLDN18 gene")+labs(y= "CLDN18 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmCLDN18_stat <- logcpmCLDN18_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmCLDN18_stat

library(cowplot)
title <- ggdraw() + draw_label("Classical subtype genes", fontface='bold')
p <- plot_grid(logcpmGATA6_nostat,logcpmTFF1_nostat,logcpmEPS8L3_nostat,logcpmCTSE_nostat,logcpmNRI2_nostat,logcpmCLDN18_nostat, labels=c("AUTO"))
plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) 

#basal:  ####

##S100a2,####
logcpmS100A2 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="S100A2")

logcpmS100A2 <- as.data.frame(t(logcpmS100A2))
logcpmS100A2$samples <- rownames(logcpmS100A2)
#Add the group
logcpmS100A2 <- logcpmS100A2 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmS100A2$SampleType <- "Normal"
tumor <- which(logcpmS100A2$group=="t")
logcpmS100A2$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmS100A2 %>%
  ggplot(aes(x = SampleType, y = S100A2, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of S100A2 gene")
logcpmS100A2_nostat <- logcpmS100A2 %>%
  ggplot(aes(x = SampleType, y = S100A2, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of S100A2 gene")+labs(y= "S100A2 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmS100A2_stat <- logcpmS100A2_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmS100A2_stat
##Ahnak2,####
logcpmAHNAK2 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="AHNAK2")

logcpmAHNAK2 <- as.data.frame(t(logcpmAHNAK2))
logcpmAHNAK2$samples <- rownames(logcpmAHNAK2)
#Add the group
logcpmAHNAK2 <- logcpmAHNAK2 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmAHNAK2$SampleType <- "Normal"
tumor <- which(logcpmAHNAK2$group=="t")
logcpmAHNAK2$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmAHNAK2 %>%
  ggplot(aes(x = SampleType, y = AHNAK2, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of AHNAK2 gene")
logcpmAHNAK2_nostat <- logcpmAHNAK2 %>%
  ggplot(aes(x = SampleType, y = AHNAK2, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of AHNAK2 gene")+labs(y= "AHNAK2 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmAHNAK2_stat <- logcpmAHNAK2_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmAHNAK2_stat
##Itga3, ####
logcpmITGA3 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="ITGA3")

logcpmITGA3 <- as.data.frame(t(logcpmITGA3))
logcpmITGA3$samples <- rownames(logcpmITGA3)
#Add the group
logcpmITGA3 <- logcpmITGA3 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmITGA3$SampleType <- "Normal"
tumor <- which(logcpmITGA3$group=="t")
logcpmITGA3$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmITGA3 %>%
  ggplot(aes(x = SampleType, y = ITGA3, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression ofITGA3 gene")
logcpmITGA3_nostat <- logcpmITGA3 %>%
  ggplot(aes(x = SampleType, y = ITGA3, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of ITGA3 gene")+labs(y= "ITGA3 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmITGA3_stat <- logcpmITGA3_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmITGA3_stat
##Krt5,####
logcpmKRT5 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="ENSG00000186081")

logcpmKRT5 <- as.data.frame(t(logcpmKRT5))
logcpmKRT5$samples <- rownames(logcpmKRT5)
#Add the group
logcpmKRT5 <- logcpmKRT5 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmKRT5$SampleType <- "Normal"
tumor <- which(logcpmKRT5$group=="t")
logcpmKRT5$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmKRT5 %>%
  ggplot(aes(x = SampleType, y = ENSG00000186081, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of KRT5 gene")
logcpmKRT5_nostat <- logcpmKRT5 %>%
  ggplot(aes(x = SampleType, y = ENSG00000186081, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of KRT5 gene")+labs(y= "KRT5 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmKRT5_stat <- logcpmKRT5_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmKRT5_stat
##Tp63, ####
logcpmTP63 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="ENSG00000073282")

logcpmTP63 <- as.data.frame(t(logcpmTP63))
logcpmTP63$samples <- rownames(logcpmTP63)
#Add the group
logcpmTP63 <- logcpmTP63 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmTP63$SampleType <- "Normal"
tumor <- which(logcpmTP63$group=="t")
logcpmTP63$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmTP63 %>%
  ggplot(aes(x = SampleType, y = ENSG00000073282, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of TP63 gene")
logcpmTP63_nostat <- logcpmTP63 %>%
  ggplot(aes(x = SampleType, y = ENSG00000073282, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of TP63 gene")+labs(y= "TP63 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmTP63_stat <- logcpmTP63_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmTP63_stat
##Muc16####
logcpmMUC16 <- logcounts %>% 
  as.data.frame %>% 
  dplyr::filter(rownames(logcounts)=="ENSG00000181143")

logcpmMUC16<- as.data.frame(t(logcpmMUC16))
logcpmMUC16$samples <- rownames(logcpmMUC16)
#Add the group
logcpmMUC16 <- logcpmMUC16 %>%
  left_join(final_matrix$samples, by = "samples")

logcpmMUC16$SampleType <- "Normal"
tumor <- which(logcpmMUC16$group=="t")
logcpmMUC16$SampleType[tumor] <- "Tumor"

#Plot the boxplot
logcpmMUC16 %>%
  ggplot(aes(x = SampleType, y = ENSG00000181143, fill = group)) +
  geom_boxplot() +geom_jitter()
theme(aspect.ratio = 1) +
  ggtitle("Expression of MUC16 gene")
logcpmMUC16_nostat <- logcpmMUC16 %>%
  ggplot(aes(x = SampleType, y = ENSG00000181143, fill = group)) +
  geom_boxplot() +geom_jitter(color="black")+
  theme(legend.position="none") +
  ggtitle("Expression of MUC16 gene")+labs(y= "MUC16 log CPM", x = "Tissue")+theme_bw()

#Add statistics 
library(ggpubr)
logcpmMUC16_stat <- logcpmMUC16_nostat+stat_compare_means(label = "p.signif", label.x.npc = "centre",method = "t.test")
logcpmMUC16_stat

title <- ggdraw() + draw_label("Basal subtype genes", fontface='bold')
p <- plot_grid(logcpmS100A2_nostat,logcpmAHNAK2_nostat,logcpmITGA3_nostat,logcpmKRT5_nostat,logcpmTP63_nostat,logcpmMUC16_nostat, labels=c("AUTO"))
plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) 
