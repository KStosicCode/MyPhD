library("DESeq2")
library("ggplot2")
library("pheatmap")

setwd("C:/Users/equertin/Desktop/R Studio for Bulk and Single Cell Analysis/Bulk")

sampleTable <- read.delim("LCM_samplesheet.txt")

counts <- read.delim("LCM_count_table.txt") #pogledaj dvadeseti red (rownames od counts)

geneTable <- read.delim("ensembl_100_eid_symbol_conversion.txt")
?read.delim

sampleTable$treatment[is.na(sampleTable$treatment)] <- "untreated"
sampleTable$radiation[is.na(sampleTable$radiation)] <- "untreated"

sampleTable[is.na(sampleTable)] <- "unknown" 

rownames(geneTable)<-geneTable[,"Gene.stable.ID"]
geneTable.sorted<-geneTable[rownames(counts),] #pogledaj deveti red (counts je assignment)


dds<-DESeqDataSetFromMatrix(countData=counts,
                            colData=sampleTable,
                            design=~treatment)


dds <- dds[rowSums(counts(dds)) >= 40, ] 

colData(dds)$treatment <- relevel(colData(dds)$treatment,
                                  ref="mFOLFIRINOX")


dds <- DESeq(dds)

plotDispEsts(dds)

res <- results(dds,contrast=c("treatment","mFOLFIRINOX","Gem"))

resultsNames(dds)

res.shr <- lfcShrink(dds=dds,coef=2)

resFix <- res.shr[!is.na(res.shr$padj), ]

resPlot <- as.data.frame(resFix)


library("plyr")
library("org.Hs.eg.db")

keys <- keys(org.Hs.eg.db)
head(keys)
s <- AnnotationDbi::select(org.Hs.eg.db,keys=keys,
                           columns=c("ENSEMBL","SYMBOL"))
head(s)
s <- s[!is.na(s$ENSEMBL),c("ENSEMBL","SYMBOL")]

resPlot$ENSEMBL <- rownames(resPlot)

resPlot <- join(resPlot,s,by=c("ENSEMBL"))

resPlot <- resPlot[!duplicated(resPlot$ENSEMBL), ]


resString <- resPlot[abs(resPlot$log2FoldChange) >= 1 &
                       (resPlot$padj < 0.01), ]
resString2 <- resPlot[abs(resPlot$log2FoldChange) >= 1 &
                        (resPlot$padj) < 0.05, ]


vsd <- vst(dds,blind=FALSE)


#head(geneTable)
#head(resPlot)
#head(counts)
#head(geneTable.sorted)
#dim(geneTable.sorted)
#dim(resPlot)
#resPlot.info<-cbind(resPlot,geneTable.sorted[rownames(resPlot),])