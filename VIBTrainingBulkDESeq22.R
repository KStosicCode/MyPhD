#Load required packages####
#DESeq2: identification of DE genes
#ggplot2: figures
#pheatmap: to create heat maps

library("DESeq2")
library("ggplot2")
library("pheatmap")


#Set working directory####
#this is the folder where R will look for data and write files
#please adapt this location to match your own computer

setwd("~/training/zelfgegeven/NGS/data")


#Load sample annotation

sampleTable <- read.table("GSE52778_metadata.txt",header=TRUE)
sampleTable


#grouping info: cells and treatment
#every cell line was untreated and treated: paired data

#Locate the data ####
#Define the folder that contains the counts files  
#Create a variable folder that contains the name of this folder
#Counts files are in the htseq_counts folder in the DE folder

folder <- "htseq_counts"


#Import the data ####
#Both counts and metadata need to be stored in a DESeqDataSet
#You also have to define the design: 
#cells is the pairing variable 
#treatment is the grouping variable

dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                  directory=folder,
                                  design=~cells+treatment)


#create background list of genes measured in experiment
#write.table(rownames(dds),"EnsemblAll.txt",quote=FALSE,
#            row.names=FALSE,col.names=FALSE)

#Importing featureCounts() output 
#counts <- read.delim("featureCounts.tabular")
#dds2 <- DESeqDataSetFromMatrix(countData=counts,
#                               colData=sampleTable,
#                               design=~cells+treatment)

#Importing Salmon/Kallisto output
#The airway folder should be placed in your working dir
#It contains a folder salmon with count files in the DE folder
#Create the names of the salmon files:
#salmon + sample name + /quant.sf

salmonfiles <- paste("salmon/",sampleTable$run_accession, 
                     "/quant.sf",sep="")


#Meta data needed for import of salmon data:
#salmon file names have to be in column called "files"
#sample names have to be in column called "names"
#Create new data frame with the 2 required columns

sampleTable2 <- data.frame(files=salmonfiles,
                           names=sampleTable$run_accession)


## Import quantifications on the transcript level
#load tximeta package
#load org.Hs.eg.db

library(tximeta)
library(org.Hs.eg.db)


#Other organisms: Check Bioconductor's annotation packages
#https://bioconductor.org/packages/3.11/data/annotation/

#!SELECT 2 (No use temp) when asked which folder to use!
#DON'T USE 1 (default) you probably have no write access there

st <- tximeta::tximeta(sampleTable2)


#will generate a table of counts per transcript

## Import quantifications on the gene level

sg <- tximeta::summarizeToGene(st)


#will generate a table of counts per gene

## Add gene symbols

sg <- tximeta::addIds(sg,"SYMBOL",gene=TRUE)
sg


#Salmon returns estimated counts, which are not integers
#Need to be converted to integers before they are passed to DESeq2

sampleTable2[3:4] <- sampleTable[3:4]
dds3 <- DESeqDataSetFromMatrix(countData=round(assay(sg)),
                               colData=sampleTable2,
                               design=~cells + treatment)
rm(dds3,sampleTable2,sg,st)


#Go back to htseq counts
#inspects the dds object

dds


#Counts can be retrieved using the counts() function on dds
#look at the manual of counts()

?counts


#The normalized argument allows to retrieve normalized counts
#if normalization was performed

#Kahoot 9: what kind of data type is counts(dds)?

class(counts(dds))
str(counts(dds))


#Kahoot 10: How many of the first 10 genes have 0 counts?

head(counts(dds),n=10)
counts(dds)[1:10,]
assay(dds)[1:10,]


#Kahoot 11: What's the count of ENSG00000004799 in SRR1039509?

counts(dds)["ENSG00000004799",]


#Kahoot 12: How many genes?

nrow(counts(dds))


#Remove genes (rows) with less than 10 counts in total

dds <- dds[rowSums(counts(dds)) >= 10,]


#Kahoot 13: How many genes now?

nrow(counts(dds))


#For gene enrichment analysis
#create background list of genes expressed in tissue

write.table(rownames(dds),"EnsemblTissue.txt",quote=FALSE,row.names=FALSE,
            col.names=FALSE)


#Kahoot 14: What data type is colData?

str(colData(dds))


#show the grouping info for treatment

colData(dds)$treatment


#show the grouping info for cell line

colData(dds)$cells


#relevel to get control as reference

colData(dds)$treatment <- relevel(colData(dds)$treatment,
                                  ref="control")


#Differential expression analysis####
#look at the manual of DESeq()

?DESeq


#Kahoot 19: What's the reduced model in our case?

#DESEq() runs the complete DE pipeline on the raw counts
#The test parameter allows you to use LRT instead of Wald
#Better when you have few replicates or you want to study main effects
#The fitType parameter allows you to use Loess (local) regression
#to estimate the dispersions when you have 3 replicates
#run the full DESeq2 analysis

dds <- DESeq(dds,
             test="LRT",
             full=~cells+treatment,
             reduced=~cells)


#DESeq() result contains the fitted statistical model obtained by:
#normalizing the counts for different library sizes
#estimating the variance for each gene included
#calculating log-fold changes for all included factors
#You can create the dispersion plot

plotDispEsts(dds)


#To see the factors used for library size normalization

sizeFactors(dds)


#If you do not want to do library normalization
#e.g. most genes are expected to be DE
#sizeFactors(dds) <- rep(1,8)
#Better solution would be to subset the DESeqDataSet by housekeeping genes
#Normalize the subset using DESeq()
#Apply the resulting normalization factors back to the full dataset
#sizeFactors(dds) <- vector with size factors of HK subset

#Kahoot 20: Normalized count ENSG00000004799 in SRR1039509?

#the default are the raw counts

counts(dds)["ENSG00000004799",]


#matrix with a normalized counts

counts(dds,normalized=TRUE)["ENSG00000004799",]


#check comparisons that were made

resultsNames(dds)


#the last comparison is comparison of treated to control over all cell types
#the last comparison is the one shown in the results table

#look at the manual of results()

?results


#results() returns log2FC and p vals for last variable in design formula 
# > 2 levels for this variable => it extracts results for 
#comparison of last level over first level
#use contrast argument to define other comparisons
#store results into a new object

res <- results(dds,name="treatment_Dex_vs_control")
head(res,n=10)


#res < metadata with info on the meaning of the columns:

mcols(res,use.names=TRUE)


#ask for a summary of res

summary(res)
summary(results(dds,alpha=0.0001))


#Kahoot 2: Which gene has lowest padj?

head(res[order(res$padj),])


#Most people are inclined to select DE genes based on lfc
#Kahoot 3: How many genes with lfc > 1 (> 2-fold up)?

sum(res$log2FoldChange > 1)


#That's not a good idea
#Select on padj or on padj and lfc
#Kahoot 4: How many genes with lfc > 1 and padj < 0.01?

sum((res$log2FoldChange > 1)&(res$padj < 0.01))
#NA su problem
sum((res$log2FoldChange > 1)&(res$padj < 0.01),na.rm=TRUE)

#Shrink lfcs ####
#in current versions of DESEq2 shrinkage is not done by the DESeq function
#you can calculate them afterwards using lfcShrink()

resultsNames(dds)
res.shr <- lfcShrink(dds=dds,coef=5)


#Kahoot 12: Original (unshrunken) lfc of ENSG00000004799?
#Kahoot 13: Shrunken lfc of ENSG00000004799?

res["ENSG00000004799",]
res.shr["ENSG00000004799",]


#Plot the counts for one gene####
#plot the counts for gene ENSG00000004799

plotCounts(dds,gene="ENSG00000004799",intgroup="treatment",
           main="ENSG00000004799")


#nicer version using ggplot2

d <- plotCounts(dds,gene="ENSG00000004799", 
                intgroup="treatment", 
                main="ENSG00000004799", 
                returnData=TRUE)
d


#ggplot2 can handle the plots generated by DESeq2

ggplot(d,aes(treatment,count)) + 
  geom_point() +
  scale_y_continuous(trans="log10",breaks=c(25,100,400)) + 
  labs(title="ENSG00000004799")


#make same plot using ggplot but this time
#color the dots according to cell line
#create d again 
#Kahoot 14: What do you need to change?

d <- plotCounts(dds,gene="ENSG00000004799", 
                intgroup=c("treatment","cells"), 
                main="ENSG00000004799",returnData=TRUE)
  
ggplot(d,aes(treatment,count,color=cells)) +
  geom_point() + 
  scale_y_log10(breaks=c(25,100,400)) + 
  labs(title="ENSG00000004799")


#Histogram of shrunken log fold changes####
#From res.shr remove rows with missing adjusted p-values
#Use ggplot2 to make a histogram with 60 bins
#Remember: the input for ggplot2 needs to be a data frame
#remove rows with missing padj
#Kahoot 15: How many genes were removed?

sum(is.na(res.shr$padj))
resFix <- res.shr[!is.na(res.shr$padj),]


#input for ggplot2 needs to be a data frame

resPlot <- as.data.frame(resFix)


#create the histogram with ggplot2

ggplot(resPlot,aes(x=log2FoldChange)) + geom_histogram(bins=60)


#Add gene symbols to data####

library(plyr)


#Make data frame with Ensembl IDs and symbols 
#for all human genes

keys <- keys(org.Hs.eg.db)
head(keys)
s <- AnnotationDbi::select(org.Hs.eg.db,keys=keys,
                           columns=c("ENSEMBL","SYMBOL"))
head(s)


#Remove rows with missing Ensembl IDs

s <- s[!is.na(s$ENSEMBL),c("ENSEMBL","SYMBOL")]
head(s)


#Create column in resPlot with Ensembl IDS (now row names)

resPlot$ENSEMBL <- rownames(resPlot)
head(resPlot)


#Join both tables based on Ensembl IDs

resPlot <- join(resPlot,s,by=c("ENSEMBL"))  


#Remove duplicated rows

resPlot <- resPlot[!duplicated(resPlot$ENSEMBL),]


#Kahoot 16: What is the gene symbol of ENSG00000004799?

head(resPlot)
resPlot[resPlot$ENSEMBL=="ENSG00000004799",]


#Select DE genes ####
#summary(res) showed that there are many DE genes at FDR of 10%
#But 10% FALSE positives is quite high
#Select genes with adjusted p-value < 0.01 (FDR = 1%)

#Often an additional threshold is set on lfc
#Expression needs to change to call gene DE
#Select genes with shrunken lfc at least 1 or -1 (2-fold)
#Note that TFs often show lower lfc but have huge impact

#Stringency depends on follow up plans:
#Experimental study of a few DE genes: very strict
#Pathway analysis: not so strict
#Kahoot 17: How many DE genes according to these thresholds?

resString <- resPlot[abs(resPlot$log2FoldChange) >= 1 &
                       (resPlot$padj > 0.01),]
nrow(resString)
head(resString)

 #Volcano plot with ggplot2 ####
#create a column to use for the coloring
#color dots with padj < e-30 and abs(lfc) > 2 red
#I am going to take extreme thresholds because I want to label

resPlot$threshold <- abs(resPlot$log2FoldChange)>2 & resPlot$padj < exp(-30)
sum(resPlot$threshold)


#Create the Volcano plot with ggplot2: use hollow points 
#change colors to black and red
#change legend title to "lfc>2 and p<1e-10"
#change plot title to "p-value versus fold change"
#add horizontal line at p=1e-10
#add labels of red genes
library(ggrepel)
ggplot(resPlot,aes(log2FoldChange,-log10(padj),color=threshold)) + 
  geom_point(shape=1) +
  scale_color_manual(name="lfc > 2 and p < 1e-10",
                     values=c("black","red")) +
  ggtitle("p-value versus fold change") + 
  geom_hline(yintercept=10,color="red",linetype=2) +
  geom_vline(xintercept=-2,color="red",linetype=2) +
  geom_vline(xintercept=2,color="red",linetype=2) +    
  geom_text_repel(aes(label=ifelse(threshold==TRUE,SYMBOL,"")),size=2)

#Different color for down and up - no labels
#Use regular thresholds: 2-fold and FDR=1%
#Create a logical operation for UP and one for DOWN
#Kahoot 18: How many down genes?
logic1 <- 
logic2 <- 

resPlot$threshold <- ifelse(logic1,'B',ifelse(logic2,'C','A'))

ggplot(resPlot,aes(log2FoldChange,-log10(padj))) + 
  geom_point(shape=1,aes(color=threshold)) +
  scale_color_manual(name="differential expression",
                     values=c("black","red","green"),
                     labels=c("|lfc|<2 or p>e-10",
                              ">4 fold UP and p<e-10",
                              ">4 fold DOWN and p<e-10")) +
  ggtitle("p-value versus fold change") + 
  geom_hline(yintercept=10,color="red",linetype=2)+
  geom_vline(xintercept=-2,color="red",linetype=2) +
  geom_vline(xintercept=2,color="red",linetype=2)

#PCA plot####
#transform the data:
#use DESeq2 variance stabilizing transformation
#Variance grows with mean
#If one performs PCA directly on size-factor-normalized read counts, 
#the result depends only on the few most strongly expressed genes 
#because they show the largest absolute differences between samples. 
#A strategy to avoid this is to take the log of normalized counts 
#plus a small pseudocount. 
#But now genes with the lowest counts tend to dominate the results 
#because the log amplifies differences for the smallest values
#Solution: DESeq2 offers transformations that stabilize the variance
#the rlog and the variance stabilizing transformation (VST)
#We will use the VST implemented with the vst() function:
vsd <- vst(dds,blind=FALSE)
class(vsd)
#Is sample annotation still there? 
colData(vsd) 
#Kahoot 19: Vst count of ENSG00000004799 in SRR1039509?

#The samples are projected onto the 2D plane such that they spread 
#out in the two directions that explain most of the differences 
data <- plotPCA(vsd,intgroup=c("treatment","cells"),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

#make the plot with ggplot using data as data:
#it's a scatter plot of PCA2 (y) vs PCA1 (x) use geom_point()
#color points according to treatment
#shape points according to cells
#label X-axis: PC1: percentVar[1] % variance: use labs()
#label Y-axis: PC1: percentVar[2] % variance 
#Kahoot 20: What factor determines PC1?
#Kahoot 21: What if cells determines PC1?

#Kahoot 22: If you see an outlier sample, can you remove it?
#The x-axis is the direction that separates the data the most. . 
#The y-axis is a direction that separates the data the second most.
#The percent of the total variance is printed in the axis label. 
#These percentages do not add to 100%
#because there are more dimensions that contain the remaining variance

#color PCA according to expression of ENSG00000165995
data$expr <- assay(vsd)["ENSG00000165995",]
ggplot(data,aes(PC1,PC2,color=expr,shape=treatment)) + geom_point(size=3) +
  labs(x=paste0("PC1: ",percentVar[1],"% variance"),
       y=paste0("PC2: ",percentVar[2],"% variance"))

#The differences between cells (different shapes) are considerable,
#though not stronger than the differences due to the treatment.
#This shows why it is important to account for the cell line effect

#Which genes contribute most to PCs?
#Do PCA yourself using prcomp
#input: vst counts only - genes are columns
pcs <- prcomp(t(assay(vsd)))
head(pcs$rotation)
#Look for genes with highest/lowest contribution = markers
#Kahoot 23: Which gene has the highest contribution to PC1?

x <- rownames(assay(vsd))[which.min(pcs$rotation[,1])]
resPlot[resPlot$ENSEMBL==x,]

#Heat map ####
#Make heat map only for DE genes 
#Simplify for training: 
#Create a vector with the ENSEMBL IDs of the 10 genes with lowest padj
#order resPlot according to increasing padj
resOrdered <- resPlot[order(resPlot$padj),]
#view 10 most significant genes
resOrdered[1:10,]
#select gene names of 10 most significant genes
select <- resOrdered$ENSEMBL[1:10]
select

#Retrieve vst values for selected genes
#vst values are already log transformed
toplot <- assay(vsd)[select,]

#for selection you need Ensembl IDs: they are unique
#for heat map you want gene symbols because they are meaningfull
rownames(toplot) <- resOrdered$SYMBOL[1:10]

#make heat map of normalized counts for these genes
df <- as.data.frame(colData(vsd)[,c("treatment","cells")])
ann_colors <- list(treatment=c(control="blue",Dex="red"),
                   cells=c(N61311="aquamarine",N052611="blue",
                           N080611="darkgrey",N061011="azure"))

#use pheatmap() to make the heat map of the toplot data
#show gene symbols for the rows
#use df as annotation for the columns
#use ann_colors as colors for the column annotations
#for the heat map low expression=green, high expression=red, medium=black
pheatmap(toplot,annotation_col=df,annotation_colors=ann_colors,
         color=colorRampPalette(c("green","black","red"))(50))
#Kahoot 24: Which gene is the least like the others?

toplot
#You look at absolute expression (normalized + variance stabilized + log)
#ADAMTS1 is indeed the gene with the highest vst expression
#You can also look at expression pattern (up/down) by scaling the genes
#mean = 0 and SD = 1 to make genes comparable
?pheatmap
#Kahoot 25: Which gene is the odd one out if you scale the genes?
pheatmap(toplot,annotation_col=df,annotation_colors=ann_colors,
         ?,border_color="transparent",
         clustering_method="average",
         cutree_cols=2,cutree_rows=2,fontsize=7,
         color=colorRampPalette(c("green","black","red"))(50))
toplot

#What if you want to know which genes belong to each cluster?
#Do the clustering yourself
#Data = matrix with scaled genes: use scale()
#Calculate distances between scaled proteins: use dist()
?scale
?dist
toplots <- t(scale(t(toplot)))
#Calculate Euclidean distances
diss <- dist(toplots)
#Hierarchical clustering based on Euclidean distances
h <- hclust(diss,method="average")
#pheatmap performs by default complete clustering
#average often gives better results than complete
#change clustering of pheatmap by clustering_method="average"
#Make dendrogram
library(dendextend)
dend <- as.dendrogram(h)
#Find clusters
cutree(dend,k=2)

#Plot expression patterns and color per cluster ####
#If you want to plot unscaled vsts use toplot instead of toplots
#Add cluster numbers + symbols to data
toplotdf <- as.data.frame(toplots)
toplotdf$cluster <- cutree(dend,k=2)
toplotdf$symbol <- row.names(toplots)

#Convert data with cluster info to long format
library(reshape2)
toplotl <- melt(toplotdf,id.vars=9:10,measure.vars=1:8)
#Make a line plot: you need to tell R how to draw the lines
ggplot(toplotl,aes(variable,value,color=factor(cluster),group=symbol)) + 
  geom_line() + theme(axis.text.x=element_text(angle=90)) +
  theme(axis.title.x=element_blank()) +
  labs(y="scaled vst",color="cluster")

#if you want to plot FPKM values = normalized for transcript length
#you need salmon/kallisto data or gene annotation in rowRanges(dds)
#robust: DESeq library size (TRUE) or total count normalization (FALSE)
#fpkm(dds,robust=TRUE)
#DO NOT USE THEM FOR DESeq() !

#Save results ####
#write resString to a file outputRNASeq.txt
write.table(resString,file="outputRNASeq.txt",quote=FALSE)

#Write Ensembl IDs of DE genes
write.table(rownames(resString),file="outputNames.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

#ENSEMBL IDs of upregulated genes
EnsemblIDs <- rownames(resString[resString$log2FoldChange > 0,])
write.table(EnsemblIDs,file="outputupEnsembl.txt",quote=FALSE,
            col.names=FALSE,row.names=FALSE)
#Gene symbols of upregulated genes
Symbols <- na.omit(resPlot[resPlot$ENSEMBL %in% EnsemblIDs,"SYMBOL"])
write.table(Symbols,file="outputupNames.txt",quote=FALSE,col.names=FALSE,
            row.names=FALSE)