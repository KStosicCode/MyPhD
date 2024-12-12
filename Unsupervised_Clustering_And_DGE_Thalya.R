library(edgeR)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(stats)
library(dplyr)
library(ggplot2)
library(ggfortify)

setwd("Documents/RNA-seq_ULB/Christelle_Data/mFOLFIRINOX/")
setwd("C:/Users/equertin/Desktop/R Studio for Bulk and Single Cell Analysis/Bulk")

#Importing table with counts
table <- read.table("LCM_count_table.txt")
pattern <- "X"
str_detect(colnames(table), pattern)
colnames(table) <- str_replace(colnames(table),pattern,"")

#Importing Thalya's table
sampleTable2 <- read.csv2("REP-NON-REP-binarisé.csv")

#Changing the names of the values 0 and 1 into NR (non-responders) and R (responders)
index=which(sampleTable2$Conclusion==0)
sampleTable2$Conclusion[index]="NR"
index=which(sampleTable2$Conclusion==1)
sampleTable2$Conclusion[index]="R"
#Changing the name of the first column of Thalya's table
colnames(sampleTable2)[1] <- "patient_id"
#Discarding empty rows in Thalya's table
sampleTable2 <- sampleTable2[1:62,]

#17H14853 (1R,1?), DISCARD
index <- which(sampleTable2$patient_id=="17H14853")
sampleTable2 <- sampleTable2[-index,]
#19H07206 only take 1
index <- which(duplicated(sampleTable2$patient_id)==TRUE)
sampleTable2 <- sampleTable2[-index,]

#Importing Christelle's table
#metadata <- read.csv2("LCM_samplesheet.csv",header=T)
metadata <- read.delim("LCM_samplesheet.txt")
index <- which(metadata$responder==levels(as.factor(metadata$responder))[1])
metadata$responder[index]<-   "R"
index <- which(is.na(metadata$responder))
metadata$responder[index]<-   "NoData"
index <- which(is.na(metadata$treatment))
metadata$treatment[index]<-   "No_treated"
colnames(metadata)[4] <- "treatment"
colnames(metadata)[6] <- "sample_type"
metadata <- dplyr::rename(metadata,samples = sample)

#Creating DGEList-Object
count_matrix <- DGEList(table,samples=colnames(table),group=metadata$responder)
#Joining/merging together count_matrix$samples with the metadata
count_matrix$samples <- merge(metadata,count_matrix$samples,by="samples")

#Filtering
keep <- filterByExpr(count_matrix)
count_matrix <- count_matrix[keep,]

#See quantification values and filter
myCPM <- cpm(count_matrix)
lcpm <- cpm(count_matrix,log=T)

#counts per million reads bigger than 10
threshold <- myCPM>=10
head(threshold)

#True summary
table(rowSums(threshold))
keep <- rowSums(threshold)>=2
summary(keep)
dim(count_matrix$counts)
count_matrix <- count_matrix[keep,keep.lib.sizes=T]
dim(count_matrix$counts)

#Normalization of expression distributions
count_matrix  <- calcNormFactors(count_matrix)

#Removal of the samples with high rRNA share
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
rRNApersample <- data.frame()
expr_data <- count_matrix$counts
view(expr_data)
genes <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                       "transcript_biotype"), 
                        values = rownames(count_matrix), mart = mart)
for(index in 1:ncol(expr_data)){
  sample <- expr_data[,index]
  keep<- as.data.frame(sample)>0
  sample  <- sample[keep]
  #Gene classifications 
  keep <- which(rownames(as.data.frame(sample))%in%genes$ensembl_gene_id)
  gene_post<- genes[keep,]
  sum <- length(which(gene_post$transcript_biotype=="rRNA"))+
    length(which(gene_post$transcript_biotype=="Mt_rRNA"))+
    length(which(gene_post$transcript_biotype=="rRNA_pseudogene"))
  rRNApersample[index,1] <- sum
}
View(rRNApersample)

#Outliers
rRNApersample <- as.numeric(rRNApersample$V1)
mean <- mean(rRNApersample)
std <- sd(rRNApersample)

#Threshold values for outliers
Tmax <- mean+(2*std)

#Find outliers
outliers_values <- rRNApersample[which(rRNApersample>Tmax)]
outlier_pos <- which(rRNApersample>Tmax)

#Remove outliers
count_matrix <- count_matrix[,-outlier_pos]

#Giving genes their names instead of Ensembl IDs
ordb <- org.Hs.eg.db
rownames(count_matrix) <- mapIds(ordb,
                                 keys=rownames((count_matrix)),
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first")
head(rownames(count_matrix$counts))
view(count_matrix$counts)
view(count_matrix$samples)
#Getting rid of the genes without name (missing values)
index <- which(is.na(rownames(count_matrix)))
count_matrix$counts <- count_matrix$counts[-index,]
nrow(count_matrix$counts)
#dim(count_matrix$counts)
#dim(count_matrix$samples)
ncol(count_matrix$counts)
nrow(count_matrix$samples)

#See quantification values 
myCPM <- cpm(count_matrix)
lcpm <- cpm(count_matrix,log=T)

#Graphs
plot(myCPM[,1],count_matrix$counts[,1])

#QC
count_matrix$samples$lib.size
barplot(count_matrix$samples$lib.size/1e06,names=colnames(count_matrix),las=2,ann=F,cex.names=0.75)
mtext(side=1,text="Samples",line=4)
mtext(side=2,text="Library size (millions)",line=3)
title("Barplot library size")

#log2 counts per million 
logcounts <- cpm(count_matrix,log=T)
boxplot(logcounts,xlab="samples",ylab="log2 counts per million",las=2,outline=F)
abline(h=median(logcounts),col="blue")
title("Logcounts boxplot")

#Joining together Thalya's table and count_matrix$samples 
count_matrix$samples <- left_join(count_matrix$samples, sampleTable2, by = "patient_id",keep=F,multiple="all")
index <- which(count_matrix$samples$Conclusion=="?")
count_matrix$samples$Conclusion[index] <- "NotKnown"
index <- which(is.na(count_matrix$samples$Conclusion))
count_matrix$samples$Conclusion[index] <- "NotKnown"

#Filtering
#Filtering out rows that have "n" in the sample_type column (and keep those that have "t")
count_matrix$samples$sample_type
index = which(count_matrix$samples$sample_type=="t")
count_matrix <- count_matrix[,index]
count_matrix$samples$sample_type
ncol(count_matrix$counts)
nrow(count_matrix$samples)
length(count_matrix$samples$sample_type)

#Filtering out rows that don't have "gratter" in the type column (and keep those that have "punch")
count_matrix$samples$type
index = which(count_matrix$samples$type=="gratter")
count_matrix <- count_matrix[,index]
count_matrix$samples$type
ncol(count_matrix$counts)
nrow(count_matrix$samples)
length(count_matrix$samples$type)

#Filtering out rows that have "normal" in the Zone column (and keeping the rest)
#count_matrix$samples$zone
#index=which(count_matrix$samples$zone=="normal")
#count_matrix <- count_matrix[,-index]
#count_matrix$samples$zone

#Filtering out rows that don't have "R" or "NR" in the Conclusion column 
#count_matrix is a DGEList-Object, count_matrix$samples is a dataframe and
#count_matrix$samples$Conclusion is a column (vector)
count_matrix$samples$Conclusion
index=which(count_matrix$samples$Conclusion=="R" | count_matrix$samples$Conclusion=="NR")
count_matrix <- count_matrix[,index]
count_matrix$samples$Conclusion
ncol(count_matrix$counts)
nrow(count_matrix$samples)
length(count_matrix$samples$Conclusion)

#Filtering out rows that don't have mFOLFIRINOX in the Type.chimio table
#count_matrix$samples$Type.chimio
#index = which(count_matrix$samples$Type.chimio=="mFOLFIRINOX")
#count_matrix <- count_matrix[,index]
#count_matrix$samples$Type.chimio
#ncol(count_matrix$counts)
#nrow(count_matrix$samples)
#length(count_matrix$samples$Type.chimio)

#Filtering out rows that have "Gem" or "gemox" in the treatment column (and keep
#mFOLFIRINOX and FOL)
#count_matrix$samples$treatment
#index=which(count_matrix$samples$treatment!="Gem" & count_matrix$samples$treatment!="gemox")
#count_matrix <- count_matrix[,index]
#count_matrix$samples$treatment
#ncol(count_matrix$counts)
#nrow(count_matrix$samples)
#length(count_matrix$samples$treatment)

#Filtering out rows that belong to the same patient (keep only one sample for each patient)
count_matrix$samples$patient_id
unique(count_matrix$samples$patient_id)
#Removing the duplicates for NR (Finding out the index position for each of the  
#five Non-Responders in order to see if they have any duplicates and if yes, to 
#remove those duplicates)
which(count_matrix$samples$patient_id=="18H00665")
which(count_matrix$samples$patient_id=="19H10569")
index <- which(count_matrix$samples$samples=="269" | count_matrix$samples$samples=="270"
               | count_matrix$samples$samples=="272")
count_matrix <- count_matrix[,-index]
which(count_matrix$samples$patient_id=="19H10569")
which(count_matrix$samples$patient_id=="19H10871")
index <- which(count_matrix$samples$samples=="235" | count_matrix$samples$samples=="238"
               | count_matrix$samples$samples=="239")
count_matrix <- count_matrix[,-index]
which(count_matrix$samples$patient_id=="19H10871")
which(count_matrix$samples$patient_id=="20H00014")
which(count_matrix$samples$patient_id=="20H04238")
index <- which(count_matrix$samples$samples=="307")
count_matrix <- count_matrix[,-index]
which(count_matrix$samples$patient_id=="20H04238")
#Removing the rest of the duplicates (the duplicates for Responders)
unique(count_matrix$samples$patient_id)
count_matrix$samples$patient_id
index=which(duplicated(count_matrix$samples$patient_id)==TRUE)
count_matrix <- count_matrix[,-index]
count_matrix$samples$patient_id
ncol(count_matrix$counts)
nrow(count_matrix$samples)
length(count_matrix$samples$patient_id)
count_matrix$samples$Conclusion
which(count_matrix$samples$Conclusion=="NR")

#Creating a design (model) matrix
group <- factor(count_matrix$samples$Conclusion) 
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
view(design)

#Creating a contrast matrix according to Responders VS Non-Responders  
cont.matrix <- makeContrasts(RespVsNresp=R - NR,
                             levels=design)
dim(cont.matrix)

#Voom
par(mfrow=c(1,1))
v <- voom(count_matrix, design, plot = TRUE)
names(v)
v$targets
v$E
view(v$E)
dim(v$E)
v$weights
v$design
v

#Data for the unsupervised learning
mydata=t(v$E)
mydata<-scale(mydata)

#Choosing optimal number of clusters - option 1 - Silhouette
library(factoextra)
fviz_nbclust(mydata, kmeans, method = "silhouette", print.summary =FALSE,verbose =FALSE)

#Choosing optimal number of clusters - option 2 - WSS Plot 
#WSS Plot - option A - same as Silhouette, only put "wss" for the method argument
#fviz_nbclust(mydata, kmeans, method = "wss", print.summary =FALSE,verbose =FALSE)

#WSS Plot - option B - Defining/Creating wssplot() function for use 
#wssplot <- function(data, nc=15, seed=1234)
#{
#wss <- (nrow(data)-1)*sum(apply(data,2,var))
#for (i in 2:nc){
#set.seed(seed)
#wss[i] <- sum(kmeans(data, centers=i)$withinss)}
#plot(1:nc, wss, type="b", xlab="Number of Clusters",
#ylab="Within groups sum of squares")
#wss
#}

#wssplot(mydata)

#K-Means Cluster
KM=kmeans(mydata, 2)

#Evaluating Cluster Analysis - option 1 - Cluster Plot
autoplot(KM, mydata, frame=TRUE)

#Evaluating Cluster Analysis - option 2 - Cluster Centres
#KM$centers
#head(KM$centers)

#Getting the information which sample belong to which cluster
summary(KM)
head(KM$cluster)
KM$cluster
index1=which(KM$cluster==1)
index2=which(KM$cluster==2)

#Assigning the values/results from Unsupervised clustering to a new column inside
#the count_matrix$samples
count_matrix$samples$Cluster <- KM$cluster
view(count_matrix$samples)

# Scatter plot of mydata
#plot(mydata, 
#col = KM$cluster,
#main = "k-means with 2 clusters",
#xlab = "",
#ylab = "")

#Fitting the data into lm model
fit <- lmFit(v,design)
names(fit)
coef.fit <- fit$coefficients
head(coef(fit))
head(coef.fit)

#Creating a contrast matrix
fit.cont <- contrasts.fit(fit, cont.matrix)

#Empirical Bayes Statistics for Differential Expression and removal of the false positives,
#by borrowing information across all the genes to obtain more precise estimates of
#gene-wise variability
fit.cont <- eBayes(fit.cont)
plotSA(fit.cont, main="Final model: Mean-variance trend")
dim(fit.cont)
summa.fit <- decideTests(fit.cont, p.value  = 0.05, lfc=1)
summary(summa.fit)

#Extracting a general table of the top-ranked genes from a linear model fit and sorting them by p-values
top.table <- topTable(fit.cont, sort.by = "P", n = Inf) 
top.table<- topTable(fit.cont,coef=NULL,number=Inf,sort.by = "B",resort.by="logFC")
head(top.table, 5) 
#Getting rid of NAs in top.table
index <- which(is.na(top.table$ID))
top.table <- top.table[-index,]
head(top.table, 5) 
class(top.table)

#Saving the general table (top.table)
save(top.table,file="R_vs_NR_top.table_Thalya_DEG.Rdata")
write.csv(x=top.table,"R_vs_NR_top.table_Thalya_DEG.csv")
write.table(top.table, file = "R_vs_NR_top.table_Thalya_DEG.txt", row.names = F, sep = "\t", quote = F)
#Importing the saved table
w <- read.delim("R_vs_NR_top.table_Thalya_DEG.txt")

#Up and down regulated genes in R vs NR
R_vs_NR <- topTable(fit.cont,coef=NULL,number=Inf,sort.by = "logFC")
library(dplyr)
R_vs_NR <- dplyr::filter(R_vs_NR, (logFC>1|logFC< -1)&P.Value<0.05)
#Take only those that have no missing values
R_vs_NR <-  R_vs_NR[complete.cases(R_vs_NR), ]
#Up regulated
R_vs_NR_up <- R_vs_NR[R_vs_NR$logFC>0,]
#Down regulated
R_vs_NR_down <- R_vs_NR[R_vs_NR$logFC<0,]

save(R_vs_NR_up,file="R_vs_NR_up_Thalya.Rdata")
write.csv(x=R_vs_NR_up,"R_vs_NR_up_Thalya.csv")
write.table(R_vs_NR_up, file = "R_vs_NR_up_Thalya.txt", row.names = F, sep = "\t", quote = F)
save(R_vs_NR_down,file="R_vs_NR_down_Thalya.Rdata")
write.csv(x=R_vs_NR_down,"R_vs_NR_down_Thalya.csv")
write.table(R_vs_NR_down, file = "R_vs_NR_down_Thalya.txt", row.names = F, sep = "\t", quote = F)

#sampleTable2 <- read.csv2("REP-NON-REP-binarisé.csv")
#sampleTable2 <- sampleTable2[1:62,]
#view(sampleTable2)
#metadata <- read.delim("LCM_samplesheet.txt")
#view(metadata)
#sampleTableJ <- read.csv2("JulieTable4.csv",header=T,dec=c(",","."),
#skip=1)
#view(sampleTableJ)

#index=which(sampleTable2$Num.anapath%in%metadata$patient_id)
#o <- sampleTable2[-index,]
#which(o$patient_id%in%metadata$patient_id)

#which(o$Num.Erasme%in%sampleTableJ$ID.Erasme)
