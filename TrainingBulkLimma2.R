if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")

library(edgeR)
library(limma)

title: "limma PipeLine"
author: "Brandon YEO"
date: '2022-04-08'
output:
  md_document:
  variant: gfm
---
  
  Chapters 
1. Data import and understanding of data 
2. Data normalization and making comparison 
3. DEGs isolation from groups 

In General, the pipeline looks like this
Step 1 - DGEList() 
Step 2 - model.matrix()
Step 3 - voom()
Step 4 - lmfit()
Step 5 - makeContrasts()  + Contrasts.fit()

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# Based on #https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# PSA, don't use single variable object name 
```

```{r}
library(edgeR)
library(limma)
```

# Chapter 1 Data import and Understanding of Data 

Importovanje matrice
```{r}
#Postoje dva pristupa - DGEList() i readDGE() - DGEList() omogucava ucitavanje matrice,
#tj. koristi se kada vec postoji matrica gde su sva ocitavanja (counts) ubacena zajedno i zelis
#da takvu matricu koja vec sadrzi sve (, a u prethodnom koraku je importovana uz pomoc 
#read.delim()) konvertujes u DGEList-object, dok u slucaju da su ocitavanja (counts) zasebna
#(kao sto je Nikola poslao za Zulin projekat), onda se svi ucitaju tako sto se iskoristi 
#readDGE()

#Odnosno, ako su ocitavanja (counts) poreklom od svih uzoraka "uskladistena" u jedan fajl,
#ta jedna tabela koja ih sve sadrzi moze da se ucita sa read.delim() i potom konvertuje u
#DGEList-object funkcijom DGEList(), dok ako postoji zaseban counts fajl za svaki uzorak,
#koristi se readDGE() i direktno nastaje DGEList-object
counts <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt")
head(counts)
```

Create DGE object 
```{r}
#DGEList() je kao DESeqDataSetFromMatrix()
d0 <- DGEList(counts) 

#Ovde dolazi do izracunavanja faktora normalizacije, koji se potom koriste kao skalirajuci
#faktor za razlike velicini biblioteka, kao kod DESeq()
d0 <- calcNormFactors(d0) 

Filter data on lower count rate
```{r}
#prag u vidu ekspresije jednom na milion counts
cutoff <- 1 

head(d0)
#cpm znaci counts per million i predstavlja jednu od mogucih transformacija raw counts prema
#skali koja uzima u obzir razlike u velicinama biblioteka, a jos neke od mogucih transformacija
#koje su u upotrebi su log-cpm, rpkm i fpktm
head(cpm(d0)) 

#1 kod apply() znaci da se max primeni samo na redove
drop <- which(apply(cpm(d0),1,max) < cutoff) 

#otarasice se indeksa koji su smesteni u drop u prethodnom koraku, a tada su izabrani oni redovi
#koji imaju jako slabu ekspresiju gena i dodeljeni su kao vrednost za drop;
#iako je drugaciji pristup, rezultat je isti kao dds[rowSums(counts(dds)) >= od nekog broja, ] 
d <- d0[-drop,] 

# number of genes left   
dim(d)              
```

set sample names             
```{r}
#Sample names
snames <- colnames(counts) 
snames
cultivar <- substr(snames, 1, nchar(snames) - 2)
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar 
```

```{r}
time
```

```{r}
#ovako oformljena grupa koja kombinuje cultivar i time moze da se iskoristi za update grupe 
#unutar samples unutar DGEList-objekta po imenu d, a to bi izgledalo ovako
#d[["samples"]][["group"]] <- group ili ovako 
#d$samples$group <- group 
group <- interaction(cultivar, time) 
group 
```

```{r}
plotMDS(d, col = as.numeric(group),cex = 1.5)
```
- Checking with PCA 
- https://www.quora.com/Whats-the-difference-between-MDS-and-PCA?share=1
- https://stats.stackexchange.com/questions/14002/whats-the-difference-between-principal-component-analysis-and-multidimensional
```{r}
library(factoextra)
library(FactoMineR)
pca.raw.d <- log2(d$counts+0.5)
pca.d <- PCA(t(pca.raw.d),graph = F)
fviz_pca_ind(pca.d, col.ind = group)
```
# Chapter 2 - Data Normalization and making models 

- Voom Transformation and calculation of variance 
- Transform RNA-Seq Data Ready for Linear Modelling

https://www.montana.edu/rotella/documents/502/DesignMatricesR.pdf
```{r}
library(pheatmap)
#model.matrix() kreira matricu modelovanu prema odgovarajucem faktoru ili faktorima, ali tako da
#kao vrednosti dodeljuje samo 0 i 1 (dummy variable, a i ima smisla, ako se uzme u obzir da se
#vektor sa faktorizovanim vrednostima/nivoima/levelima konvertuje u matricu), u zavisnosti od
#toga da li odgovarajuci uzorak pripada ili ne pripada doticnoj grupi
mm <- model.matrix(~ 0 + group)
#Drugi naziv za model matricu je design matrix
pheatmap(mm,cluster_rows = FALSE,cluster_cols = FALSE)
```

```{r}
#voom() vrsi transformaciju cpm data u logcpm, procenjuje odnos izmedju srednje vrednosti i 
#varijanse i to koristi da izracuna odgovarajuci opservacioni nivo, tj. prikladne
#observational-level weights, na taj nacin formirajuci podatke da budu spremni za
#linearno modelovanje
voom.y.d <- voom(d, mm, plot = T)
```

# fitting data into lm model 
```{r}
#lmFit() ima za cilja da uklopi podatke u linearni model, tj. da primeni linearni model na svaki
#gen u seriji uzoraka
fit <- lmFit(voom.y.d, mm)

#koeficijenti mogu da se razumeju kao prosecna ekspresija (izmedju) biological replicates; sto 
#je > broj > je i nivo ekspresije, a kako je u pitanju log skala, razlika za 1 predstavlja
#dvostruko po>nje
coef.fit <- fit$coefficients 
head(coef(fit)) #coef() of varijable fit
head(coef.fit) #head() od varijable coef.fit
```

# Chapter 3 - Establish sample group for DEGs analysis 

```{r}
# An example on how to use make contrast (OVO JE SAMO PRIMER KAKO SE KORISTI!)
x <- c("B-A","C-B","C-A") 
makeContrasts(contrasts=x,levels=c("A","B","C"))
```

- Using Contrast and contrast fit
```{r}
#makeContrasts() se koristi samo da se iz vektora sa imenima grupa odabere koje grupe zelis da
#medjusobno uporedis
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef.fit))
contr
```

- toptable(): Extract a table of the top-ranked genes from a linear model fit.    
- eBayes -  Empirical Bayes Statistics for Differential Expression
```{r}
#Za razliku od makeContrasts() gde je cilj da se samo selektuju grupe od interesa (ne porede se
#vrednosti ocitavanja/counts, vec se samo vrsi odabir grupa koje zelis da uporedis u narednom
#koraku), contrasts.fit() ce iz rezultata linearnog modela (fit) da "izvuce" podatke vezane za
#grupe (contr) koje si obelezio da zelis da uporedis (zahvaljujuci makeContrasts()) i
#uporedice ih
tmp <- contrasts.fit(fit, contr)

#Empirijski Bayes metod "pozajmljuje" informacije od svih gena kako bi obezbedio stabilnu 
#procenu varijanse i tipicno ima veze sa modelovanjem veze izmedju ekspresije gena i 
#varijanse gena;
#Multiple testing errors, tj. uklanja lazno pozitivne rezultate iz analize
tmp <- eBayes(tmp) 

#Omogucava ekstrakciju DEG, uz to da ih sortira po p-vrednosti (sort.by = "P")
top.table <- topTable(tmp, sort.by = "P", n = Inf) 
head(top.table, 5) 
```


```{r}
#compared to coef
coef.fit
head(coef.fit)
#logFC iz top.table odgovara (makar bi trebalo da odgovara) razlici izmedju vrednosti za 5. i 2.
#kolonu coef.fit (, jer koeficijenti u sustini predstavljaju prosecnu ekspresiju bioloskih
#replikata za dati gen date grupe), pa ako se uporede trebalo bi da razlika u vrednosti 
#ekspresije izmedju 5. i 2. kolone (u log skali) odgovara vrednosti za logFC za taj gen u toj
#grupi u slucaju top.table
coef_DEG <- coef.fit[rownames(coef.fit) %in% rownames(top.table)[1:5],]
coef_DEG
coef_DEG[,c(5,2)]
```

```{r}
length(which(top.table$adj.P.Val < 0.05))
library(dplyr)
DEGs <- top.table %>%  arrange(logFC) %>% filter(adj.P.Val <0.05) 
head(DEGs)
```
Export gene list to csv
```{r}
top.table$Gene <- rownames(top.table)

#promena redosleda kolona, tako da kolona "Gene" umesto sedma bude prva, a prvih sest zauzmu
#pozicije od druge do sedme
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
head(top.table)
write.table(top.table, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
```

#  Chpater 3.1  - Change Groupings

```{r}
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05)) # number of DE genes
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "I5_v_C_time6.txt", row.names = F, sep = "\t", quote = F)
head(top.table,20)
```

```{r}
p_data <-   top.table %>% filter(adj.P.Val < 0.05) 

p_data %>%  ggplot(aes(x=adj.P.Val,y=logFC)) + 
  geom_text(label=rownames(p_data), size=4.0,alpha=0.7, aes(col=AveExpr)) 
