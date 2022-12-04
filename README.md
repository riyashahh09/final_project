# Final Project Outline

## Title

Differential Gene Expression in Liver and Intraheptic Bile Duct Cancers. 

## Author

Riya Shah

## Overview of project

I have identified diffrentially expressed genes for Liver and Intraheptic Bile Duct Cancers within a specific cohort. My cohort includeed males diagnosed with this cancer and compared those who received radiation treatment v/s those who did not. My controlling factor was age - I chose to analysize males between 0-60 years of age. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. 

## Data

I will use the data from [https://portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository). For this analysis, I used the TCGA-LIHC cohort that had 115 ht-seq counts files for tumors in males aged of 0-60. Out of those, 10 agreed to receive radiation while 105 did not. All individual case files are availabe here  https://drive.google.com/drive/folders/17NrVJrYvQVCE2wTqheMrbA38HoVIPhRE?usp=sharing and the clinical files can be found here https://github.com/riyashahh09/final_project/tree/main/Clinical%20Files. 

## Milestone 1 

I created 2 separate files in occordance to the vignette specifications. 
My first file was called genes_per_file.csv, it can be found here . Using an additional file - key.tsv - that I downloaded from https://portal.gdc.cancer.gov/repository, I I merged the gene_id from each individual case file to its case_id.

I also extracted case_submitter_id, treatment_or_therapy, vital_status, and treatment_type from the clinical.tsv file and called it clinical_processed.csv. These files can be found here https://github.com/riyashahh09/final_project/tree/main/Final%20Files. 

The code I used to combine the data can be accessed here https://drive.google.com/drive/folders/1NVYUoYYVc0WLlX9QwJGDQngBfvwq3bdV?usp=sharing. 

*Known Issues:*

The first issue I faced was while downloading the files. The number of cases in my Clinical Files did not match those in my Case Files, i.e. certain case files were downloaded more than once. Secondly, when I first created my dataframe, I  did not realize that the vignette required a specific layout. So, I had to go back and modify my dataframe in to a format that was readable and processable by the vignette. 

## Milestone 2

**To load datasets into the vignette**

```R
pasCts <- system.file("genes_per_file.csv")
pasAnno <- system.file("clinical_processed.csv")
cts <- as.matrix(read.csv(cts,sep=",",row.names="gene_id"))
coldata <- read.csv(anno, row.names=1)
coldata <- coldata[,c("treatment_or_therapy","vital_status")]
coldata$treatment_or_therapy <- factor(coldata$treatment_or_therapy)
coldata$vital_status <- factor(coldata$vital_stauts)
head(cts, 2)
```

*Known Issues:* 

Loading the dataframe in to the matrix was most difficult. This was becuse the vignette had very strict and specific rules for the file that it would process. While running my vignette, one of the commands was returning a FALSE when it should have returned TRUE. It took me some time to figure out that there was 1 extra case in my genes_per_file.csv than in my clinical_processed.csv. This was because I had selected against "not reported" in my clinical file but forgot to remove that case file from my folder Case Files.  

## Deliverable

**Loading the libraries and our dataset**

Before we can run the vignette on our data, we will have to load our libraries that will help us perform the DESeq2 Analysis and our Dataset. 
To install any package for DESeq2 Analysis, we can use the following command. 
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("name_of_library")
```
Once we have installed all libraries, we can load them using the following command. 
```R 
library("name_of_library")
```
A few important libraries include DESeq2, IHW, vsn, ggplot2, pheatmeap, and RColorBrewer. 

**Creating a count matrix input**

We can make the DESeq2 dataset using the code below. 

```R
coldata[coldata=="Not"] <- "No"
coldata$treatment_or_therapy <- factor(coldata$treatment_or_therapy)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment_or_therapy)
                              
dds
```
**Pre-filtering**

This step reduced the memory size of the dds data object, and increases the speed of testing functions within DESeq2. We perform a pre filtering to keep only those rows containing at least 10 reads total.

```R 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
**Differential expression analysis**

```R
dds <- DESeq(dds)
res <- results(dds)
res
```
<img width="658" alt="Screen Shot 2022-12-02 at 12 11 50 AM" src="https://user-images.githubusercontent.com/112148797/205476995-1f6c5a49-84e9-4c56-9eae-f98227495077.png">

We can select for the condition of interest and build a results table accordingly. 
```R
res <- results(dds, name="treatment_or_therapy_yes_vs_no")
```

**Results**

**MA plot**

```R 
plotMA(res, ylim=c(-10,10)) 
```

<img width="889" alt="plotMA(res, ylim=c(-10,10))" src="https://user-images.githubusercontent.com/112148797/205472927-6f2f9ad8-44d3-4141-ac7f-6a8c96c95e9d.png">

**Alternate shrinkage methods**

```R
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
```

I used the following command to produce a graph which took significantly less time to run using the data shrink method ashr. It has been compared to the normal graph. 

```R
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

<img width="476" alt="Screenshot 2022-12-03 at 10 38 43 PM" src="https://user-images.githubusercontent.com/112148797/205478063-b89eac1a-0c70-444e-bd8b-c64e6723632c.png">


**Plot counts** 

To analyze the counts per reads for a single gene across the groups, we use plot counts.

```R
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment_or_therapy")
```


<img width="243" alt="Screenshot 2022-12-03 at 7 35 42 PM" src="https://user-images.githubusercontent.com/112148797/205473016-ceea02bb-38fd-48e2-8709-b0a7693ba9d7.png">

The returnData argument only returns a data.frame for plotting with ggplot.
<img width="881" alt="Screen Shot 2022-12-02 at 11 46 26 PM" src="https://user-images.githubusercontent.com/112148797/205473081-52cea3f3-9472-40d0-baef-70f41abc9989.png">

**More information on results column** 

```R
mcols(res)$description
```
<img width="873" alt="Screen Shot 2022-12-02 at 12 10 39 AM" src="https://user-images.githubusercontent.com/112148797/205478026-d750accb-6c51-4975-9833-a5e52a3bf022.png">



**Extracting transformed values**

```R 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
```

**Graphs based on extracted values** 

```R 
ntd <- normTransform(dds)
meanSdPlot(assay(ntd)) 
```
<img width="877" alt="meanSdPlot(assay(ntd))" src="https://user-images.githubusercontent.com/112148797/205473141-114d885e-ee37-442e-a968-8aa0ff34e74f.png">

```R 
meanSdPlot(assay(vsd)) 
```
<img width="893" alt="vsd meansd plot" src="https://user-images.githubusercontent.com/112148797/205473146-b75974cc-72e9-4751-b453-451d07792615.png">

**Sample clustering and visualization** 

```R
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("treatment_or_therapy","vital_status")])
```

**Generating heatmap of count matrix**

```R 
BiocManager::install("pheatmap")
library("pheatmap")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
<img width="888" alt="heatmap1" src="https://user-images.githubusercontent.com/112148797/205473256-825847a9-d643-4b35-888a-cb4c2953d653.png">

**Heatmap of sample to sample analysis**

```R
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment_or_therapy, vsd$vital_status, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
<img width="885" alt="Screen Shot 2022-12-03 at 12 08 07 AM" src="https://user-images.githubusercontent.com/112148797/205473343-171c6e9e-52ee-4b2d-8530-1f2cff3aad8b.png">


**Principal component analysis** 

```R 
pcaData <- plotPCA(vsd, intgroup=c("treatment_or_therapy", "vital_status"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```
<img width="879" alt="Screen Shot 2022-12-03 at 12 24 43 AM" src="https://user-images.githubusercontent.com/112148797/205473345-4d703e8d-e4a4-4d06-a68a-a86484eddc17.png">

*Known Issues:* 

The commnd below was taking very long to run. In the time that it was running, my RStudio kept shutting down. Becasue of this, I was not able to use rlog to plot any graphs.  

```R 
rld <- rlog(dds, blind=FALSE)
```
