# Final Project Outline

## Title
Differential Gene Expression in Liver and Intraheptic Bile Duct Cancers. 

## Author
Riya Shah

## Overview of Project
I have identified diffrentially expressed genes for Liver and Intraheptic Bile Duct Cancers within a specific cohort. My cohort includeed males diagnosed with this cancer and compared those who received radiation treatment v/s those who did not. My controlling factor was age - I chose to analysize males between 0-60 years of age. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. 

## Data
I will use the data from [https://portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository). For this analysis, I used the TCGA-LIHC cohort that had 115 ht-seq counts files for tumors in males aged of 0-60. Out of those, 10 agreed to receive radiation while 105 did not. All individual case files are availabe here  https://drive.google.com/drive/folders/17NrVJrYvQVCE2wTqheMrbA38HoVIPhRE?usp=sharing and the clinical files can be found here https://github.com/riyashahh09/final_project/tree/main/Clinical%20Files. 

## Milestone 1 
I created 2 separate files in occordance to the vignette specifications. 
My first file was called genes_per_file.csv, it can be found here . Using an additional file - key.tsv - that I downloaded from https://portal.gdc.cancer.gov/repository, I I merged the gene_id from each individual case file to its case_id.
I also extracted case_submitter_id, treatment_or_therapy, vital_status, and treatment_type from the clinical.tsv file and called it clinical_processed.csv. These files can be found here https://github.com/riyashahh09/final_project/tree/main/Final%20Files. 

Known Issues: The first issue I faced was while downloading the files. The number of cases in my Clinical Files did not match those in my Case Files, i.e. certain case files were downloaded more than once. Secondly, when I first created my dataframe, I  did not realize that the vignette required a specific layout. So, I had to go back and modify my dataframe in to a format that was readable and processable by the vignette. 

## Milestone 2
I started by loading my new dataframe into the vignette using count matrix input. 
Then I ran a first analysis of my previosuly loaded dataset. I was able to generate a few graphs that I have attached below. 
![WhatsApp Image 2022-11-30 at 2 30 30 AM](https://user-images.githubusercontent.com/112148797/204779937-3fab5524-37f7-4b9b-a91d-47e0fc47f062.jpeg)

Known Issues: Loading the dataframe in to the matrix was most difficult. This was becuse the vignette had very strict and specific rules for the file that it would process. While running my vignette, one of the commands was returning a FALSE when it should have returned TRUE. It took me some time to figure out that there was 1 extra case in my genes_per_file.csv than in my clinical_processed.csv. This was because I had selected against "not reported" in my clinical file but forgot to remove that case file from my folder Case Files.  

## Deliverable
Loading the Libraries and our Dataset:
Before we can run the vignette on our data, we will have to load our libraries that will help us perform the DESeq2 Analysis and our Dataset. 
To install any package for DESeq2 Analysis, we can use the following command. 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("name_of_library")

Once we have installed all libraries, we can load them using the following command. 
library(name_of_library)

A few important libraries include DESeq2, IHW, vsn, ggplot2, pheatmeap, and RColorBrewer. 

**Results**
**MA Plot 

~~~ R plotMA(res, ylim=c(-10,10)) ~~~
<img width="889" alt="plotMA(res, ylim=c(-10,10))" src="https://user-images.githubusercontent.com/112148797/205472927-6f2f9ad8-44d3-4141-ac7f-6a8c96c95e9d.png">

**Plot Counts** 
To analyze the counts per reads for a single gene across the groups, we use plot counts.
~~~ R plotCounts(dds, gene=which.min(res$padj), intgroup="treatment_or_therapy") ~~~


<img width="243" alt="Screenshot 2022-12-03 at 7 35 42 PM" src="https://user-images.githubusercontent.com/112148797/205473016-ceea02bb-38fd-48e2-8709-b0a7693ba9d7.png">

The returnData argument only returns a data.frame for plotting with ggplot.
<img width="881" alt="Screen Shot 2022-12-02 at 11 46 26 PM" src="https://user-images.githubusercontent.com/112148797/205473081-52cea3f3-9472-40d0-baef-70f41abc9989.png">

