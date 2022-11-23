# Final Project Outline

## Title
Differential Gene Expression in Liver and Intraheptic Bile Duct Cancers. 

## Author
Riya Shah

## Overview of project
I will identify differentially expressed genes for Liver and Intraheptic Bile Duct Cancers. My analysis will compare the vital status of those who received treatment or therpay and those who did not, while grouping for age - below 50 and above 50. 
My analysis will utilize the DeSEQ2 package and follow this vignette - http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html.

## Data
I will use the data from [https://portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository). I will use the TCGA-LIHC cohort which contains more than 200 cases. For my study, I will randomly select 40 cases and perform my analysis on those. 

## Milestone 1
I will create a dataframe for my analysis by extracting useful data from multiple case and clinical files. This will be completed by November 22, 2022 - end of day. 
For easy analysis, I chose the first 100 files from the multiple case files. I merged all necessary data into 2 files - clinical_data and genes_per_file. The former contains the case id, case submitter id, treatment or therapy, and vital status; the latter is a compilation of case submitter id, gene id, gene name, and unstranded data. Both files can be found in my google drive - https://drive.google.com/drive/folders/1NVYUoYYVc0WLlX9QwJGDQngBfvwq3bdV?usp=sharing. 

## Milestone 2
My entire data will be loaded into vignette (through htseq), for seeking feedback. I will perform a deeper analysis of my data and correct for any known issues. This will happen by November 29, 2022 - end of day.

## Deliverable
My final project which is due on December 3, 2022 - end of day - will include a complete repository with clear documentation and description of my analysis and results.
