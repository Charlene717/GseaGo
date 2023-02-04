##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)


##### Import setting and Import* #####
  ## File setting*
  InFOLName_GE <- "Input_TCGA"  # Input Folder Name
  SampleName <- "Xena_TCGA_LGG_GE"
  SamplePhenoName <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

  ## Import genetic data file
  GeneExp.df <- read.table(paste0(InFOLName_GE,"/",SampleName), header=T, row.names = 1, sep="\t")
  Anno.df <- read.table(paste0(InFOLName_GE,"/",SamplePhenoName), header=T, row.names = 1, sep="\t")

##### Export #####
  ## Export tsv files
  GeneExp.df <- data.frame(Gene=row.names(GeneExp.df),GeneExp.df)
  write.table(GeneExp.df,paste0(InFOLName_GE,"/",SampleName,".tsv"), col.names =T, row.names = F, sep="\t", quote = F)

  Anno.df <- data.frame(Gene=row.names(Anno.df),Anno.df)
  write.table(Anno.df,paste0(InFOLName_GE,"/",SamplePhenoName,".tsv"), col.names =T, row.names = F, sep="\t", quote = F)


  GeneExpS.df <- GeneExp.df[1:1000,]
  write.table(GeneExpS.df,paste0(InFOLName_GE,"/",SampleName,"_S.tsv"), col.names =T, row.names = F, sep="\t", quote = F)



