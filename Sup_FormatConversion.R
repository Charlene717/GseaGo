##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


##### Import setting and Import* #####
## Import setting*
InFOLName_GE <- "Input_TCGA"  # Input Folder Name
FileName_GE <- "Xena_TCGA_LGG_GE"
FileName_Pheno <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

## Import file
GeneExp.df <- read.table(paste0(InFOLName_GE,"/",FileName_GE), header=T, row.names = 1, sep="\t")
Anno.df <- read.table(paste0(InFOLName_GE,"/",FileName_Pheno), header=T, row.names = 1, sep="\t")

##### Export #####
## Export tsv files
GeneExp.df <- data.frame(Gene=row.names(GeneExp.df),GeneExp.df)
write.table(GeneExp.df,paste0(InFOLName_GE,"/",FileName_GE,".tsv"), col.names =T, row.names = F, sep="\t", quote = F)

Anno.df <- data.frame(Gene=row.names(Anno.df),Anno.df)
write.table(Anno.df,paste0(InFOLName_GE,"/",FileName_Pheno,".tsv"), col.names =T, row.names = F, sep="\t", quote = F)

## Create small dataset for test
GeneExpS.df <- GeneExp.df[,sample(ncol(GeneExp.df),100)]
write.table(GeneExpS.df,paste0(InFOLName_GE,"/",FileName_GE,"_S.tsv"), col.names =T, row.names = F, sep="\t", quote = F)



