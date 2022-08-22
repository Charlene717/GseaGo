##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("SeuratData")) install.packages("SeuratData")
  if(!require("patchwork")) install.packages("patchwork")
  if(!require("plyr")) install.packages("plyr")
  if(!require("eoffice")) install.packages("eoffice")

  library(tidyverse)
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  library(plyr)
  library(eoffice)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

##### Import setting and Import* #####
  ## File setting*
  InFOLName_GE <- "Input_TCGA"  # Input Folder Name
  SampleName <- "Xena_TCGA_LGG_GE"

  ## Import genetic data file
  GeneExp.df <- read.table(paste0(InFOLName_GE,"/",SampleName), header=T, row.names = 1, sep="\t")

  ## Import GSEA gene sets
  InputGSEA <- "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
  InFOLName_GSEA <- "Input_Genesets"
  Pathway.all <- read.delim2(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InputGSEA))),
                             header = F,sep = "\t")

##### Conditions setting* #####
  ## Group by gene expression
  TarGene_name <- "TP53"
  Mode_Group <- list(Mode="Mean",SD=1)
  # Mode_Group <- list(Mode=c("Mean","Quartile","Tumor2Normal"),Q2="Only")

  ## Group by phenotype


##### Current path and new folder setting* #####
  ProjectName = "TCGA"
  Sampletype = "LGG"

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype,"_",TarGene_name)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Update the genename ####
  # ## Test
  # UpdateSymbolList("SEPT1")
  # A <- UpdateSymbolList(row.names(GeneExp.df))
  # B <- row.names(GeneExp.df)
  # # sum(c("a","c")==c("a","b"))
  # sum(A==B)
  # summary(A==B)
  #
  ## Error: Timeout was reached: [rest.genenames.org] Operation timed out after 10005 milliseconds with 0 bytes received


  ## Update the genename
  UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
  if(UpdateGene == "Yes"){
    row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
  }



#************************************************************************************************************************#
##### Grouping #####

  ##### Group by gene expression 1: CutOff by total  #####


  ##### Group by gene expression 2: CutOff by Comparison #####


  ##### Group by phenotype #####


#************************************************************************************************************************#
##### Run Enrichment analysis in R #####



#************************************************************************************************************************#
##### Build files for GSEA official input #####
  source("FUN_GSEA_ForOFFL.R")

  FUN_GSEA_ForOFFL(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
                   TarGeneName = TarGene_name, GroupMode = Mode_Group,
                   Save.Path = Save.Path, SampleName = SampleName)

##### Build files for Metascape official input #####



#### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo_",SampleName,".RData"))




