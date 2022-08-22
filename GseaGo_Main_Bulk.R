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

##### Extract Target gene and Statistics ####
  # Extract data with TarGene_name
  TarGene_Mean <- GeneExp.df[TarGene_name,] %>%
                  as.numeric() %>%
                  mean()

  #rowMeans(data.matrix(TarGene))
  TarGene_SD <- GeneExp.df[TarGene_name,] %>%
                as.numeric() %>%
                sd()

  # Quartile
  TarGene_Q <- GeneExp.df[TarGene_name,] %>%
               as.numeric() %>%
               quantile()


##### Group the expression matrix according to the expression level of Target gene ####
  if(Mode_Group$Mode=="Mean"){
    if(Mode_Group$SD==0){
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] >= TarGene_Mean+TarGene_SD*(Mode_Group$SD)]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] < TarGene_Mean-TarGene_SD*(Mode_Group$SD)]
    }else{
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] >= TarGene_Mean+TarGene_SD*(Mode_Group$SD)]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] <= TarGene_Mean-TarGene_SD*(Mode_Group$SD)]
    }
    #rm(TarGene_Mean, TarGene_SD)

  }else{
    if(Mode_Group$Q2=="Only"){ # Mode="Quartile"
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] >= TarGene_Q[3]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] < TarGene_Q[3]]
    }else{
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] >= TarGene_Q[4]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGene_name,] <= TarGene_Q[2]]
    }
    #rm(TarGene_Q)

  }


##### Visualization #####
  ## https://www.jianshu.com/p/9e5b7ffcf80f

  data <- reshape2::melt(GeneExp.df[TarGene_name,]%>%
                           as.numeric())
  TGeneDen.p <- ggplot(data,aes(value,fill=value, color=value)) +
    xlab("Expression level") +
    geom_density(alpha = 0.6, fill = "lightgray") +
    geom_rug() + theme_bw()

  ## Set the color
  Mean_SD.clr <- list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" )
  Mean_Q.clr <- list(rect="#abede1", line="#12705f",text="#12705f" )

  ## Plot Mean and SD
  TGeneDen.p
  TGeneDen_SD.p <- ggPlot_vline(TGeneDen.p,data)
  TGeneDen_SD.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= TarGene_name, x ="Expression level", y = "Density") -> TGeneDen_SD.p

## Plot Quartiles
  TGeneDen_Q.p <- ggPlot_vline(TGeneDen.p,data,
                               Line.clr = Mean_Q.clr,
                               Line1 = TarGene_Q[2],
                               Line2 = TarGene_Q[3],
                               Line3 = TarGene_Q[4],
                               Text.set = c("Q1","Q2","Q3"),
                               rectP = list(xWidth=0.015, yminP=0.45, ymaxP=0.55,alpha=0.8)
  )

  TGeneDen_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= TarGene_name, x ="Expression level", y = "Density") -> TGeneDen_Q.p

## Plot Quartiles & Mean and SD
  TGeneDen_SD_Q.p <- ggPlot_vline(TGeneDen_SD.p,data,
                                  Line.clr = Mean_Q.clr,
                                  Line1 = TarGene_Q[2],
                                  Line2 = TarGene_Q[3],
                                  Line3 = TarGene_Q[4],
                                  Text.set = c("Q1","Q2","Q3"),
                                  Text.yPos = 0.35,
                                  rectP = list(xWidth=0.015, yminP=0.3, ymaxP=0.4,alpha=0.8)
  )

  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= TarGene_name, x ="Expression level", y = "Density") -> TGeneDen_SD_Q.p




  pdf(
    file = paste0(Save.Path,"/",SampleName,"_",TarGene_name,"_DensityPlot.pdf"),
    width = 10,  height = 8
  )
  print(TGeneDen_SD.p)
  print(TGeneDen_Q.p)
  print(TGeneDen_SD_Q.p)

  dev.off()

  # Export PPT
  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
                                      OL_Thick = 1.5) +
    labs(title= TarGene_name,
         x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p

  topptx(TGeneDen_SD_Q2.p,paste0(Save.Path,"/",SampleName,"_",TarGene_name,"_DensityPlot.pptx"))

  rm(TGeneDen_SD_Q2.p)


##### Note #####
  ## Finding Peak Values For a Density Distribution
  # http://ianmadd.github.io/pages/PeakDensityDistribution.html
  which.max(density(data$value)$y)
  max(density(data$value)$y)

  ## Plot multiple gene




##### Build files for GSEA official input #####
  source("FUN_GSEA_ForOFFL.R")

  FUN_GSEA_ForOFFL(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
                   TarGeneName = TarGene_name, GroupMode = Mode_Group,
                   Save.Path = Save.Path, SampleName = SampleName)



# #### Save RData ####
#   save.image(paste0(Save.Path,"/GseaGo.RData"))




