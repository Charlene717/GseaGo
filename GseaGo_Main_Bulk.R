##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("SeuratData")) install.packages("SeuratData")
  if(!require("patchwork")) install.packages("patchwork")
  if(!require("eoffice")) install.packages("eoffice")

  library(tidyverse)
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  library(eoffice)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

##### Current path and new folder setting* #####
  ProjectName = "ifnb"
  Sampletype = "PBMC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

  ## Import information
  InputFolder = "Input_files_10x"
  InputAnno = "PBMC_Ano.csv"

  InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"

  ##### Files setting and import * #####
  ## File setting*
  FileName <- "Xena_TCGA_LGG_GE"

  ## Import genetic data file
  GeneExp.df <- read.table(FileName, header=T, row.names = 1, sep="\t")

  ##### Conditions setting* #####
  Target_gene_name <- "TP53"
  Mode_Group <- list(Mode="Mean",SD=1) # Mode_Group <- list(Mode="Quartile",Q2="Only")

  ##### Current path and new folder setting #####
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)



  ##### Extract Target gene and Statistics ####
  # Extract data with Target_gene_name
  Target_gene_Mean <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>%
    mean()

  #rowMeans(data.matrix(Target_gene))
  Target_gene_SD <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>%
    sd()

  # Quartile
  Target_gene_Q <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>%
    quantile()


  ##### Group the expression matrix according to the expression level of Target gene ####
  if(Mode_Group$Mode=="Mean"){
    if(Mode_Group$SD==0){
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
    }else{
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
    }
    #rm(Target_gene_Mean, Target_gene_SD)

  }else{
    if(Mode_Group$Q2=="Only"){ # Mode="Quartile"
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[3]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Q[3]]
    }else{
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[4]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Q[2]]
    }
    #rm(Target_gene_Q)

  }

  ##### Build Expression matrix for GSEA #####
  GeneExp_GSEA.df <- cbind(
    NAME=row.names(GeneExp.df),
    Description = rep("na", nrow(GeneExp.df)),
    GeneExp.df[,c(GeneExp_high.set, GeneExp_low.set)]
  )

  GSEA_SampleCol.df <- data.frame(t(colnames(GeneExp_GSEA.df)), stringsAsFactors=FALSE)
  colnames(GSEA_SampleCol.df) <- GSEA_SampleCol.df

  GeneExp_GSEA.df <- rbind(GSEA_SampleCol.df,GeneExp_GSEA.df)

  GeneExp_GSEA.df <- data.frame(
    "NAME" = c("#1.2",nrow(GeneExp.df)),
    "Description" = c('',length(c(GeneExp_high.set, GeneExp_low.set)))
  ) %>%
    rbind.fill(GeneExp_GSEA.df)
  rm(GSEA_SampleCol.df)

  ##### Build Group Files #####
  ## Set the group array
  Pheno_sum.df <- c(ncol(GeneExp_GSEA.df)-2,2,1) %>% t() %>% data.frame() %>%
    rbind.fill(c(paste0("#",Target_gene_name,"_high"),paste0(Target_gene_name,"_Low")) %>% t() %>% data.frame(stringsAsFactors=FALSE)) %>%
    rbind.fill(c(rep(0,length(GeneExp_high.set)),rep(1,length(GeneExp_low.set))) %>% t() %>% data.frame())


  ##### Export Result #####
  if(Mode_Group$Mode=="Mean"){
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,Mode_Group$SD,"SD_",
                  Target_gene_name,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,Mode_Group$SD,"SD_",
                  Target_gene_name,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )
  }else{
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,"Q2",Mode_Group$Q2,"_",
                  Target_gene_name,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,"Q2",Mode_Group$Q2,"_",
                  Target_gene_name,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )

  }

  ##### Visualization #####
  ## https://www.jianshu.com/p/9e5b7ffcf80f

  data <- reshape2::melt(GeneExp.df[Target_gene_name,]%>%
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
    labs(title= Target_gene_name, x ="Expression level", y = "Density") -> TGeneDen_SD.p

  ## Plot Quartiles
  TGeneDen_Q.p <- ggPlot_vline(TGeneDen.p,data,
                               Line.clr = Mean_Q.clr,
                               Line1 = Target_gene_Q[2],
                               Line2 = Target_gene_Q[3],
                               Line3 = Target_gene_Q[4],
                               Text.set = c("Q1","Q2","Q3"),
                               rectP = list(xWidth=0.015, yminP=0.45, ymaxP=0.55,alpha=0.8)
  )

  TGeneDen_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= Target_gene_name, x ="Expression level", y = "Density") -> TGeneDen_Q.p

  ## Plot Quartiles & Mean and SD
  TGeneDen_SD_Q.p <- ggPlot_vline(TGeneDen_SD.p,data,
                                  Line.clr = Mean_Q.clr,
                                  Line1 = Target_gene_Q[2],
                                  Line2 = Target_gene_Q[3],
                                  Line3 = Target_gene_Q[4],
                                  Text.set = c("Q1","Q2","Q3"),
                                  Text.yPos = 0.35,
                                  rectP = list(xWidth=0.015, yminP=0.3, ymaxP=0.4,alpha=0.8)
  )

  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= Target_gene_name, x ="Expression level", y = "Density") -> TGeneDen_SD_Q.p




  pdf(
    file = paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,"_DensityPlot.pdf"),
    width = 10,  height = 8
  )
  print(TGeneDen_SD.p)
  print(TGeneDen_Q.p)
  print(TGeneDen_SD_Q.p)

  dev.off()

  # Export PPT
  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
                                      OL_Thick = 1.5) +
    labs(title= Target_gene_name,
         x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p

  topptx(TGeneDen_SD_Q2.p,paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,"_DensityPlot.pptx"))

  rm(TGeneDen_SD_Q2.p)


  ##### Note #####
  ## Finding Peak Values For a Density Distribution
  # http://ianmadd.github.io/pages/PeakDensityDistribution.html
  which.max(density(data$value)$y)
  max(density(data$value)$y)

  ## Plot multiple gene

  #### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo.RData"))




