##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  # if(!require("tidyverse")) install.packages("tidyverse")
  # library(tidyverse)

  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  #### BiocManager installation ####
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("clusterProfiler","enrichplot","pathview") # c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  # options(stringsAsFactors = FALSE)


  # Sys.setlocale(category = "LC_ALL", locale = "UTF-8")


##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_ggPlot_vline.R")
  source("FUN_GSEA_ANAL.R")
  source("FUN_DistrPlot.R")

##### Load RData* #####
  # load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-08-13_PBMC_Main/06_Cell_type_annotation.RData")
  load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")

  # Clean up Object
  scRNA.SeuObj <- PBMC.combined
  rm(list=setdiff(ls(), "scRNA.SeuObj"))

  ## Save Ori
  scRNA_Ori.SeuObj <- scRNA.SeuObj

  ## Clean up data (Delete other type)
  scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]

  ## Clean up data (Delete other type)


##### Import setting and Import #####
  ## Import GSEA gene sets
  # InputGSEA <- "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
  # InputGSEA <- "m5_go_bp_v0_3_symbols.gmt"  # InputGSEA <- "m2.all.v0.3.symbols.gmt"

  InputGSEA <- "m5_go_bp_v0_3_symbols.gmt"
  InFOLName_GSEA <- "Input_Genesets"
  Pathway.all <- read.delim2(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA))),
                             header = F,sep = "\t")

##### Conditions setting* #####
  Group_Mode <- "GoupByGeneExp"   # c("GoupByPheno","GoupByGeneExp")
  TarGene_name <- "Chil3"

  GeneExpSet.lt <- list(GeneExpMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                        UpCutoff = 1, LowerCutoff = 1)

  if(Group_Mode == "GoupByGeneExp"){
    ## Group by GeneExp
    AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("High","Low") )   ## DEG by GeneExp group
  }else{
    ## Group by Pheno
    AnnoSet.lt <- list(GroupType = "sample_type", GroupCompare = c("Primary Tumor","Recurrent Tumor") )
  }

  Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

##### Current path and new folder setting* #####
  ProjectName = "CC10X"
  Sampletype = "PBMC"

  ExportAnno2 = "EO_Neu"
  if(Group_Mode == "GoupByGeneExp"){
    ExportAnno = paste0(TarGene_name,"_",GeneExpSet.lt$GeneExpMode,"_",ExportAnno2)

  }else{
    ExportAnno = paste0(Group_Mode,"_",ExportAnno2)

  }

  # ExportAnno = "Chil3Mean_PathM2"
  # ExportAnno = "Recur2Prim"

  ExportName = paste0(ProjectName,"_",Sampletype,"_",ExportAnno)


  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype,"_", ExportAnno)
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


  ## Update the genename ##* Take very long time
  UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
  if(UpdateGene == "Yes"){
    row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
  }

#************************************************************************************************************************#
##### Data preprocess setting #####
  ## Extract data from scRNA.SeuObj
  GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  Anno.df <- scRNA.SeuObj@meta.data
  Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)
  row.names(Anno.df) <- Anno.df[,1]

  ## Select Pheno column
  Anno_Ori.df <- Anno.df
  colnames(Anno.df)

  # PhenoColKeep.set <- c("X_INTEGRATION","X_PATIENT","histological_type","sample_type","gender")
  # Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  # colnames(Anno.df)
  #
  # head(Anno.df)

  ## Select Pheno row
  PhenoRowKeep.set <- list(col="Cachexia" ,row=c("EO"))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]

  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]

  ## Select Pheno row2
  PhenoRowKeep.set <- list(col="celltype" ,row=c("Neu"))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]

  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]


  # ## Delete specific cell type
  # ## Clean up data
  # Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
  # GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]



#************************************************************************************************************************#
##### Visualization #####
  source("FUN_DistrPlot.R")
  ##### Group by gene expression 1: CutOff by total  #####
  Plot.DistrPlot <- FUN_DistrPlot(GeneExp.df,
                                  TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                  Save.Path = Save.Path, ExportName = ExportName)
  Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
  Plot.DistrPlot_SD_Q


#************************************************************************************************************************#
##### Grouping #####
  source("FUN_Group_GE.R")
  ##### Group by gene expression 1: CutOff by total  #####
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                    TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                    Save.Path = Save.Path, ExportName = ExportName)
  Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
  GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
  GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

  ##### Group by gene expression 2: CutOff by Comparison #####
  ## FUN Comparison (Visualization and value)

  ##### Group by phenotype #####


#************************************************************************************************************************#
##### Run Enrichment analysis in R #####
  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                              GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                              ThrSet = Thr.lt,
                              TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID = "ID",
                              Save.Path = Save.Path, ExportName = ExportName, AnnoName = "AvB")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


  #### Run GSEA ####
  source("FUN_GSEA_ANAL.R")

  GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, CMGeneSet = Pathway.all,
                                  NumGenesetsPlt=15,
                                  TarGeneName = TarGene_name,
                                  ThrSet = Thr.lt, Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                                  Save.Path = Save.Path, ExportName = ExportName, AnnoName = "Path")

  #### Run ORA ####
  ## FUN ORA

#************************************************************************************************************************#
##### Build files for GSEA official input #####
  source("FUN_GSEA_ForOFFL.R")

  FUN_GSEA_ForOFFL(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
                   GroupMode = Group_Mode,
                   TarGeneName = TarGene_name, GeneExpSet = GeneExpSet.lt,
                   Save.Path = Save.Path, SampleName = ExportName,
                   AnnoName="Recur2Prim")

##### Build files for Metascape official input #####



#### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo_",ExportName,".RData"))




