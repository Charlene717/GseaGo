##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

  Rec_Time_Point.lt <- list()
  Rec_Time_Spend.lt <- list()

  Rec_Time_Point.lt[["Start_Time"]] <- Sys.time() # %>% as.character()

##### Load Packages #####
  source("FUN_Package_InstLoad.R")
  PKG_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
  PKG_BiocManager.set <- c("clusterProfiler","enrichplot","pathview","limma")

  FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)
  rm(PKG_Basic.set,PKG_BiocManager.set,FUN_Package_InstLoad)

##### Function setting #####
  ## Call function
  source("FUN_DistrPlot_GE.R")
  source("FUN_Beautify_ggplot.R")
  # source("FUN_Find_Markers.R")
  # source("FUN_VolcanoPlot.R")
  # source("FUN_GSEA_LargeGeneSet.R")
  # source("FUN_GSEA_ggplot.R")
  source("FUN_ggPlot_vline.R")
  source("FUN_GSEA_ANAL.R")

##### Import setting and data loading* #####
  Rec_Time_Point.lt[["Input_Start_Time"]] <- Sys.time() # %>% as.character()

  #### (Required) Set input data ####
  ## Set Import path
  SetInputPath_FOL <- "Input_TCGA"  # Input Folder Name
  SetInput_GE <- "Xena_TCGA_LGG_GE"
  SetInput_Meta <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

  ## Load Gene expression file
  GeneExp.df <- read.table(paste0(SetInputPath_FOL,"/",SetInput_GE), header=T, row.names = 1, sep="\t")
  colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
  ## Load Annotation file
  Metadata.df <- read.table(paste0(SetInputPath_FOL,"/",SetInput_Meta), header=T, sep="\t")
  row.names(Metadata.df) <- Metadata.df[,1]

  ## Reorder the Metadata.df
  Metadata.df <- left_join(data.frame("sampleID"=colnames(GeneExp.df)),
                           Metadata.df)
  row.names(Metadata.df) <- Metadata.df[,1]

  #### (Optional) Set GSEA genesets ####
  ## Set Import GSEA genesets path
  SetInputPath_Genesets_FOL <- "Input_Genesets/Gsea_Genesets_Hs"
  SetInput_GSEAGeneSet <- "msigdb.v2022.1.Hs.symbols.gmt"
  SetInput_GSEAGeneSet_MetaData <- "msigdb_v2022.1.Hs.txt"

  ## Load GSEA genesets file
  GSEAGeneset.df <- read.delim2(paste0(getwd(),"/",SetInputPath_Genesets_FOL,"/",SetInput_GSEAGeneSet),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",SetInputPath_Genesets_FOL,"/",SetInput_GSEAGeneSet))),
                             header = F,sep = "\t")

  GSEAGeneSet_MetaData.df <- read.delim2(paste0(getwd(),"/",SetInputPath_Genesets_FOL,"/",SetInput_GSEAGeneSet_MetaData),sep = "\t")
  GSEAGeneSet_MetaData.df <- GSEAGeneSet_MetaData.df[,c("STANDARD_NAME","SYSTEMATIC_NAME","CATEGORY_CODE","DESCRIPTION_BRIEF","DESCRIPTION_FULL")]

Rec_Time_Point.lt[["Input_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Input"]] <- Rec_Time_Point.lt[["Input_End_Time"]] - Rec_Time_Point.lt[["Input_Start_Time"]]


##### Conditions setting* #####
  Set_Species <- "Homo sapiens"
  Set_GroupMode <- "GoupByPheno"  # "GoupByPheno", "GoupByGeneExp"

  if(Set_GroupMode == "GoupByPheno"){
    Set_GroupCond <-  list(GroupType = "sample_type",
                           GroupPair = c("Primary Tumor","Recurrent Tumor"))
  }else if(Set_GroupMode == "GoupByGeneExp"){
    Set_TarGene_name = "TP53"
    Set_TarGene <-  list(TarGeneName = Set_TarGene_name,
                         GEGroupMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                         UpCutoff = 1, LowerCutoff = 1)

    Set_GroupCond <- list(GroupType = Set_TarGene_name, GroupPair = c("High","Low") )   ## DEG by GeneExp group
  }else{
    print("Please set Set_GroupMode by GoupByPheno or GoupByGeneExp")
  }

  ## Set DEG Analysis
  Set_DEGThr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )
  ## - [] Add Metric for ranking gene
  ## - [] Modify DGE


##### Current path and new folder setting* #####
  SetExport_ProjectName = "TCGA"
  SetExport_Sampletype = "LGG"
  SetExport_Anno = "_Recur2Prim" # SetExport_Anno = "_Name"

  if(Set_GroupMode == "GoupByGeneExp"){
    if(Set_TarGene$GEGroupMode == "Customize"){
      SetExport_Cond = paste0(Set_GroupMode,"_",Set_TarGene_name,"_",Set_TarGene$GEGroupMode,"_Up", Set_TarGene$UpCutoff,
                          "_Low_" ,Set_TarGene$LowerCutoff,SetExport_Anno)
    }else{
      SetExport_Cond = paste0(Set_GroupMode,"_",Set_TarGene_name,"_",Set_TarGene$GEGroupMode,SetExport_Anno)
    }

  }else{
    SetExport_Cond = paste0(Set_GroupMode,SetExport_Anno)
  }


  SetExport_Name = paste0(SetExport_ProjectName,"_",SetExport_Sampletype,"_",SetExport_Cond)
  # rm(SetExport_ProjectName,SetExport_Sampletype,SetExport_Anno, SetExport_Cond)

  Save_Path = paste0(getwd(),"/",Sys.Date(),"_",SetExport_Name)
  ## Create new folder
  if (!dir.exists(Save_Path)){dir.create(Save_Path)}

  ## -[] Add setting record

#************************************************************************************************************************#
##### Update the genename ####
Rec_Time_Point.lt[["Update_Genename_Start_Time"]] <- Sys.time() # %>% as.character()

  ## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
  Set_UpdateGene <- "Yes"  # Set_UpdateGene <- c("Yes","No")

  if(Set_UpdateGene == "Yes"){

    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
    library(limma)

    source("RUN_UpdateGeneName.R")
  }

Rec_Time_Point.lt[["Update_Genename_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Update_Genename"]] <- Rec_Time_Point.lt[["Update_Genename_End_Time"]] - Rec_Time_Point.lt[["Update_Genename_Start_Time"]]

#************************************************************************************************************************#
##### Data preprocess* #####
  ## Select Pheno column
  # colnames(Metadata.df)
  PhenoColKeep.set <- c("sampleID","X_PATIENT","histological_type","sample_type","gender")
  Metadata.df <- Metadata.df[,c(PhenoColKeep.set)]
  colnames(Metadata.df)

  # head(Metadata.df)

  ## Select Pheno row
  PhenoRowKeep.set <- list(col="sample_type",row=c("Primary Tumor","Recurrent Tumor"))
  Metadata.df <- Metadata.df[Metadata.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Metadata.df$sampleID]

  rm(PhenoColKeep.set,PhenoRowKeep.set)

  # ## Replace
  # Metadata.df[,"sample_type"] <- gsub("Primary Tumor", "PrimTu", Metadata.df[,"sample_type"])

  ## -[] Normalization or Standardization


#************************************************************************************************************************#
##### Visualization for Exploratory Data Analysis(EDA) #####
Rec_Time_Point.lt[["EDA_Start_Time"]] <- Sys.time() # %>% as.character()

  source("FUN_DistrPlot_GE.R")

  if(Set_GroupMode == "GoupByGeneExp"){
    Plot.DistrPlot <- FUN_DistrPlot_GE(GeneExp.df,
                                       TarGeneName = Set_TarGene_name, GroupSet = Set_TarGene,
                                       Save.Path = Save_Path, ExportName = SetExport_Name)
    Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
    Plot.DistrPlot_SD_Q

  }else if(Set_GroupMode == "GoupByPheno"){
    Plot.Barplot <- ggplot(Metadata.df, aes(x=as.factor(Metadata.df[,Set_GroupCond$GroupType]), fill=as.factor(Metadata.df[,Set_GroupCond$GroupType]))) + geom_bar()
    # Plot.Barplot <- ggplot(Metadata.df, aes(x=as.factor(gender), fill=as.factor(gender))) + geom_bar()

    Plot.Barplot + labs(fill=Set_GroupCond$GroupType, x=Set_GroupCond$GroupType, y = "count")+
        theme_classic() %>% FUN_BeautifyggPlot(AxisTitleSize=2,LegPos = c(0.75, 0.85))+
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) -> Plot.Barplot1
    Plot.Barplot1

    pdf(file = paste0(Save_Path,"/BarPlot_",SetExport_Name,".pdf"),
        width = 10,  height = 8)
    Plot.Barplot1 %>% print()

    dev.off()
    rm(Plot.Barplot, Plot.Barplot1)

  }else{
    print("Please set the GroupMode as GoupByPheno or GoupByGeneExp")
  }

Rec_Time_Point.lt[["EDA_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["EDA"]] <- Rec_Time_Point.lt[["EDA_End_Time"]] - Rec_Time_Point.lt[["EDA_Start_Time"]]


#************************************************************************************************************************#
##### (Optional) Grouping by GeneExp #####
  if(Set_GroupMode == "GoupByGeneExp"){
    Rec_Time_Point.lt[["GrpGeneExp_Start_Time"]] <- Sys.time() # %>% as.character()

    source("FUN_Group_GE.R")
    GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Metadata.df,
                                      TarGeneName = Set_TarGene_name, GroupSet = Set_TarGene,
                                      Save.Path = Save_Path, ExportName = SetExport_Name)
    Metadata.df <- GeneExp_group.set[["AnnoNew.df"]]
    GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
    GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

    Rec_Time_Point.lt[["GrpGeneExp_End_Time"]] <- Sys.time() # %>% as.character()
    Rec_Time_Spend.lt[["GrpGeneExp"]] <- Rec_Time_Point.lt[["GrpGeneExp_End_Time"]] - Rec_Time_Point.lt[["GrpGeneExp_Start_Time"]]

  }

#************************************************************************************************************************#
##### Run Differential Expression Gene(DEG) analysis in R #####
Rec_Time_Point.lt[["DEG_Start_Time"]] <- Sys.time() # %>% as.character()

  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Metadata.df,
                                  GroupType = Set_GroupCond[["GroupType"]], GroupCompare = Set_GroupCond[["GroupPair"]],
                                  ThrSet = Set_DEGThr.lt,
                                  SampleID = "sampleID",
                                  Save.Path = Save_Path, ExportName = SetExport_Name, AnnoName = "")
  DEG_Extract.df <- DEG_ANAL.lt[["DEG_Extract.df"]]

Rec_Time_Point.lt[["DEG_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["DEG"]] <- Rec_Time_Point.lt[["DEG_End_Time"]] - Rec_Time_Point.lt[["DEG_Start_Time"]]

##### Run Enrichment analysis in R #####
Rec_Time_Point.lt[["GSEA_Start_Time"]] <- Sys.time() # %>% as.character()
  #### Run GSEA ####
  source("FUN_GSEA_ANAL.R")
  Int_Path.set <- c(
                    "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
                    #"REACTOME_SEMA4D_MEDIATED_INHIBITION_OF_CELL_ATTACHMENT_AND_MIGRATION",
                    "REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY",
                    "HALLMARK_E2F_TARGETS",
                    "REACTOME_DNA_REPLICATION",
                    "REACTOME_G2_M_CHECKPOINTS",
                    "KEGG_OLFACTORY_TRANSDUCTION",
                    "KEGG_CALCIUM_SIGNALING_PATHWAY"
  )

  GSEA_Result.lt <- FUN_GSEA_ANAL(DEG_Extract.df, CMGeneSet = GSEAGeneset.df,
                                  DefaultGeneSet = "C2", Species = Set_Species, # Speices type can check by msigdbr_species()
                                  NumGenesetsPlt = 15,
                                  TarGeneName = Set_TarGene_name,
                                  ThrSet = Set_DEGThr.lt,
                                  Save.Path = Save_Path, ExportName = SetExport_Name, AnnoName = "Path",
                                  Keyword = "HALLMARK",
                                  Int_Path =  Int_Path.set,
                                  pAdjustMethod = "BH",  # pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                  nPerm = 100000,
                                  minGSSize = 15, maxGSSize = 500)

Rec_Time_Point.lt[["GSEA_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["GSEA"]] <- Rec_Time_Point.lt[["GSEA_End_Time"]] - Rec_Time_Point.lt[["GSEA_Start_Time"]]


  #### Run ORA ####
Rec_Time_Point.lt[["ORA_Start_Time"]] <- Sys.time() # %>% as.character()
  ## FUN ORA
  source("RUN_ORA.R")
Rec_Time_Point.lt[["ORA_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["ORA"]] <- Rec_Time_Point.lt[["ORA_End_Time"]] - Rec_Time_Point.lt[["ORA_Start_Time"]]

#************************************************************************************************************************#
##### Build files for official input #####
#### Build files for GSEA official input ####
Rec_Time_Point.lt[["OFFL_GSEA_Start_Time"]] <- Sys.time() # %>% as.character()

  source("FUN_GSEA_ForOFFL.R")
  if(Set_GroupMode == "GoupByGeneExp"){
     Group1.set <- GeneExp_high.set
     Group2.set <- GeneExp_low.set
     Group1_Name <- paste0(Set_GroupCond$GroupType,"_high")
     Group2_Name <- paste0(Set_GroupCond$GroupType,"_low")


  }else if(Set_GroupMode == "GoupByPheno"){
    Group1.set <- Metadata.df[Metadata.df[,Set_GroupCond[["GroupType"]] ]%in% Set_GroupCond[["GroupPair"]][1],][,1]
    Group2.set <- Metadata.df[Metadata.df[,Set_GroupCond[["GroupType"]] ]%in% Set_GroupCond[["GroupPair"]][2],][,1]
    Group1_Name <- paste0(Set_GroupCond$GroupPair[1]) %>% gsub(" ", "_", .)
    Group2_Name <- paste0(Set_GroupCond$GroupPair[2]) %>% gsub(" ", "_", .)


  }else{
    print("Please Check Set_GroupMode which should be GoupByPheno or GoupByGeneExp")
  }


  FUN_GSEA_ForOFFL(GeneExp.df,
                   Group1 = Group1.set, Group2 = Group2.set,
                   Group1Name = Group1_Name,Group2Name = Group2_Name,
                   SavePath = Save_Path, ExportName = SetExport_Name,
                   AnnoName = "") # AnnoName = "_Name"

Rec_Time_Point.lt[["OFFL_GSEA_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["OFFL_GSEA"]] <- Rec_Time_Point.lt[["OFFL_GSEA_End_Time"]] - Rec_Time_Point.lt[["OFFL_GSEA_Start_Time"]]

#### Build files for Metascape/ORA official input ####




#### Save RData ####
  Rec_Time_Point.lt[["END_Time"]] <- Sys.time() # %>% as.character()
  Rec_time_diff <- Rec_Time_Point.lt[["END_Time"]] - Rec_Time_Point.lt[["Start_Time"]]
  Rec_time_diff
  save.image(paste0(Save_Path,"/GseaGo_",SetExport_Name,".RData"))

#### Record ####
  ## Record time log
  Rec_Time_Point.lt[["END_Time2"]] <- Sys.time()

  Rec_SaveTime_diff <- Rec_Time_Point.lt[["END_Time2"]] - Rec_Time_Point.lt[["END_Time"]]
  Rec_SaveTime_diff

  write(paste(" Program running time：", as.character(Rec_time_diff), "mins\n",
              "Save RData time：", as.character(Rec_SaveTime_diff), "mins"), file = paste0(Save_Path,"/Rec_time_log.txt"))

  ## Record parameter

