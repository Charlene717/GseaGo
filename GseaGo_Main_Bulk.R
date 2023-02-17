##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  source("FUN_Package_InstLoad.R")
  PKG_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
  PKG_BiocManager.set <- c("clusterProfiler","enrichplot","pathview","limma")

  FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

##### Function setting #####
  ## Call function
  source("FUN_DistrPlot.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")
  source("FUN_ggPlot_vline.R")

##### SetImport setting* #####
  ## File setting*
  SetImportPath_FOL <- "Input_TCGA"  # Input Folder Name
  SetImport_GE <- "Xena_TCGA_LGG_GE"
  SetImport_Anno <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

  ## SetImport genetic data file
  GeneExp.df <- read.table(paste0(SetImportPath_FOL,"/",SetImport_GE), header=T, row.names = 1, sep="\t")
  colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
  # GeneExp_Ori.df <- GeneExp.df

  Anno.df <- read.table(paste0(SetImportPath_FOL,"/",SetImport_Anno), header=T, sep="\t")
  # Anno_Ori.df <- Anno.df
  row.names(Anno.df) <- Anno.df[,1]

  ## Reorder the Anno.df
  Anno.df <- left_join(data.frame("sampleID"=colnames(GeneExp.df)),
                       Anno.df)

  ## SetImport GSEA gene sets
  SetImportPath_Genesets_FOL <- "Input_Genesets/Gsea_Genesets_Hs"
  SetImport_GSEAGeneSet <- "msigdb.v2022.1.Hs.symbols.gmt"
  SetImport_GSEAGeneSet_MetaData <- "msigdb_v2022.1.Hs.txt"

  GSEAGeneset.df <- read.delim2(paste0(getwd(),"/",SetImportPath_Genesets_FOL,"/",SetImport_GSEAGeneSet),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",SetImportPath_Genesets_FOL,"/",SetImport_GSEAGeneSet))),
                             header = F,sep = "\t")

  GSEAGeneSet_MetaData.df <- read.delim2(paste0(getwd(),"/",SetImportPath_Genesets_FOL,"/",SetImport_GSEAGeneSet_MetaData),sep = "\t")
  GSEAGeneSet_MetaData.df <- GSEAGeneSet_MetaData.df[,c("STANDARD_NAME","SYSTEMATIC_NAME","CATEGORY_CODE","DESCRIPTION_BRIEF","DESCRIPTION_FULL")]

##### Conditions setting* #####
  Set_Species <- "Homo sapiens"
  Set_GroupMode <- "GoupByPheno"  # "GoupByPheno", "GoupByGeneExp"

  TarGene_name = "TP53"
  Set_TarGene <-  list(TarGeneName = TarGene_name,
                       GEGroupMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                       UpCutoff = 1, LowerCutoff = 1)


  if(Set_GroupMode == "GoupByPheno"){
    Set_GroupCond <-  list(GroupType = "sample_type",
                           GroupPair = c("Primary Tumor","Recurrent Tumor"))
  }else if(Set_GroupMode == "GoupByGeneExp"){
    Set_GroupCond <- list(GroupType = TarGene_name, GroupPair = c("High","Low") )   ## DEG by GeneExp group
  }else{
    print("Please set Set_GroupMode by GoupByPheno or GoupByGeneExp")
  }


  ## - [] Add Metric for ranking gene
  ## - [] Modify DGE
  DEGThr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )



##### Current path and new folder setting* #####
  Export_ProjectName = "GSEA_TCGA"
  Export_Sampletype = "LGG"

  Export_Anno = "Recur2Prim"

  if(Set_GroupMode == "GoupByGeneExp"){
    if(Set_TarGene$GEGroupMode == "Customize"){
      Export_Cond = paste0(Set_GroupMode,"_",TarGene_name,"_",Set_TarGene$GEGroupMode,"_Up", Set_TarGene$UpCutoff,
                          "_Low_" ,Set_TarGene$LowerCutoff,"_",Export_Anno)
    }else{
      Export_Cond = paste0(Set_GroupMode,"_",TarGene_name,"_",Set_TarGene$GEGroupMode,"_",Export_Anno)
    }

  }else{
    Export_Cond = paste0(Set_GroupMode,"_",Export_Anno)
  }


  Export_Name = paste0(Export_ProjectName,"_",Export_Sampletype,"_",Export_Cond)

  Save.Path = paste0(getwd(),"/",Sys.Date(),"_",Export_Name)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}

  ## -[] Add setting record

#************************************************************************************************************************#
##### Update the genename ####
  ## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
  Set_UpdateGene <- "Yes"  # Set_UpdateGene <- c("Yes","No")

  if(Set_UpdateGene == "Yes"){

    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
    library(limma)

    source("RUN_UpdateGeneName.R")
  }

#************************************************************************************************************************#
##### Data preprocess* #####
  ## Select Pheno column
  colnames(Anno.df)

  PhenoColKeep.set <- c("sampleID","X_PATIENT","histological_type","sample_type","gender")
  Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  colnames(Anno.df)

  head(Anno.df)

  ## Select Pheno row
  PhenoRowKeep.set <- list(col="sample_type",row=c("Primary Tumor","Recurrent Tumor"))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$sampleID]

  # ## Replace
  # Anno.df[,"sample_type"] <- gsub("Primary Tumor", "PrimTu", Anno.df[,"sample_type"])

  ## -[] Normalization or Standardization


#************************************************************************************************************************#
##### Visualization for Exploratory Data Analysis(EDA) #####
  source("FUN_DistrPlot.R")

  if(Set_GroupMode == "GoupByGeneExp"){
    Plot.DistrPlot <- FUN_DistrPlot_GE(GeneExp.df,
                                       TarGeneName = TarGene_name, GroupSet = Set_TarGene,
                                       Save.Path = Save.Path, ExportName = Export_Name)
    Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
    Plot.DistrPlot_SD_Q

  }else if(Set_GroupMode == "GoupByPheno"){
    p <- ggplot(Anno.df, aes(x=as.factor(Anno.df[,Set_GroupCond$GroupType]), fill=as.factor(Anno.df[,Set_GroupCond$GroupType]))) + geom_bar()
    # p <- ggplot(Anno.df, aes(x=as.factor(gender), fill=as.factor(gender))) + geom_bar()

    p + labs(fill=Set_GroupCond$GroupType, x=Set_GroupCond$GroupType, y = "count")+
        theme_classic() %>% BeautifyggPlot(AxisTitleSize=2,LegPos = c(0.75, 0.85))+
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  }else{
    print("Please set the GroupMode as GoupByPheno or GoupByGeneExp")
  }



#************************************************************************************************************************#
##### Grouping by GeneExp #####
  source("FUN_Group_GE.R")
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                    TarGeneName = TarGene_name, GroupSet = Set_TarGene,
                                    Save.Path = Save.Path, ExportName = Export_Name)
  Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
  GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
  GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

#************************************************************************************************************************#
##### Run Enrichment analysis in R #####
  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                                  GroupType = Set_GroupCond[["GroupType"]], GroupCompare = Set_GroupCond[["GroupPair"]],
                                  ThrSet = DEGThr.lt,
                                  TarGeneName = TarGene_name, GroupMode = Set_TarGene, SampleID = "sampleID",
                                  Save.Path = Save.Path, ExportName = Export_Name, AnnoName = "")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


    # #### Test: DEG by GeneExp group ####
    # Set_GroupCond <- list(GroupType = TarGene_name, GroupPair = c("High","Low") )
    # source("FUN_DEG_Analysis.R")
    # DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
    #                                 GroupType = Set_GroupCond[["GroupType"]], GroupCompare = Set_GroupCond[["GroupPair"]],
    #                                 ThrSet = DEGThr.lt,
    #                                 TarGeneName = TarGene_name, GroupMode = Set_TarGene, SampleID = "sampleID",
    #                                 Save.Path = Save.Path, SetImport_GE = SetImport_GE, AnnoName = "AvB")
    # DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]



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

  GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, CMGeneSet = GSEAGeneset.df,
                                  DefaultGeneSet = "C2", Species = Set_Species, # Speices type can check by msigdbr_species()
                                  NumGenesetsPlt = 15,
                                  TarGeneName = TarGene_name,
                                  ThrSet = DEGThr.lt,
                                  Save.Path = Save.Path, ExportName = Export_Name, AnnoName = "Path",
                                  Keyword = "HALLMARK",
                                  Int_Path =  Int_Path.set,
                                  pAdjustMethod = "BH",  # pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                  nPerm = 100000,
                                  minGSSize = 15, maxGSSize = 500)

  #### Run ORA ####
  ## FUN ORA

#************************************************************************************************************************#
##### Build files for GSEA official input #####
  source("FUN_GSEA_ForOFFL.R")
  if(Set_GroupMode == "GoupByGeneExp"){
     Group1.set <- GeneExp_high.set
     Group2.set <- GeneExp_low.set

  }else{
    Group1.set <- Anno.df[Anno.df[,Set_GroupCond[["GroupType"]] ]%in% Set_GroupCond[["GroupPair"]][1],][,1]
    Group2.set <- Anno.df[Anno.df[,Set_GroupCond[["GroupType"]] ]%in% Set_GroupCond[["GroupPair"]][2],][,1]

  }


  FUN_GSEA_ForOFFL(GeneExp.df,
                   Group1 = Group1.set, Group2 = Group2.set,
                   GroupMode = Set_GroupMode,
                   TarGeneName = TarGene_name, GeneExpSet = Set_TarGene,
                   Save.Path = Save.Path, ExportName = Export_Name,
                   AnnoName = "GSEA")


##### Build files for Metascape official input #####




#### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo_",Export_Name,".RData"))




