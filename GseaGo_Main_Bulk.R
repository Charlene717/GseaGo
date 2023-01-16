##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  source("FUN_Package_InstLoad.R")
  Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")

  BiocManager.set <- c("clusterProfiler","enrichplot","pathview","limma")
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
  # c(organism,"fgsea")

  Package_InstLoad(Basic.set = Basic.set, BiocManager.set = BiocManager.set)

##### Function setting #####
  ## Call function
  source("FUN_DistrPlot.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")
  source("FUN_ggPlot_vline.R")

##### Import setting and Import* #####
  ## File setting*
  InFOLName_GE <- "Input_TCGA"  # Input Folder Name
  SampleName <- "Xena_TCGA_LGG_GE"
  SamplePhenoName <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

  ## Import genetic data file
  GeneExp.df <- read.table(paste0(InFOLName_GE,"/",SampleName), header=T, row.names = 1, sep="\t")
  colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
  GeneExp_Ori.df <- GeneExp.df

  Anno.df <- read.table(paste0(InFOLName_GE,"/",SamplePhenoName), header=T, sep="\t")
  Anno_Ori.df <- Anno.df
  row.names(Anno.df) <- Anno.df[,1]

  ## Reorder the Anno.df
  Anno.df <- left_join(data.frame("sampleID"=colnames(GeneExp.df)),
                       Anno.df)

  ## Import GSEA gene sets
  InFOLName_Genesets <- "Input_Genesets"
  InputGSEA <- "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
  Pathway.all <- read.delim2(paste0(getwd(),"/",InFOLName_Genesets,"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InFOLName_Genesets,"/",InputGSEA))),
                             header = F,sep = "\t")

##### Conditions setting* #####
  SpeciesSet = "Homo sapiens"
  DEGThr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

  Group_Mode <- "GoupByPheno"   # c("GoupByPheno","GoupByGeneExp")

  ## GoupByGeneExp Setting
  TarGene_name <- "TP53"
  GeneExpSet.lt <- list(GeneExpMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                        UpCutoff = 1, LowerCutoff = 1)

  ## GoupByPheno Setting
  GrpCompare_Pheno.lt <- list(Type = "sample_type", GroupPair = c("Primary Tumor","Recurrent Tumor"))

  if(Group_Mode == "GoupByGeneExp"){
    ## Group by GeneExp
    AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("High","Low") )   ## DEG by GeneExp group
  }else{
    ## Group by Pheno
    AnnoSet.lt <- list(GroupType = GrpCompare_Pheno.lt[["Type"]],
                       GroupCompare = GrpCompare_Pheno.lt[["GroupPair"]] )
  }

##### Current path and new folder setting* #####
  ProjectName = "GSEA_TCGA"
  Sampletype = "LGG"

  ExportAnno2 = "Recur2Prim"

  if(Group_Mode == "GoupByGeneExp"){
    if(GeneExpSet.lt$GeneExpMode == "Customize"){
      ExportAnno = paste0(Group_Mode,"_",TarGene_name,"_",GeneExpSet.lt$GeneExpMode,"_Up", GeneExpSet.lt$UpCutoff,
                          "_Low_" ,GeneExpSet.lt$LowerCutoff,"_",ExportAnno2)
    }else{
      ExportAnno = paste0(Group_Mode,"_",TarGene_name,"_",GeneExpSet.lt$GeneExpMode,"_",ExportAnno2)
    }

  }else{
    ExportAnno = paste0(Group_Mode,"_",ExportAnno2)
  }


  ExportName = paste0(ProjectName,"_",Sampletype,"_",ExportAnno)

  Version = paste0(Sys.Date(),"_",ExportName)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}


##### Update the genename ####
  ## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
  library(limma)

  UpdateGene <- "Yes"  # UpdateGene <- c("Yes","No")

  if(UpdateGene == "Yes"){
    source("RUN_FUN_UpdateGeneName.R")
  }

#************************************************************************************************************************#
##### Data preprocess setting #####
  ## Select Pheno column
  colnames(Anno.df)

  PhenoColKeep.set <- c("sampleID","X_PATIENT","histological_type","sample_type","gender")
  Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  colnames(Anno.df)

  head(Anno.df)

  ## Select Pheno row
  PhenoRowKeep.set <- list(col="sample_type",row=c("Primary Tumor","Recurrent Tumor"))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% row.names(Anno.df)]

  # ## Replace
  # Anno.df[,"sample_type"] <- gsub("Primary Tumor", "PrimTu", Anno.df[,"sample_type"])



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
##### Grouping by GeneExp #####
  source("FUN_Group_GE.R")
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                    TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                    Save.Path = Save.Path, ExportName = ExportName)
  Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
  GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
  GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

#************************************************************************************************************************#
##### Run Enrichment analysis in R #####
  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                                  GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                                  ThrSet = DEGThr.lt,
                                  TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID = "sampleID",
                                  Save.Path = Save.Path, ExportName = ExportName, AnnoName = "")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


    # #### Test: DEG by GeneExp group ####
    # AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("High","Low") )
    # source("FUN_DEG_Analysis.R")
    # DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
    #                                 GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
    #                                 ThrSet = DEGThr.lt,
    #                                 TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID = "sampleID",
    #                                 Save.Path = Save.Path, SampleName = SampleName, AnnoName = "AvB")
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

  GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, CMGeneSet = Pathway.all,
                                  DefaultGeneSet = "C2", Species = SpeciesSet, # Speices type can check by msigdbr_species()
                                  NumGenesetsPlt = 15,
                                  TarGeneName = TarGene_name,
                                  ThrSet = DEGThr.lt,
                                  Save.Path = Save.Path, ExportName = ExportName, AnnoName = "Path",
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
  if(Group_Mode == "GoupByGeneExp"){
     Group1.set <- GeneExp_high.set
     Group2.set <- GeneExp_low.set

  }else{
    Group1.set <- Anno.df[Anno.df[,GrpCompare_Pheno.lt[["Type"]] ]%in% GrpCompare_Pheno.lt[["GroupPair"]][1],][,1]
    Group2.set <- Anno.df[Anno.df[,GrpCompare_Pheno.lt[["Type"]] ]%in% GrpCompare_Pheno.lt[["GroupPair"]][2],][,1]

  }


  FUN_GSEA_ForOFFL(GeneExp.df,
                   Group1 = Group1.set, Group2 = Group2.set,
                   GroupMode = Group_Mode,
                   TarGeneName = TarGene_name, GeneExpSet = GeneExpSet.lt,
                   Save.Path = Save.Path, ExportName = ExportName,
                   AnnoName = "GSEA")

##### Build files for Metascape official input #####



#### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo_",ExportName,".RData"))




