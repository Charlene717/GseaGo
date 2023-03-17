#### Differential Expression Gene Analysis ####
## Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
## Ref: https://www.jianshu.com/p/b6912d318de5
## Ref: https://www.jianshu.com/p/0e1ad0cc4ce6

FUN_DEG_Analysis = function(GeneExp.df, Metadata.df,
                            GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                            ThrSet = DEGThr.lt,
                            SampleID = "sampleID",
                            Save.Path = Save_Path, ExportName = ExportName, AnnoName = "AvB"
){

  ##### Load Packages  #####
  source("FUN_Package_InstLoad.R")
  Basic.set <- c("tidyverse","ggplot2")
  BiocManager.set <- c("edgeR","baySeq")
  FUN_Package_InstLoad(Basic.set = Basic.set, BiocManager.set = BiocManager.set)

  #************************************************************************************************************************#


  #### Differential Expression Gene Analysis ####
  library(edgeR)
  Anno_Ints.df <- Metadata.df[Metadata.df[,GroupType] %in% GroupCompare,]
  # Anno_Ints.df <- Metadata.df[Metadata.df$ReCluster %in% c("AD","AC"),]  # c("CoreCD00","CDOri")
  matrix_Ints.df <- GeneExp.df
  # colnames(matrix_Ints.df) <-  gsub("\\.", "-", colnames(matrix_Ints.df))
  matrix_Ints.df <- matrix_Ints.df[,colnames(matrix_Ints.df) %in% Anno_Ints.df[,SampleID]]

  matrix_Ints_ID.df <- data.frame(ID = colnames(matrix_Ints.df))
  colnames(matrix_Ints_ID.df) <- SampleID

  Anno_Ints.df <- left_join(matrix_Ints_ID.df, Anno_Ints.df)

  DGE_Ints.lt <- DGEList(counts=matrix_Ints.df, group=Anno_Ints.df[,GroupType], lib.size=rep(1000,ncol(matrix_Ints.df)))
  DEG_Extract.lt <- exactTest(DGE_Ints.lt,pair = GroupCompare, dispersion=0.2)
  DEG_Extract.df <- topTags(DEG_Extract.lt,n = nrow(matrix_Ints.df)) %>%
                           as.data.frame() %>%
                           data.frame(Gene=row.names(.),.)

  #### Add gene sets by threshold filtering ####
  # ThrSet = Thr.lt
  length(ThrSet)
  DEG_Extract_Flt.df <- DEG_Extract.df

  DEG_Extract_FltH.df <- filter(DEG_Extract.df, DEG_Extract.df[,ThrSet[["LogFC"]][1]] >= ThrSet[["LogFC"]][2] &
                               DEG_Extract.df[,ThrSet[["pVal"]][1]] < ThrSet[["pVal"]][2])

  DEG_Extract_FltL.df <- filter(DEG_Extract.df, DEG_Extract.df[,ThrSet[["LogFC"]][1]] <= as.numeric(ThrSet[["LogFC"]][2])*(-1) &
                               DEG_Extract.df[,ThrSet[["pVal"]][1]] < ThrSet[["pVal"]][2])
  DEG_Extract_Flt.df <- rbind(DEG_Extract_FltH.df,DEG_Extract_FltL.df)

  DEG_Extract_Flt.set <- rownames(DEG_Extract_Flt.df)
  DEG_Extract_FltH.set <- rownames(DEG_Extract_FltH.df)
  DEG_Extract_FltL.set <- rownames(DEG_Extract_FltL.df)


  #### Export file ####
  write.table(DEG_Extract.df, file = paste0(Save.Path,"/DEGAnalysis_",ExportName,"_",AnnoName,".tsv"),
              sep="\t", row.names= F, quote = FALSE)
  write.table(DEG_Extract_Flt.df, file = paste0(Save.Path,"/DEGAnalysis_Flt_",ExportName,"_",AnnoName,".tsv"),
              sep="\t", row.names= F, quote = FALSE)


  #### Output ####
  Output <- list()
  Output[["DEG_Extract.df"]] <- DEG_Extract.df
  Output[["DEG_Extract_Flt.df"]] <- DEG_Extract_Flt.df
  Output[["DEG_Extract_FltH.set"]] <- DEG_Extract_FltH.set
  Output[["DEG_Extract_FltL.set"]] <- DEG_Extract_FltL.set
  Output[["Thr.lt"]] <- ThrSet

  return(Output)

}
