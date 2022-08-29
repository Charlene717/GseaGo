#### Differential Expression Gene Analysis ####
## Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
## Ref: https://www.jianshu.com/p/b6912d318de5

FUN_DEG_Analysis = function(GeneExp.df, Anno.df,
                            GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                            ThrSet = Thr.lt,
                            TarGeneName = TarGene_name, GroupMode = Mode_Group, SampleID = "X_INTEGRATION",
                            Save.Path = Save.Path, SampleName = SampleName, AnnoName = "AvB"
){

  ##### Parameter setting* #####
  # Set the desired organism
  organism = "org.Dm.eg.db"


  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  # BiocManager::install()
  Package.set <- c(organism,"edgeR","baySeq")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### Differential Expression Gene Analysis ####
  library(edgeR)
  Anno_Ints.df <- Anno.df[Anno.df[,GroupType] %in% GroupCompare,]
  # Anno_Ints.df <- Anno.df[Anno.df$ReCluster %in% c("AD","AC"),]  # c("CoreCD00","CDOri")
  matrix_Ints.df <- GeneExp.df
  colnames(matrix_Ints.df) <-  gsub("\\.", "-", colnames(matrix_Ints.df))
  matrix_Ints.df <- matrix_Ints.df[,colnames(matrix_Ints.df) %in% Anno_Ints.df[,SampleID]]

  matrix_Ints_ID.df <- data.frame(ID = colnames(matrix_Ints.df))
  colnames(matrix_Ints_ID.df) <- SampleID

  Anno_Ints.df <- left_join(matrix_Ints_ID.df, Anno_Ints.df)

  DGE_Ints.lt <- DGEList(counts=matrix_Ints.df, group=Anno_Ints.df[,GroupType], lib.size=rep(1000,ncol(matrix_Ints.df)))
  DE_Extract.lt <- exactTest(DGE_Ints.lt, dispersion=0.2)
  DE_Extract.df <- topTags(DE_Extract.lt,n = nrow(matrix_Ints.df)) %>%
                           as.data.frame() %>%
                           data.frame(Gene=row.names(.),.)

  #### Add gene sets by threshold filtering ####
  # ThrSet = Thr.lt
  length(ThrSet)
  DE_Extract_Flt.df <- DE_Extract.df

  DE_Extract_FltH.df <- filter(DE_Extract.df, DE_Extract.df[,ThrSet[["LogFC"]][1]] >= ThrSet[["LogFC"]][2] &
                               DE_Extract.df[,ThrSet[["pVal"]][1]] < ThrSet[["pVal"]][2])

  DE_Extract_FltL.df <- filter(DE_Extract.df, DE_Extract.df[,ThrSet[["LogFC"]][1]] <= as.numeric(ThrSet[["LogFC"]][2])*(-1) &
                               DE_Extract.df[,ThrSet[["pVal"]][1]] < ThrSet[["pVal"]][2])
  DE_Extract_Flt.df <- rbind(DE_Extract_FltH.df,DE_Extract_FltL.df)

  DE_Extract_Flt.set <- rownames(DE_Extract_Flt.df)
  DE_Extract_FltH.set <- rownames(DE_Extract_FltH.df)
  DE_Extract_FltL.set <- rownames(DE_Extract_FltL.df)


  #### Export file ####
  write.table(DE_Extract.df, file = paste0(Save.Path,"/DEGAnalysis_",SampleName,"_",AnnoName,".tsv"),
              sep="\t", row.names= F, quote = FALSE)
  write.table(DE_Extract_Flt.df, file = paste0(Save.Path,"/DEGAnalysis_Flt_",SampleName,"_",AnnoName,".tsv"),
              sep="\t", row.names= F, quote = FALSE)


  #### Output ####
  Output <- list()
  Output[["DE_Extract.df"]] <- DE_Extract.df
  Output[["DE_Extract_Flt.df"]] <- DE_Extract_Flt.df
  Output[["DE_Extract_FltH.set"]] <- DE_Extract_FltH.set
  Output[["DE_Extract_FltL.set"]] <- DE_Extract_FltL.set
  Output[["Thr.lt"]] <- ThrSet

  return(Output)

}
