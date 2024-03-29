## Build files for GSEA official input

FUN_Group_GE = function(GeneExp.df, MetaData.df,
                        TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                        Save.Path = Save.Path, ExportName = ExportName
){

  ##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  # if(!require("patchwork")) install.packages("patchwork")
  # if(!require("eoffice")) install.packages("eoffice")

  library(tidyverse)
  # library(patchwork)
  # library(eoffice)



  ##### Extract Target gene and Statistics ####
  # Extract data with TarGeneName
  TarGene_Mean <- GeneExp.df[TarGeneName,] %>%
                  as.numeric() %>%
                  mean()

  # rowMeans(data.matrix(TarGene))
  TarGene_SD <- GeneExp.df[TarGeneName,] %>%
                as.numeric() %>%
                sd()

  # Quartile
  TarGene_Q <- GeneExp.df[TarGeneName,] %>%
               as.numeric() %>%
               quantile()


  ##### Group the expression matrix according to the expression level of Target gene ####
  if(GroupSet$GEGroupMode == "Mean"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*0]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] < TarGene_Mean-TarGene_SD*0]

  }else if(GroupSet$GEGroupMode == "Mean1SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*1]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*1]

  }else if(GroupSet$GEGroupMode == "Mean2SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*2]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*2]

  }else if(GroupSet$GEGroupMode == "Mean3SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*3]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*3]

  }else if(GroupSet$GEGroupMode == "Median"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Q[3]]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] < TarGene_Q[3]]

  }else if(GroupSet$GEGroupMode == "Quartile"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Q[4]]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Q[2]]

  }else if(GroupSet$GEGroupMode == "Customize"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= GroupSet$UpCutoff]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= GroupSet$LowerCutoff]

  }else{
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*0]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] < TarGene_Mean-TarGene_SD*0]

  }


  ##### Add new anno #####
  GeneExp_high.df <- data.frame(ID = GeneExp_high.set %>% as.data.frame(), TarGene = "High")
  GeneExp_low.df <- data.frame(ID = GeneExp_low.set %>% as.data.frame(), TarGene = "Low")
  GeneExpAnno.df <- rbind(GeneExp_high.df, GeneExp_low.df)
  colnames(GeneExpAnno.df) <- c(colnames(MetaData.df)[1], TarGeneName)

  AnnoNew.df <- left_join(MetaData.df, GeneExpAnno.df)


  ## Export tsv
  write.table( AnnoNew.df,
    file=paste0(Save.Path,"/AnnoNew_",ExportName,".tsv"),
    quote = FALSE,row.names = FALSE, na = "",col.names = TRUE,sep = '\t')


  ## Set Output
  Output <- list()
  Output[["GeneExp_high.set"]] <- GeneExp_high.set
  Output[["GeneExp_low.set"]] <- GeneExp_low.set
  Output[["AnnoNew.df"]] <- AnnoNew.df


  return(Output)

}

