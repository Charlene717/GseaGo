## Build files for GSEA official input

FUN_Group_GE = function(GeneExp.df, Anno.df,
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
  if(GroupSet$GeneExpMode == "Mean"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*0]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] < TarGene_Mean-TarGene_SD*0]

  }else if(GroupSet$GeneExpMode == "Mean1SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*1]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*1]

  }else if(GroupSet$GeneExpMode == "Mean2SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*2]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*2]

  }else if(GroupSet$GeneExpMode == "Mean3SD"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Mean+TarGene_SD*3]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Mean-TarGene_SD*3]

  }else if(GroupSet$GeneExpMode == "Median"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Q[3]]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] < TarGene_Q[3]]

  }else if(GroupSet$GeneExpMode == "Quartile"){
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] >= TarGene_Q[4]]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[TarGeneName,] <= TarGene_Q[2]]

  }else if(GroupSet$GeneExpMode == "Customize"){
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
  colnames(GeneExpAnno.df) <- c(colnames(Anno.df)[1], TarGene_name)

  AnnoNew.df <- left_join(Anno.df, GeneExpAnno.df)


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

