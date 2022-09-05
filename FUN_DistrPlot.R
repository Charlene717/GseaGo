## Build files for GSEA official input

FUN_DistrPlot = function(GeneExp.df,
                         TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                         Save.Path = Save.Path, ExportName = ExportName
){

  ##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("patchwork")) install.packages("patchwork")
  if(!require("eoffice")) install.packages("eoffice")

  library(tidyverse)
  library(patchwork)
  library(eoffice)

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


  ##### Basic DistPlt Function ####
  ## reshape df
  data <- reshape2::melt(GeneExp.df[TarGeneName,] %>% as.numeric())

  ## Set the color
  Custom.clr <- list(rect="#ffd5b5", line="#c95f22",text="#c95f22")

  ## Line.Set
  Line1V = GroupSet[["LowerCutoff"]]
  Line2V = GroupSet[["LowerCutoff"]]
  Line3V = GroupSet[["UpCutoff"]]

  DistPlt_Ori <- function(data,Line1V,Line2V,Line3V,Custom.clr,TarGene = TarGeneName ,Text_Basic.set = c("L1","L2","L3")) {
    TGeneDen.p <- ggplot(data,aes(value,fill=value, color=value)) +
      xlab("Expression level") +
      geom_density(alpha = 0.6, fill = "lightgray") +
      geom_rug() + theme_bw()

    ## Plot Mean and SD
    TGeneDen_SD.p <- ggPlot_vline(TGeneDen.p,
                                  data,
                                  Line.clr = Custom.clr,
                                  Line1 = Line1V,
                                  Line2 = Line2V,
                                  Line3 = Line3V,
                                  Text.set = Text_Basic.set)
    TGeneDen_SD.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
      labs(title= TarGene, x ="Expression level", y = "Density") -> TGeneDen_SD.p

    return(TGeneDen_SD.p)
  }

  DistPlt_Ori(data,Line1V,Line2V,Line3V,Custom.clr)

  ##### Set group conditions ####
    if(GroupSet$GeneExpMode == "Mean1SD"){
      Line1V = TarGene_Mean+TarGene_SD
      Line2V = TarGene_Mean
      Line3V = TarGene_Mean-TarGene_SD
      Text.set = c("Mean+1SD","Mean","Mean-1SD")

    }else if(GroupSet$GeneExpMode == "Mean2SD"){
      Line1V = TarGene_Mean+2*TarGene_SD
      Line2V = TarGene_Mean
      Line3V = TarGene_Mean-2*TarGene_SD
      Text.set = c("Mean+2SD","Mean","Mean-2SD")

    }else if(GroupSet$GeneExpMode == "Mean3SD"){
      Line1V = TarGene_Mean+3*TarGene_SD
      Line2V = TarGene_Mean
      Line3V = TarGene_Mean-3*TarGene_SD
      Text.set = c("Mean+3SD","Mean","Mean-3SD")

    }else if(GroupSet$GeneExpMode == "Mean"){
      Line1V = TarGene_Mean
      Line2V = TarGene_Mean
      Line3V = TarGene_Mean
      Text.set = c("Mean","Mean","Mean")

    }else if(GroupSet$GeneExpMode == "Quartile"){
      Line1V = TarGene_Q[4]
      Line2V = TarGene_Q[3]
      Line3V = TarGene_Q[2]
      Text.set = c("Q3","Q2","Q1")

    }else if(GroupSet$GeneExpMode == "Median"){
      Line1V = TarGene_Q[3]
      Line2V = TarGene_Q[3]
      Line3V = TarGene_Q[3]
      Text.set = c("Q2","Q2","Q2")

    }else if(GroupSet$GeneExpMode == "Customize"){
      Line1V = Line1V
      Line2V = Line2V
      Line3V = Line3V
      Text.set = c("LHigh","LHigh","LLow")

    }else{
      Line1V = TarGene_Mean+TarGene_SD
      Line2V = TarGene_Mean
      Line3V = TarGene_Mean-TarGene_SD
      Text.set = c("Mean+1SD","Mean","Mean-1SD")
    }

    TGeneDenR.p <- DistPlt_Ori(data,Line1V,Line2V,Line3V,Custom.clr,Text_Basic.set = Text.set)
    TGeneDenR.p


  ##### Visualization #####
  ## https://www.jianshu.com/p/9e5b7ffcf80f

  ## Set the color
  Mean_SD.clr <- list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" )
  Mean_Q.clr <- list(rect="#abede1", line="#12705f",text="#12705f" )

  ## Plot Mean and SD
  TGeneDen_SD.p <- DistPlt_Ori(data,
                               Line1V = TarGene_Mean+TarGene_SD,
                               Line2V = TarGene_Mean,
                               Line3V = TarGene_Mean-TarGene_SD,
                               Mean_SD.clr,
                               Text_Basic.set = c("Mean+1SD","Mean","Mean-1SD"))
  ## Plot Quartiles
  TGeneDen_Q.p <- DistPlt_Ori(data,
                               Line1V = TarGene_Q[4],
                               Line2V = TarGene_Q[3],
                               Line3V = TarGene_Q[2],
                               Mean_Q.clr,
                               Text_Basic.set = c("Q3","Q2","Q1"))



  ## Plot Quartiles & Mean and SD
  TGeneDen_SD_Q.p <- ggPlot_vline(TGeneDen_SD.p,data,
                                  Line.clr = Mean_Q.clr,
                                  Line1 = TarGene_Q[4],
                                  Line2 = TarGene_Q[3],
                                  Line3 = TarGene_Q[2],
                                  Text.set = c("Q3","Q2","Q1"),
                                  Text.yPos = 0.35,
                                  rectP = list(xWidth=0.015, yminP=0.3, ymaxP=0.4,alpha=0.8)
  )

  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
    labs(title= TarGeneName, x ="Expression level", y = "Density") -> TGeneDen_SD_Q.p



  #### Export PDF ####
  pdf(
    file = paste0(Save.Path,"/DensityPlot_",ExportName,".pdf"),
    width = 10,  height = 8
  )
    print(TGeneDenR.p)
    print(TGeneDen_SD.p)
    print(TGeneDen_Q.p)
    print(TGeneDen_SD_Q.p)

  dev.off()

  # #### Export PPT ####
  # TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
  #                                     OL_Thick = 1.5) +
  #   labs(title= TarGeneName,
  #        x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p
  #
  # topptx(TGeneDen_SD_Q2.p,paste0(Save.Path,"/DensityPlot_",ExportName,"_",TarGeneName,".pptx"))
  #
  # rm(TGeneDen_SD_Q2.p)


  ##### Note #####
  ## Finding Peak Values For a Density Distribution
  # http://ianmadd.github.io/pages/PeakDensityDistribution.html
  which.max(density(data$value)$y)
  max(density(data$value)$y)

  ## Plot multiple gene

  ## Set Output
  Output <- list()
  Output[["TGeneDenR.p"]] <- TGeneDenR.p
  Output[["TGeneDen_SD.p"]] <- TGeneDen_SD.p
  Output[["TGeneDen_Q.p"]] <- TGeneDen_Q.p
  Output[["TGeneDen_SD_Q.p"]] <- TGeneDen_SD_Q.p

  return(Output)

}
