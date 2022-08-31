## Build files for GSEA official input

FUN_Group_GE = function(GeneExp.df, Anno.df,
                        TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                        Save.Path = Save.Path, SampleName = SampleName
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


  ## Set Output
  Output <- list()
  Output[["GeneExp_high.set"]] <- GeneExp_high.set
  Output[["GeneExp_low.set"]] <- GeneExp_low.set
  Output[["AnnoNew.df"]] <- AnnoNew.df


  return(Output)

}



# ##### Visualization #####
# ## https://www.jianshu.com/p/9e5b7ffcf80f
#
# data <- reshape2::melt(GeneExp.df[TarGeneName,]%>%
#                          as.numeric())
# TGeneDen.p <- ggplot(data,aes(value,fill=value, color=value)) +
#   xlab("Expression level") +
#   geom_density(alpha = 0.6, fill = "lightgray") +
#   geom_rug() + theme_bw()
#
# ## Set the color
# Mean_SD.clr <- list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" )
# Mean_Q.clr <- list(rect="#abede1", line="#12705f",text="#12705f" )
#
# ## Plot Mean and SD
# TGeneDen.p
# TGeneDen_SD.p <- ggPlot_vline(TGeneDen.p,data,
#                               Line1 = TarGene_Mean+TarGene_SD,
#                               Line2 = TarGene_Mean,
#                               Line3 = TarGene_Mean-TarGene_SD,)
# TGeneDen_SD.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
#   labs(title= TarGeneName, x ="Expression level", y = "Density") -> TGeneDen_SD.p
#
# ## Plot Quartiles
# TGeneDen_Q.p <- ggPlot_vline(TGeneDen.p,data,
#                              Line.clr = Mean_Q.clr,
#                              Line1 = TarGene_Q[2],
#                              Line2 = TarGene_Q[3],
#                              Line3 = TarGene_Q[4],
#                              Text.set = c("Q1","Q2","Q3"),
#                              rectP = list(xWidth=0.015, yminP=0.45, ymaxP=0.55,alpha=0.8)
# )
#
# TGeneDen_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
#   labs(title= TarGeneName, x ="Expression level", y = "Density") -> TGeneDen_Q.p
#
# ## Plot Quartiles & Mean and SD
# TGeneDen_SD_Q.p <- ggPlot_vline(TGeneDen_SD.p,data,
#                                 Line.clr = Mean_Q.clr,
#                                 Line1 = TarGene_Q[2],
#                                 Line2 = TarGene_Q[3],
#                                 Line3 = TarGene_Q[4],
#                                 Text.set = c("Q1","Q2","Q3"),
#                                 Text.yPos = 0.35,
#                                 rectP = list(xWidth=0.015, yminP=0.3, ymaxP=0.4,alpha=0.8)
# )
#
# TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) +
#   labs(title= TarGeneName, x ="Expression level", y = "Density") -> TGeneDen_SD_Q.p
#
#
#
#
# pdf(
#   file = paste0(Save.Path,"/",SampleName,"_",TarGeneName,"_DensityPlot.pdf"),
#   width = 10,  height = 8
# )
# print(TGeneDen_SD.p)
# print(TGeneDen_Q.p)
# print(TGeneDen_SD_Q.p)
#
# dev.off()
#
# # Export PPT
# TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
#                                     OL_Thick = 1.5) +
#   labs(title= TarGeneName,
#        x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p
#
# topptx(TGeneDen_SD_Q2.p,paste0(Save.Path,"/",SampleName,"_",TarGeneName,"_DensityPlot.pptx"))
#
# rm(TGeneDen_SD_Q2.p)
#
#
# ##### Note #####
# ## Finding Peak Values For a Density Distribution
# # http://ianmadd.github.io/pages/PeakDensityDistribution.html
# which.max(density(data$value)$y)
# max(density(data$value)$y)
#
# ## Plot multiple gene
