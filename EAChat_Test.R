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

##### Function setting #####
## Call function
source("FUN_DistrPlot_GE.R")
source("FUN_Beautify_ggplot.R")
# source("FUN_Find_Markers.R")
# source("FUN_GSEA_LargeGeneSet.R")
# source("FUN_GSEA_ggplot.R")
source("FUN_ggPlot_vline.R")
source("FUN_GSEA_ANAL.R")
source("FUN_VolcanoPlot.R")

##### Import setting and data loading* #####
Rec_Time_Point.lt[["Input_Start_Time"]] <- Sys.time() # %>% as.character()

## Set Import path
SetInputPath_FOL <- "Input_TCGA"  # Input Folder Name
SetInput_GE <- "Xena_TCGA_LGG_GE"
SetInput_Meta <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

## Load Gene expression file
GeneExp.df <- read.table(paste0(SetInputPath_FOL,"/",SetInput_GE), header=T, row.names = 1, sep="\t")
colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
## Load Metadata file
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
if (!dir.exists(Save_Path)){dir.create(Save_Path)}   ## Create new folder

## -[] Add setting record

#************************************************************************************************************************#
##### Update the genename ####
Rec_Time_Point.lt[["Update_Genename_Start_Time"]] <- Sys.time() # %>% as.character()

## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
Set_UpdateGene <- "No"  # Set_UpdateGene <- c("Yes","No")

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
  ## Distribution Plot
  Plot.DistrPlot <- FUN_DistrPlot_GE(GeneExp.df,
                                     TarGeneName = Set_TarGene_name, GroupSet = Set_TarGene,
                                     Save.Path = Save_Path, ExportName = SetExport_Name)
  Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
  Plot.DistrPlot_SD_Q


}else if(Set_GroupMode == "GoupByPheno"){
  Plot.Barplot <- ggplot(Metadata.df, aes(x=as.factor(Metadata.df[,Set_GroupCond$GroupType]), fill=as.factor(Metadata.df[,Set_GroupCond$GroupType]))) + geom_bar()
  ## BarPlot
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
source("RUN_DEG_Analysis.R")
DEG_Extract.df <- DEG_ANAL.lt[["DEG_Extract.df"]]

#### Volcano Plot ####
source("FUN_VolcanoPlot.R")
Plot.Volcano <- FUN_VolcanoPlot(DEG_Extract.df,
                                DiffThr = list("logFC",-1,1),
                                StatsTestThr = list("PValue",0.05),
                                color = c(High = "#ef476f",Mid = "gray",Low = "#0077b6"),
                                ShowGeneNumPos = 7, ShowGeneNumNeg = 7,
                                SizePoint = 3,  SizeAxisTitle = 16, SizeAxisText = 14, SizeLableText = 5,
                                ThkFrameLine = 2 , ThkThrLine = 0.8)
Plot.Volcano

pdf(file = paste0(Save_Path,"/DEG_VolcanoPlot_",SetExport_Name,".pdf"),
    width = 8,  height = 8)
Plot.Volcano %>% print()

dev.off()



Rec_Time_Point.lt[["DEG_End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["DEG"]] <- Rec_Time_Point.lt[["DEG_End_Time"]] - Rec_Time_Point.lt[["DEG_Start_Time"]]




##### Run Enrichment analysis in R #####
Rec_Time_Point.lt[["GSEA_Start_Time"]] <- Sys.time() # %>% as.character()
#### Run GSEA ####
source("FUN_GSEA_ANAL.R")
GSEAGeneSet_Int.set <- c("REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
                         #"REACTOME_SEMA4D_MEDIATED_INHIBITION_OF_CELL_ATTACHMENT_AND_MIGRATION",
                         "REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY",
                         "HALLMARK_E2F_TARGETS",
                         "REACTOME_DNA_REPLICATION",
                         "REACTOME_G2_M_CHECKPOINTS",
                         "KEGG_OLFACTORY_TRANSDUCTION",
                         "KEGG_CALCIUM_SIGNALING_PATHWAY"
)

Result_GSEA.lt <- FUN_GSEA_ANAL(DEG_Extract.df, CMGeneSet = GSEAGeneset.df,
                                DefaultGeneSet = "C2", Species = Set_Species, # Speices type can check by msigdbr_species()
                                NumGenesetsPlt = 15,
                                TarGeneName = Set_TarGene_name,
                                ThrSet = Set_DEGThr.lt,
                                Save.Path = Save_Path, ExportName = SetExport_Name, AnnoName = "Path",
                                Keyword = "HALLMARK",
                                Int_Path =  GSEAGeneSet_Int.set,
                                pAdjustMethod = "BH",  # pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                nPerm = 100000,
                                minGSSize = 15, maxGSSize = 500)

Result_GSEA <- Result_GSEA.lt[["Result_GSEA"]]
rm(Result_GSEA.lt)

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




##### Save RData #####
Rec_Time_Point.lt[["Save_RData_Start_Time"]] <- Sys.time() # %>% as.character()

save.image(paste0(Save_Path,"/GseaGo_",SetExport_Name,".RData"))

Rec_Time_Point.lt[["Save_RData_End_Time"]] <- Sys.time()
Rec_Time_Spend.lt[["Save_RData"]] <- Rec_Time_Point.lt[["Save_RData_End_Time"]] - Rec_Time_Point.lt[["Save_RData_Start_Time"]]

##### Record #####
#### Record time log ####
Rec_Time_Point.lt[["END_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Total_Time"]] <- Rec_Time_Point.lt[["END_Time"]] - Rec_Time_Point.lt[["Start_Time"]]


# ## Redundant method
# Rec_time_diff <- Rec_Time_Point.lt[["END_Time"]] - Rec_Time_Point.lt[["Start_Time"]]
# Rec_SaveTime_diff <- Rec_Time_Point.lt[["Save_RData_End_Time"]] - Rec_Time_Point.lt[["Save_RData_Start_Time"]]
# write(paste(" Program total time：", as.character(Rec_time_diff), "mins\n",
#             "Save RData time：", as.character(Rec_SaveTime_diff), "mins"), file = paste0(Save_Path,"/Rec_time_log.txt"))

# # ## Bug
# # write(paste0(Rec_Time_Spend.lt[["Total_Time(Without save R.Data)"]]) , file = paste0(Save_Path,"/Rec_time_log.txt"))

#### Record parameter ####


#### ChatGPT ####
# 載入 openai 套件
library(openai)
library(httr)
library(jsonlite)

# 設置 ChatGPT API 金鑰
source("OpenAI_API_Ch.txt")


# 定義分類函數
classify_terms <- function(terms, prompt, api_key) {
  # 將輸入的Terms轉換成一個字串
  terms_str <- paste(terms, collapse = " ")

  # 構建與 ChatGPT 對話的輸入
  input_text <- paste(prompt, terms_str)

  response <- POST(
    "https://api.openai.com/v1/engines/davinci/completions",
    # "https://api.openai.com/v1/engines/text-davinci-002/completions", # 更新引擎為text-davinci-002
    add_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type" = "application/json"
    ),
    body = toJSON(list(prompt = input_text, max_tokens = 50), auto_unbox = TRUE), # 將 auto_unbox 設置為 TRUE
    encode = "json"
  )

  # # 在API调用之后
  # if (!is.null(response) && !is.null(status_code(response)) && status_code(response) == 429) {
  #   stop("You have exceeded your API quota.")
  # }

  content <- fromJSON(content(response, as="text"))
  content[["choices"]][["text"]]
  return(content[["choices"]][["text"]])
}




# 定義一個示例的提示
# prompt <- "Please group the following Terms into groups based on relevance (You need to show the classification and their terms together): " # "請將以下Terms進行分類："
# prompt <- "Can you help me group the following Terms into groups based on relevance?" # "請將以下Terms進行分類："
# prompt <- "Please group the following Terms into groups based on relevance (Please output the corresponding classification results according to the order): " # "請將以下Terms進行分類："
prompt <- "Please add the classification label into each term based on relevance: " # "請將以下Terms進行分類："

# 要分類的Terms
terms_to_classify <- c(
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_COAGULATION",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

# terms_to_classify <- GSEAGeneSet_MetaData.df[1000:1500,]$STANDARD_NAME

# 使用函數進行分類
classification_result <- classify_terms(terms_to_classify, prompt, api_key)

# 將分類結果整合到 GSEAGeneSet_MetaData.df 中
GSEAGeneSet_MetaData.df[1000:1500, "Classification"] <- classification_result



# #################################################################################
# # 定義一個函數來處理每個富集名詞及其富集資訊
# classify_terms_with_info <- function(terms_with_info, prompt, api_key) {
#   # 初始化一個向量來存儲分類結果
#   classification_results <- character(length(terms_with_info))
#
#   # 使用迴圈處理每個富集名詞及其富集資訊
#   for (i in 1:length(terms_with_info)) {
#     term_info <- terms_with_info[i]
#     # 構建與 ChatGPT 對話的輸入，包括富集名詞及其富集資訊
#     input_text <- paste(prompt, term_info)
#
#     response <- POST(
#       "https://api.openai.com/v1/engines/text-davinci-002/completions",
#       add_headers(
#         "Authorization" = paste("Bearer", api_key),
#         "Content-Type" = "application/json"
#       ),
#       body = toJSON(list(prompt = input_text, max_tokens = 50), auto_unbox = TRUE),
#       encode = "json"
#     )
#
#     content <- fromJSON(content(response, as="text"))
#     classification_results[i] <-content[["choices"]][["text"]]
#     # classification_results[i] <- content$choices[[1]]$text
#   }
#
#   return(classification_results)
# }
#
# # 定義一個示例的提示
# prompt <- "Please add the classification label into each term based on relevance: "
#
# # 要分類的Terms及其富集資訊
# terms_with_info_to_classify <- c(
#   "SAKAI_TUMOR_INFILTRATING_MONOCYTES_UP - Additional information about this term",
#   "SAKAI_CHRONIC_HEPATITIS_VS_LIVER_CANCER_DN - Additional information about this term",
#   "JEON_SMAD6_TARGETS_UP - Additional information about this term",
#   "WANG_PROSTATE_CANCER_ANDROGEN_INDEPENDENT - Additional information about this term",
#   "SUNG_METASTASIS_STROMA_DN - Additional information about this term",
#   "SHIN_B_CELL_LYMPHOMA_CLUSTER_2 - Additional information about this term"
# )
#
# # 使用函數進行分類
# classification_results <- classify_terms_with_info(terms_with_info_to_classify, prompt, api_key)
#
# # 最後，你可以將每個富集名詞的分類結果與原始詞匯合成一個數據框或其他結構進行分析。
#
#
#

##################################################################################

# # 定義一個函數來分類富集名詞
# classify_terms_into_categories <- function(terms, api_key) {
#   # 初始化一個命名的空數據框，用於存儲分類結果
#   result_df <- data.frame(Term = character(0), Category = character(0))
#
#   # 要分類的生物學過程類別
#   categories <- c(
#     "Cell cycle related",
#     "Gene regulation related",
#     "Cellular responses and immunity related",
#     "Metabolism related",
#     "Other biological processes"
#   )
#
#   # 使用迴圈處理每個富集名詞
#   for (term in terms) {
#     # 構建與 ChatGPT 對話的輸入，包括富集名詞
#     input_text <- paste("Please classify the following biological processes：", term)
#
#     response <- POST(
#       "https://api.openai.com/v1/engines/davinci/completions",
#       add_headers(
#         "Authorization" = paste("Bearer", api_key),
#         "Content-Type" = "application/json"
#       ),
#       body = toJSON(list(prompt = input_text, max_tokens = 50), auto_unbox = TRUE),
#       encode = "json"
#     )
#
#     content <- fromJSON(content(response, as="text"))
#     classification <- content[["choices"]][["text"]]
#
#     # 將富集名詞及其分類結果添加到結果數據框中
#     result_df <- rbind(result_df, data.frame(Term = term, Category = classification))
#   }
#
#   return(result_df)
# }
#
# # 要分類的富集名詞
# terms_to_classify <- c(
#   "HALLMARK_MITOTIC_SPINDLE",
#   "HALLMARK_G2M_CHECKPOINT",
#   "HALLMARK_E2F_TARGETS",
#   "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
#   "HALLMARK_TGF_BETA_SIGNALING",
#   "HALLMARK_MYC_TARGETS_V1",
#   "HALLMARK_MYC_TARGETS_V2",
#   "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
#   "HALLMARK_NOTCH_SIGNALING",
#   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   "HALLMARK_APOPTOSIS",
#   "HALLMARK_IL2_STAT5_SIGNALING",
#   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
#   "HALLMARK_ALLOGRAFT_REJECTION",
#   "HALLMARK_GLYCOLYSIS",
#   "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
#   "HALLMARK_PANCREAS_BETA_CELLS",
#   "HALLMARK_DNA_REPAIR",
#   "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
#   "HALLMARK_ANGIOGENESIS",
#   "HALLMARK_MTORC1_SIGNALING",
#   "HALLMARK_COAGULATION",
#   "HALLMARK_SPERMATOGENESIS"
# )
#
# # 使用函數進行分類
# classification_results <- classify_terms_into_categories(terms_to_classify, api_key)
#
# # 最後，你可以使用 classification_results 數據框來查看每個富集名詞的分類結果。
# classification_results %>% View()

# ###############################################################################
# # 定義一個函數來分類富集名詞並填入數據框中的相應欄位
# classify_terms_and_fill_columns <- function(terms, prompt, api_key, df) {
#   # 使用迴圈處理每個富集名詞
#   for (i in 1:length(terms)) {
#     term <- terms[i]
#     # 構建與 ChatGPT 對話的輸入，包括富集名詞和提示
#     input_text <- paste(prompt, term)
#
#     response <- POST(
#       "https://api.openai.com/v1/engines/davinci/completions",
#       add_headers(
#         "Authorization" = paste("Bearer", api_key),
#         "Content-Type" = "application/json"
#       ),
#       body = toJSON(list(prompt = input_text, max_tokens = 50), auto_unbox = TRUE),
#       encode = "json"
#     )
#
#     content <- fromJSON(content(response, as="text"))
#     classification <- content[["choices"]][["text"]]
#
#     # 填入數據框中的相應欄位
#     df[i, "Term"] <- term
#     df[i, "Category"] <- classification
#   }
#
#   return(df)
# }
#
# # 要分類的富集名詞
# terms_to_classify <- c(
#   "HALLMARK_MITOTIC_SPINDLE",
#   "HALLMARK_G2M_CHECKPOINT",
#   "HALLMARK_E2F_TARGETS",
#   "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
#   "HALLMARK_TGF_BETA_SIGNALING",
#   "HALLMARK_MYC_TARGETS_V1",
#   "HALLMARK_MYC_TARGETS_V2",
#   "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
#   "HALLMARK_NOTCH_SIGNALING",
#   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   "HALLMARK_APOPTOSIS",
#   "HALLMARK_IL2_STAT5_SIGNALING",
#   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
#   "HALLMARK_ALLOGRAFT_REJECTION",
#   "HALLMARK_GLYCOLYSIS",
#   "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
#   "HALLMARK_PANCREAS_BETA_CELLS",
#   "HALLMARK_DNA_REPAIR",
#   "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
#   "HALLMARK_ANGIOGENESIS",
#   "HALLMARK_MTORC1_SIGNALING",
#   "HALLMARK_COAGULATION",
#   "HALLMARK_SPERMATOGENESIS"
# )
#
# # 創建一個空數據框，用於存儲結果
# result_df <- data.frame(Term = character(length(terms_to_classify)), Category = character(length(terms_to_classify)))
#
# # 定義提示
# prompt <- "Please add the classification label into each term based on relevance: "
#
# # 使用函數進行分類並填入數據框
# result_df <- classify_terms_and_fill_columns(terms_to_classify, prompt, api_key, result_df)
#
# # 最後，result_df 中將包含每個富集名詞及其分類結果
# result_df %>% View()
