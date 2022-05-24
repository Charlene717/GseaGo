## 比較 StageI 高 比上 StageII 低

## 待修改
# Stage設定(If條件式)
# 匯入多個基因簡化與自動化
# 排列組合設定
# 縮寫簡化與通用化
# 備註說明
# 參數設定說明與匯出
# 檔名縮寫

## Clean variables
rm(list = ls())

######## Path setting ########
setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy

FileNameGE <- c("Xena_TCGA_PanCancer_GE.xena")
FileNamePheno1 <- c("Xena_TCGA_PanCancer_Pheno")
FileNamePheno2 <- c("Xena_TCGA_PanCancer_Pheno2.tsv")

######## Parameter setting ########
# Fold_Change = 1.5
Comb = 2
SDM <- 0 # Multiples of standard deviation
# FCName1 = trunc(Fold_Change)
# FCName2 = (Fold_Change-FCName1)*10
######## Filename setting ########
Target_gene_name_Multi <- paste0("RHOA_Comb",Comb,"_SDM",SDM,"_AjccSTGI&II_STGCompare")
# Target_gene_name_Multi <- paste0("RHOA_Comb",Comb,"_FC",FCName1,"p",FCName2,"_AjccSTGI&II_STGCompare")
# Comb:Combination # FC:Fold_Change # STG:Stage

Target_gene_name <- c("RHOA","ROCK1","ROCK2","FN1")
SetVersion <- c("_20210525")

FileName <-  c("Xena_TCGA_Pan")
ResultFolderName <- paste0("/",Target_gene_name_Multi,SetVersion) ## Generate output folder automatically
dir.create(paste0(PathName,ResultFolderName))

######## Import Raw Data ########

##Import gene expression data
GeneExp_Ori <- read.delim(paste0(PathName,"/",FileNameGE),  # 資料檔名
                          header=T,          # 資料中的第一列，作為欄位名稱
                          sep="\t")

## Import phenotype data
Phenotype1_Ori <- read.delim(paste0(PathName,"/",FileNamePheno1),  # 資料檔名
                            header=T,          # 資料中的第一列，作為欄位名稱
                            sep="\t")

Phenotype2_Ori <- read.delim(paste0(PathName,"/",FileNamePheno2),  # 資料檔名
                            header=T,          # 資料中的第一列，作為欄位名稱
                            sep="\t")

#*************************************** Note ***************************************#
# paste0 ==> concatenate strings without any separation/delimiter
# paste("Hello", "World", sep = "-") ==> concatenate strings with seperator "-"
#*************************************** Note ***************************************#


######## Data sorting ########
GeneExp <- GeneExp_Ori 
GeneExp[,1] <- as.character(GeneExp[,1])

GeneExp <- GeneExp[!duplicated(GeneExp$sample),]
# GeneExp[!duplicated(GeneExp[,c('sample')]),] 
# GeneExp[!duplicated(GeneExp[,1]),] 
# library("data.table") 
# setDT(GeneExp)[, .SD[1], by = .(sample)] 

# unique(GeneExp, by=c("sample"))

row.names(GeneExp) <- GeneExp[,1]
#GeneExp <- GeneExp[1:length(GeneExp[,1]), 2:length(GeneExp[1,])]
GeneExp <- GeneExp[, -1]

# colnames(GeneExp) <- GeneExp[1,]
# GeneExp <- GeneExp[-1, ]

######## Filter Data by Phenotype ########
PT = "Primary Tumor"
STN = "Solid Tissue Normal"

##### Phenotype_StageI #####
Phenotype_StageI <-  Phenotype1_Ori[Phenotype1_Ori$ajcc_pathologic_tumor_stage %in% c("Stage I"),]
Phenotype_StageI$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_StageI$X_PATIENT)
Phenotype_StageI$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_StageI$sample)
### Phenotype_StageI_PT <- Phenotype_StageI[Phenotype_StageI$sample_type %in% PT,]

##### Phenotype_StageII #####
Phenotype_StageII <-  Phenotype1_Ori[Phenotype1_Ori$ajcc_pathologic_tumor_stage %in% c("Stage II"),]
Phenotype_StageII$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_StageII$X_PATIENT)
Phenotype_StageII$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_StageII$sample)
### Phenotype_StageII_PT <- Phenotype_StageII[Phenotype_StageII$sample_type %in% PT,]


##### Phenotype Primary Tumor #####
PT = "Primary Tumor"
Phenotype_PT <-  Phenotype2_Ori[Phenotype2_Ori$sample_type %in% PT,]
Phenotype_PT$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_PT$sample)

##### Phenotype Solid Tissue Normal #####
STN = "Solid Tissue Normal"
Phenotype_STN <-  Phenotype2_Ori[Phenotype2_Ori$sample_type %in% STN,]
Phenotype_STN$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_STN$sample)

##### GE Stage #####
GeneExp_StageI <-  GeneExp[,colnames(GeneExp) %in% Phenotype_StageI$sample]
GeneExp_StageII <-  GeneExp[,colnames(GeneExp) %in% Phenotype_StageII$sample]

##### GE Stage & Primary Tumor #####
GeneExp_StageI_PT <-  GeneExp_StageI[,colnames(GeneExp_StageI) %in% Phenotype_PT$sample]
GeneExp_StageII_PT <-  GeneExp_StageII[,colnames(GeneExp_StageII) %in% Phenotype_PT$sample]


##### Target_gene ##### 
library(data.table)

# Extract data with Target_gene_name
Target_gene_StageI_PT <- GeneExp_StageI_PT[Target_gene_name,]
Target_gene_StageII_PT <- GeneExp_StageII_PT[Target_gene_name,]

library(dplyr) # Basic data manupilation tools
Target_gene_StageI_PT_Mean <- rowMeans(data.matrix(Target_gene_StageI_PT))
Target_gene_StageII_PT_Mean <- rowMeans(data.matrix(Target_gene_StageII_PT))

# Target_gene_SD <- sd(data.matrix(Target_gene))
library(matrixStats)
Target_gene_StageI_PT_SD <- rowSds(data.matrix(Target_gene_StageI_PT))
Target_gene_StageII_PT_SD <- rowSds(data.matrix(Target_gene_StageII_PT))

# Abbreviation of Mean and SD
StageI_MpSD <- as.data.frame(Target_gene_StageI_PT_Mean + SDM*Target_gene_StageI_PT_SD)
colnames(StageI_MpSD) <- "StageI_MpSD"
StageI_MmSD <- as.data.frame(Target_gene_StageI_PT_Mean - SDM*Target_gene_StageI_PT_SD)
colnames(StageI_MmSD) <- "StageI_MmSD"

StageII_MpSD <- as.data.frame(Target_gene_StageII_PT_Mean + SDM*Target_gene_StageII_PT_SD)
colnames(StageII_MpSD) <- "StageII_MpSD"
StageII_MmSD <- as.data.frame(Target_gene_StageII_PT_Mean - SDM*Target_gene_StageII_PT_SD)
colnames(StageII_MmSD) <- "StageII_MmSD"


# GeneExp_High <- GeneExp[,GeneExp[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD]
# GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD]


###### Try different Combination ######
if(Comb==1){          # c("RHOA")
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[1],] >= StageI_MpSD[Target_gene_name[1],]]
  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[1],] >= StageII_MmSD[Target_gene_name[1],]]
                                                
}else if(Comb==2){   # c("RHOA","ROCK1","ROCK2","FN1")
  
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[1],] >= StageI_MpSD[Target_gene_name[1],]&
                                                GeneExp_StageI_PT[Target_gene_name[2],] >= StageI_MpSD[Target_gene_name[2],]&
                                                GeneExp_StageI_PT[Target_gene_name[3],] >= StageI_MpSD[Target_gene_name[3],]&
                                                GeneExp_StageI_PT[Target_gene_name[4],] >= StageI_MpSD[Target_gene_name[4],]]
  
  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[1],] >= StageII_MmSD[Target_gene_name[1],]&
                                                GeneExp_StageII_PT[Target_gene_name[2],] >= StageII_MmSD[Target_gene_name[2],]&
                                                GeneExp_StageII_PT[Target_gene_name[3],] >= StageII_MmSD[Target_gene_name[3],]&
                                                GeneExp_StageII_PT[Target_gene_name[4],] >= StageII_MmSD[Target_gene_name[4],]]
  
}else if(Comb==3){  # c("RHOA","ROCK1","ROCK2")
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[1],] >= StageI_MpSD[Target_gene_name[1],]&
                                                GeneExp_StageI_PT[Target_gene_name[2],] >= StageI_MpSD[Target_gene_name[2],]&
                                                GeneExp_StageI_PT[Target_gene_name[3],] >= StageI_MpSD[Target_gene_name[3],]]
  
  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[1],] >= StageII_MmSD[Target_gene_name[1],]&
                                                GeneExp_StageII_PT[Target_gene_name[2],] >= StageII_MmSD[Target_gene_name[2],]&
                                                GeneExp_StageII_PT[Target_gene_name[3],] >= StageII_MmSD[Target_gene_name[3],]]
                                                 
}else if(Comb==4){  # c("RHOA","ROCK1")
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[1],] >= StageI_MpSD[Target_gene_name[1],]&
                                                GeneExp_StageI_PT[Target_gene_name[2],] >= StageI_MpSD[Target_gene_name[2],]]
  
  
  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[1],] >= StageII_MmSD[Target_gene_name[1],]&
                                                 GeneExp_StageII_PT[Target_gene_name[2],] >= StageII_MmSD[Target_gene_name[2],]]
}else if(Comb==5){  # c("RHOA","ROCK2")
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[1],] >= StageI_MpSD[Target_gene_name[1],]&
                                                GeneExp_StageI_PT[Target_gene_name[3],] >= StageI_MpSD[Target_gene_name[3],]]
  
  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[1],] >= StageII_MmSD[Target_gene_name[1],]&
                                                 GeneExp_StageII_PT[Target_gene_name[3],] >= StageII_MmSD[Target_gene_name[3],]]

}else if(Comb==6){  # c("ROCK1","FN1")
  GeneExp_StageI_PT_High <- GeneExp_StageI_PT[,GeneExp_StageI_PT[Target_gene_name[2],] >= StageI_MpSD[Target_gene_name[2],]&
                                                GeneExp_StageI_PT[Target_gene_name[4],] >= StageI_MpSD[Target_gene_name[4],]]

  GeneExp_StageII_PT_Low <- GeneExp_StageII_PT[,GeneExp_StageII_PT[Target_gene_name[2],] >= StageII_MmSD[Target_gene_name[2],]&
                                                 GeneExp_StageII_PT[Target_gene_name[4],] >= StageII_MmSD[Target_gene_name[4],]]
  
}  else{            # c("RHOA","ROCK2")
  "Error Comb!"
  
}



GeneExp_Group1 <- GeneExp_StageI_PT_High
GeneExp_Group2 <- GeneExp_StageII_PT_Low 


###### Count the sample numbers ######
GeneExp_Group1_Num <- length(GeneExp_Group1[1,])
GeneExp_Group2_Num <- length(GeneExp_Group2[1,])
GeneExp_Gene_Num <- length(GeneExp[,1])
Sample_Num <- GeneExp_Group1_Num + GeneExp_Group2_Num



GeneExp_Group1_List <- matrix(c(0), nrow =1, ncol = GeneExp_Group1_Num)
GeneExp_Group2_List <- matrix(c(1), nrow =1, ncol = GeneExp_Group2_Num)


GeneExp_Sum <- cbind(NAME=row.names(GeneExp),Description=matrix(c("na"), nrow =GeneExp_Gene_Num, ncol =1),GeneExp_Group1,GeneExp_Group2)
GSEA_GeneExp <- data.frame(t(colnames(GeneExp_Sum)), stringsAsFactors=FALSE)
colnames(GSEA_GeneExp) <- GSEA_GeneExp

#https://blog.csdn.net/sinat_26917383/article/details/50676894
library("plyr")  
GSEA_GeneExp <- rbind.fill(GSEA_GeneExp,GeneExp_Sum)


###### Export GSEA expression table ######
GSEASetting <- data.frame(NAME = c("#1.2",GeneExp_Gene_Num),Description = c('',Sample_Num))
GSEA_GeneExp<-rbind.fill(GSEASetting,GSEA_GeneExp)

write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.csv"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = ',')
write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.gct"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t')

###### Export GSEA Phenotype info ######
Pheno_Line1 <- c(Sample_Num,2,1)
Pheno_Line2 <- c(paste0("#",Target_gene_name_Multi,"_StageI_PT_High"),paste0(Target_gene_name_Multi,"_StageII_PT_Low"))
Pheno_Line3 <- c(GeneExp_Group1_List,GeneExp_Group2_List)
Pheno_sum <- rbind.fill(data.frame(t(Pheno_Line1)),data.frame(t(Pheno_Line2), stringsAsFactors=FALSE),data.frame(t(Pheno_Line3)))
write.table(Pheno_sum,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,".cls"),quote = FALSE,row.names = FALSE, na = "",col.names = FALSE)

###### Export parameter settings ######
