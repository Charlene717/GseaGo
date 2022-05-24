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
Fold_Change = 1.5
Comb = 5

FCName1 = trunc(Fold_Change)
FCName2 = (Fold_Change-FCName1)*10
######## Filename setting ########
Target_gene_name_Multi <- paste0("RHOA_Comb",Comb,"_FC",FCName1,"p",FCName2,"_AjccSTGI&II_TvsN")
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
Phenotype_Stage <-  Phenotype1_Ori[Phenotype1_Ori$clinical_stage %in% c("Stage I","Stage II"),]
Phenotype_Stage <-  Phenotype1_Ori[Phenotype1_Ori$ajcc_pathologic_tumor_stage %in% c("Stage I","Stage II"),]
Phenotype_Stage2 <- Phenotype_Stage
Phenotype_Stage2$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage2$X_PATIENT)
Phenotype_Stage2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage2$sample)

Phenotype_Stage_All <- Phenotype1_Ori
# Try stage
Phenotype_Stage_All <- Phenotype_Stage2 
Phenotype_Stage_All$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage_All$X_PATIENT)
Phenotype_Stage_All$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage_All$sample)


PT = "Primary Tumor"
STN = "Solid Tissue Normal"
Phenotype_PT <-  Phenotype2_Ori[Phenotype2_Ori$sample_type %in% PT,]
Phenotype_PT2 <- Phenotype_PT
Phenotype_PT2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_PT2$sample)

Phenotype_STN <-  Phenotype2_Ori[Phenotype2_Ori$sample_type %in% STN,]
Phenotype_STN2 <- Phenotype_STN
Phenotype_STN2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_STN2$sample)

#
GeneExp_Stage_PT <-  GeneExp[,colnames(GeneExp) %in% Phenotype_PT2$sample]
GeneExp_Stage_STN <-  GeneExp[,colnames(GeneExp) %in% Phenotype_STN2$sample]
# GeneExp_Stage_STN2 <-  GeneExp_Stage_STN[,colnames(GeneExp_Stage_STN) %in% Phenotype_Stage2$sample]

# PSList_Test <- intersect(Phenotype_PT2$sample_type_id,Phenotype_STN2$sample_type_id)
# PSList <- intersect(Phenotype_PT2$sample_type_id,Phenotype_STN2$sample_type_id)

## 
Phenotype_STN_fromGE <- Phenotype_Stage_All[Phenotype_Stage_All$sample %in% colnames(GeneExp_Stage_STN),]
Phenotype_PT_fromGE <- Phenotype_Stage_All[Phenotype_Stage_All$sample %in% colnames(GeneExp_Stage_PT),]

PNList <- intersect(Phenotype_PT_fromGE$X_PATIENT,Phenotype_STN_fromGE$X_PATIENT)

Phenotype_STN_fromGE_intersect <- Phenotype_STN_fromGE[Phenotype_STN_fromGE$X_PATIENT %in% PNList,]
Phenotype_PT_fromGE_intersect <- Phenotype_PT_fromGE[Phenotype_PT_fromGE$X_PATIENT %in% PNList,]
Phenotype_Stage_All_PNintersect <- Phenotype_Stage_All[Phenotype_Stage_All$X_PATIENT %in% PNList, ]
  
# Final Normal and PrimaryTumor
GeneExp_Stage_PT2 <-  GeneExp_Stage_PT[,colnames(GeneExp_Stage_PT) %in% Phenotype_PT_fromGE_intersect$sample]
GeneExp_Stage_STN2 <-  GeneExp_Stage_STN[,colnames(GeneExp_Stage_STN) %in% Phenotype_STN_fromGE_intersect$sample]


######## Find Gene
GeneExp_RhoA_PT <- GeneExp_Stage_PT2[rownames(GeneExp_Stage_PT2) %in% c("RHOA","ROCK1","ROCK2","FN1"),]
GeneExp_RhoA_PT[5,] <- colnames(GeneExp_RhoA_PT)
GeneExp_RhoA_PT <- as.data.frame(t(GeneExp_RhoA_PT))
colnames(GeneExp_RhoA_PT)[5] <- c("sample")

GeneExp_RhoA_STN <- GeneExp_Stage_STN2[rownames(GeneExp_Stage_STN2) %in% c("RHOA","ROCK1","ROCK2","FN1"),]
GeneExp_RhoA_STN[5,] <- colnames(GeneExp_RhoA_STN)
GeneExp_RhoA_STN <- as.data.frame(t(GeneExp_RhoA_STN))
colnames(GeneExp_RhoA_STN)[5] <- c("sample")

library("dplyr")
GeneExp_RhoA_PT <- left_join(GeneExp_RhoA_PT,Phenotype_Stage_All_PNintersect,by="sample")
GeneExp_RhoA_PT <- GeneExp_RhoA_PT[,1:6]
GeneExp_RhoA_STN <- left_join(GeneExp_RhoA_STN,Phenotype_Stage_All_PNintersect,by="sample")
GeneExp_RhoA_STN <- GeneExp_RhoA_STN[,1:6]
GeneExp_RhoA_All <- left_join(GeneExp_RhoA_PT,GeneExp_RhoA_STN,by="X_PATIENT")
GeneExp_RhoA_All$RHOA.x <- as.numeric(GeneExp_RhoA_All$RHOA.x)
GeneExp_RhoA_All$RHOA.y <- as.numeric(GeneExp_RhoA_All$RHOA.y)
GeneExp_RhoA_All$ROCK1.x <- as.numeric(GeneExp_RhoA_All$ROCK1.x)
GeneExp_RhoA_All$ROCK1.y <- as.numeric(GeneExp_RhoA_All$ROCK1.y)
GeneExp_RhoA_All$ROCK2.x <- as.numeric(GeneExp_RhoA_All$ROCK2.x)
GeneExp_RhoA_All$ROCK2.y <- as.numeric(GeneExp_RhoA_All$ROCK2.y)
GeneExp_RhoA_All$FN1.x <- as.numeric(GeneExp_RhoA_All$FN1.x)
GeneExp_RhoA_All$FN1.y <- as.numeric(GeneExp_RhoA_All$FN1.y)


###### Try different Combination ######
if(Comb==1){
  GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > Fold_Change*GeneExp_RhoA_All$RHOA.y, ]
}else if(Comb==2){
  GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > Fold_Change*GeneExp_RhoA_All$RHOA.y &
                                                  GeneExp_RhoA_All$ROCK1.x > Fold_Change*GeneExp_RhoA_All$ROCK1.y &
                                                  GeneExp_RhoA_All$ROCK2.x > Fold_Change*GeneExp_RhoA_All$ROCK2.y &
                                                  GeneExp_RhoA_All$FN1.x > Fold_Change*GeneExp_RhoA_All$FN1.y,]
}else if(Comb==3){
  GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > Fold_Change*GeneExp_RhoA_All$RHOA.y &
                                                    GeneExp_RhoA_All$ROCK1.x > Fold_Change*GeneExp_RhoA_All$ROCK1.y &
                                                    GeneExp_RhoA_All$ROCK2.x > Fold_Change*GeneExp_RhoA_All$ROCK2.y,]
}else if(Comb==4){
  GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > Fold_Change*GeneExp_RhoA_All$RHOA.y &
                                                    GeneExp_RhoA_All$ROCK1.x > Fold_Change*GeneExp_RhoA_All$ROCK1.y,]
}  else{  
  GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > Fold_Change*GeneExp_RhoA_All$RHOA.y &
                                                    GeneExp_RhoA_All$ROCK2.x > Fold_Change*GeneExp_RhoA_All$ROCK2.y,]
}


GeneExp_Stage_PT_Final <-  GeneExp_Stage_PT2[,colnames(GeneExp_Stage_PT2) %in% GeneExp_RhoA_All_Candidates$sample.x]
GeneExp_Stage_STN_Final <-  GeneExp_Stage_STN2[,colnames(GeneExp_Stage_STN2) %in% GeneExp_RhoA_All_Candidates$sample.y]


GeneExp_PT <- GeneExp_Stage_PT_Final 
GeneExp_STN <- GeneExp_Stage_STN_Final 


###### Count the sample numbers ######
GeneExp_PT_Num <- length(GeneExp_PT[1,])
GeneExp_STN_Num <- length(GeneExp_STN[1,])
GeneExp_Gene_Num <- length(GeneExp[,1])
Sample_Num <- GeneExp_PT_Num+GeneExp_STN_Num



GeneExp_PT_Group <- matrix(c(0), nrow =1, ncol =GeneExp_PT_Num)
GeneExp_STN_Group <- matrix(c(1), nrow =1, ncol =GeneExp_STN_Num)


GeneExp_Sum <- cbind(NAME=row.names(GeneExp),Description=matrix(c("na"), nrow =GeneExp_Gene_Num, ncol =1),GeneExp_PT,GeneExp_STN)


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
Pheno_Line2 <- c(paste0("#",Target_gene_name_Multi,"_PT"),paste0(Target_gene_name_Multi,"_STN"))
Pheno_Line3 <- c(GeneExp_PT_Group,GeneExp_STN_Group)
Pheno_sum <- rbind.fill(data.frame(t(Pheno_Line1)),data.frame(t(Pheno_Line2), stringsAsFactors=FALSE),data.frame(t(Pheno_Line3)))
write.table(Pheno_sum,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,".cls"),quote = FALSE,row.names = FALSE, na = "",col.names = FALSE)

###### Export parameter settings ######
