## Clear variables
rm(list = ls())

## Initialization

setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy

FileName <- c("Xena_TCGA_SARC_GE")
Target_gene_name_Multi <- c("RHOA_related")
Target_gene_name <- c("RHOA","ROCK1","ROCK2","FN1")
SetVersion <- c("_20210512V1")

ResultFolderName <- paste0("/",Target_gene_name_Multi,SetVersion) ## Generate output folder automatically
dir.create(paste0(PathName,ResultFolderName))

## Import genetic data file
GeneExp_Ori <- read.table(paste0(PathName,"/",FileName),  # 資料檔名
                          header=T,          # 資料中的第一列，作為欄位名稱
                          sep="")
# GeneExp_Ori <- read.table(paste0(PathName,"/Xena_TCGA_LGG_GE"),  # 資料檔名
#                           header=F,          # 資料中的第一列，作為欄位名稱
#                           sep="")
###################################### Note ######################################
# paste0 ==> concatenate strings without any separation/delimiter
# paste("Hello", "World", sep = "-") ==> concatenate strings with seperator "-"
##################################################################################


GeneExp <- GeneExp_Ori 
row.names(GeneExp) <- GeneExp[,1]
#GeneExp <- GeneExp[1:length(GeneExp[,1]), 2:length(GeneExp[1,])]
GeneExp <- GeneExp[, -1]

# load package 'data.table' 
library(data.table)

# Extract data with Target_gene_name
Target_gene <- GeneExp[Target_gene_name,]

# load package 'dplyr'
library(dplyr) # Basic data manupilation tools
Target_gene_Mean <- rowMeans(data.matrix(Target_gene))
# Target_gene_SD <- sd(data.matrix(Target_gene))
library(matrixStats)
Target_gene_SD <- rowSds(data.matrix(Target_gene))

# GeneExp_High <- GeneExp[,GeneExp[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD]
# GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD]


# Test OK
# TTT <- as.data.frame(Target_gene_Mean+Target_gene_SD)
# TTT2 <- as.data.frame(Target_gene_Mean-Target_gene_SD)
# GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= TTT[Target_gene_name[1],]]
# GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name[1],] <= TTT2[Target_gene_name[1],]]

MpSD <- as.data.frame(Target_gene_Mean+Target_gene_SD)
MmSD <- as.data.frame(Target_gene_Mean-Target_gene_SD)
Mean <- as.data.frame(Target_gene_Mean)

# GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= MpSD[Target_gene_name[1],] 
#                         & GeneExp[Target_gene_name[2],] >= MpSD[Target_gene_name[2],] 
#                         & GeneExp[Target_gene_name[3],] >= MpSD[Target_gene_name[3],]
#                         & GeneExp[Target_gene_name[4],] >= MpSD[Target_gene_name[4],]]
# GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name[1],] <= MmSD[Target_gene_name[1],]]

GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= Mean[Target_gene_name[1],] 
                        & GeneExp[Target_gene_name[2],] >= Mean[Target_gene_name[2],] 
                        & GeneExp[Target_gene_name[3],] >= Mean[Target_gene_name[3],]
                        & GeneExp[Target_gene_name[4],] >= Mean[Target_gene_name[4],]]

GeneExp_Low<- GeneExp[,GeneExp[Target_gene_name[1],] <= Mean[Target_gene_name[1],] 
                        & GeneExp[Target_gene_name[2],] <= Mean[Target_gene_name[2],] 
                        & GeneExp[Target_gene_name[3],] <= Mean[Target_gene_name[3],]
                        & GeneExp[Target_gene_name[4],] <= Mean[Target_gene_name[4],]]

# Count the numbers
GeneExp_High_Num <- length(GeneExp_High[1,])
GeneExp_Low_Num <- length(GeneExp_Low[1,])
GeneExp_Gene_Num <- length(GeneExp[,1])
Sample_Num <- GeneExp_High_Num+GeneExp_Low_Num



GeneExp_High_Group <- matrix(c(0), nrow =1, ncol =GeneExp_High_Num)
GeneExp_Low_Group <- matrix(c(1), nrow =1, ncol =GeneExp_Low_Num)


GeneExp_Sum <- cbind(NAME=row.names(GeneExp),Description=matrix(c("na"), nrow =GeneExp_Gene_Num, ncol =1),GeneExp_High,GeneExp_Low)


GSEA_GeneExp <- data.frame(t(colnames(GeneExp_Sum)), stringsAsFactors=FALSE)
colnames(GSEA_GeneExp) <- GSEA_GeneExp

#https://blog.csdn.net/sinat_26917383/article/details/50676894
library("plyr")  
GSEA_GeneExp <- rbind.fill(GSEA_GeneExp,GeneExp_Sum)

#########################################


GSEASetting <- data.frame(NAME = c("#1.2",GeneExp_Gene_Num),Description = c('',Sample_Num))
GSEA_GeneExp<-rbind.fill(GSEASetting,GSEA_GeneExp)

write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.csv"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = ',')
write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.gct"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t')
#########################################


Pheno_Line1 <- c(Sample_Num,2,1)
Pheno_Line2 <- c(paste0("#",Target_gene_name_Multi,"_high"),paste0(Target_gene_name_Multi,"_Low"))
Pheno_Line3 <- c(GeneExp_High_Group,GeneExp_Low_Group)
Pheno_sum <- rbind.fill(data.frame(t(Pheno_Line1)),data.frame(t(Pheno_Line2), stringsAsFactors=FALSE),data.frame(t(Pheno_Line3)))
write.table(Pheno_sum,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,".cls"),quote = FALSE,row.names = FALSE, na = "",col.names = FALSE)
#########################################