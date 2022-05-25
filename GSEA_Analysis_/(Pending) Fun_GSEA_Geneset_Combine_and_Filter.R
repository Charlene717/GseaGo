GSEA_ComFliter = function(RVersion = "GSEA_ComFliter_20211122",
                          FolderName = c("GSEA_Geneset_Pathway") 
                          ){
## 寫法:送一串關鍵字，使用list跑迴圈，每一個都生成一個檔案，外加一個綜合版本
##### Current path and new folder setting #####
PathName = setwd(getwd())
dir.create(paste0(PathName,"/",RVersion))

##### Import files & Combine df #####
first_category_name = list.files(FolderName)            # list.files
dir = paste("/",FolderName,"/",first_category_name,sep="") 
n = length(dir)                                                     

n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)                                          

library(dplyr)
library(plyr) # merge_1<- rbind.fill(merge_1,new_1)
library(gtools) # merge_1<- smartbind(merge_1,new_1)



##### Use This #####
## (Try1) Build a large df
merge_1 <- read.delim2(paste0(PathName,dir[1]),col.names = 1:max(count.fields(paste0(PathName,dir[1]))) ,header = F,sep = "\t")

for(i in 2:n){         
  b=list.files(dir[i]) 
  n_sub[i]=length(b)   
  
 # for(j in 1:n_sub[i]){    
        new_1 <- read.delim2(paste0(PathName,dir[i]),col.names = 1:max(count.fields(paste0(PathName,dir[i]))) ,header = F,sep = "\t")
 # }
  merge_1 <- smartbind(merge_1,new_1)
  
}

unique_merge_1 <- merge_1[!duplicated(merge_1[,2]), ]

# https://www.itranslater.com/qa/details/2097979270629950464
unique_merge_1 <- unique_merge_1[rowSums(is.na(unique_merge_1))<length(unique_merge_1),]
unique_merge_1 <- unique_merge_1[,colSums(is.na(unique_merge_1))<nrow(unique_merge_1)]
# https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
unique_merge_1[is.na(unique_merge_1)] <- ""

# # https://zhuanlan.zhihu.com/p/48997338
# drop_Nas <- function (data,dim){
#        if (dim == 2){
#              na_flag <- apply(!is.na(data),2,sum) # use ! to inverse 1 (to 0 0) and 0 (to 1).
#              data <- data[,-which(na_flag == 0)]
#          }
#        else if (dim == 1){
#              na_flag <- apply(!is.na(data),1,sum)
#              data <- data[-which(na_flag == 0),]
#          }
#        else{
#              warning("dim can only equal to 1 and 2, 1 means row, 2 means column ")
#          }
#        return(data)
#    }
# 
# unique_merge_1 <- drop_Nas(unique_merge_1,2) 


write.table(unique_merge_1,paste0(PathName,'/',RVersion,'/',FolderName,'_AllIndexWithoutFilter.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')
##### Use This #####






## (Try2) Build a large list
GSEA.all.list <- list()
# GSEA.sub.list <- list()
for(i in 1:n){         
  b=list.files(dir[i]) 
  n_sub[i]=length(b)   
  
  # for(j in 1:n_sub[i]){    
    new_1 <- read.delim2(paste0(PathName,dir[i]),col.names =1:max(count.fields(paste0(PathName,dir[i]))) ,header = F,sep = "\t")

    GSEA.sub.list <- list()
    for(k in 1:length(new_1[,1])){ 
      GSEA.sub.list.ori <- as.data.frame(t(new_1[k,3:length(new_1[k,])]))
      colnames(GSEA.sub.list.ori)[[1]] <- c("Gene")
      GSEA.sub.list[[k]] <- as.character(GSEA.sub.list.ori)
      names(GSEA.sub.list)[[k]] <- new_1[k,1]
      rm(GSEA.sub.list.ori)
    }
      GSEA.all.list <- c(GSEA.all.list,GSEA.sub.list)
      rm(GSEA.sub.list,new_1)
      
 # }
  # merge_1<- rbind.fill(merge_1,new_1)
}

## Loading problem
Load.Test <- read.delim2("D:/Dropbox/##_GitHub/0-R/##_PHH_Lab/GSEA_Analysis/GSEA_Geneset/c1.all.v7.4.symbols.gmt",header = F,fill=TRUE)
Load.Test2 <- read.csv("D:/Dropbox/##_GitHub/0-R/##_PHH_Lab/GSEA_Analysis/GSEA_Geneset/c1.all.v7.4.symbols.gmt",header = F,fill=TRUE,sep=",")

##### Filter Keywords #####
## Filter List
GSEA.all.list.TTT <- GSEA.all.list[,grepl("neurons", GSEA.all.list$chr1p12, ignore.case=TRUE)]
neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]

## Filter df

# # unique_merge_1_TTT <- unique_merge_1[unique_merge_1[,1] %in% c("EMT"),]
# unique_merge_1_TTT <- unique_merge_1[grepl("EMT",unique_merge_1[,1], ignore.case=TRUE),]

# EMT
unique_merge_1_EMT1 <- unique_merge_1[grepl("EMT",unique_merge_1[,1], ignore.case=TRUE),]
unique_merge_1_EMT2 <- unique_merge_1[grepl("trans",unique_merge_1[,1], ignore.case=TRUE) 
                                      & grepl("epithelial",unique_merge_1[,1], ignore.case=TRUE),]
unique_merge_1_EMT <- smartbind(unique_merge_1_EMT1,unique_merge_1_EMT2)

rm(unique_merge_1_EMT1,unique_merge_1_EMT2)

# Zinc
unique_merge_1_Zinc <- unique_merge_1[grepl("Zinc",unique_merge_1[,1], ignore.case=TRUE),]

# DNA Repair
unique_merge_1_DNARepair <- unique_merge_1[grepl("DNA",unique_merge_1[,1], ignore.case=TRUE) 
                                      & grepl("Repair",unique_merge_1[,1], ignore.case=TRUE),]
# Methyl
unique_merge_1_Methyl <- unique_merge_1[grepl("Methyl",unique_merge_1[,1], ignore.case=TRUE),]

# Combine all index
unique_merge_1_AllIndex <- smartbind(unique_merge_1_EMT,unique_merge_1_Zinc,unique_merge_1_DNARepair,unique_merge_1_Methyl)

##### Export files #####

# write.table(unique_merge_1,paste0(PathName,'/',RVersion,'/',FolderName,'_merge.txt'),
#             row.names = FALSE,col.names= TRUE, sep = '\t')
# 
# write.table(unique_merge_1,paste0(PathName,'/',RVersion,'/',FolderName,'_merge.csv'),
#             row.names = FALSE,col.names= TRUE, sep = ',')

# EMT
write.table(unique_merge_1_EMT,paste0(PathName,'/',RVersion,'/',FolderName,'_EMT.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')

# Zinc
write.table(unique_merge_1_Zinc,paste0(PathName,'/',RVersion,'/',FolderName,'_Zinc.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')
# DNA Repair
write.table(unique_merge_1_DNARepair,paste0(PathName,'/',RVersion,'/',FolderName,'_DNARepair.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')
# Methyl
write.table(unique_merge_1_Methyl,paste0(PathName,'/',RVersion,'/',FolderName,'_Methyl.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')

# All index
write.table(unique_merge_1_AllIndex,paste0(PathName,'/',RVersion,'/',FolderName,'_AllIndex.txt'),
            row.names = FALSE,col.names= FALSE, sep = '\t')

##### Other #####

## Function setting

## Call function
filePath <- ""
# Import R files in the same folder
getFilePath <- function(fileName) {
  # Absolute path of project folder
  # path <- setwd("~")
  path <- setwd(getwd()) 
  # Combine strings without gaps
  # <<- Assigning values to global variable
  filePath <<- paste0(path ,"/" , fileName)  
  # Load file
  sourceObj <- source(filePath)
  return(sourceObj)
}
getFilePath("HSsymbol2MMsymbol.R")

## Geneset from GSEA
H.all <- read.delim(paste0(PathName,"/h.all.v7.4.symbols.gmt"),header = F)

H.all.list <- list()
for (i in c(1:length(H.all[,1]))) {
  H.all.list.ori <- as.data.frame(t(H.all[i,3:length(H.all[i,])]))
  colnames(H.all.list.ori)[[1]] <- c("Gene")
  H.all.list.ori <- HSsymbol2MMsymbol(H.all.list.ori,"Gene")
  
  # Delete NA(or 0)
  H.all.list.ori <- H.all.list.ori[H.all.list.ori$MM.symbol!=0,]

  H.all.list.ori <- unique(H.all.list.ori$MM.symbol)
  H.all.list[[i]] <- as.character(H.all.list.ori)
 
  rm(H.all.list.ori)
  names(H.all.list)[[i]] <- H.all[i,1]
}


return(P2)
}

