##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####     
  #library(dplyr)
  #library(plyr) #merge_1.df<- rbind.fill(merge_1.df,new_1)
  library(gtools) #merge_1.df<- smartbind(merge_1.df,new_1)
  library(tidyverse)

##### Current path and new folder setting* ##### 
  Res.Name = "GSEA_Geneset_Pathway_3Database"
  ResFld.Name <- paste0(Res.Name,"_",Sys.Date(),"_CF") ## Generate output folder automatically
  dir.create(ResFld.Name)

##### Import files & Combine df #####
  # target.dir <- list.dirs(Res.Name)[-1]
  list.files <- list.files(Res.Name,full.names = T)
  Nfiles = length(list.files)
  
  for(i in 1:Nfiles){
    if(i==1){
      # Deal with different number of columns
      merge_1.df <- read.delim2(list.files[1],
                             col.names = 1:max(count.fields(list.files[1])),
                             header = F,sep = "\t")
    }else{
    new_1 <- read.delim2(paste0(list.files[i]),
                         col.names = 1:max(count.fields(list.files[i])),
                         header = F,sep = "\t")
    merge_1.df <- smartbind(merge_1.df,new_1)
    }
    
  }
  rm(new_1,i)
  
  # #### Alternate ####
  #   ## Remove duplicated
  #     merge_1.df <- merge_1.df[!duplicated(merge_1.df[,2]), ]
  #   
  #   ## Remove NA
  #     # Ref: https://www.itranslater.com/qa/details/2097979270629950464
  #     merge_1.df <- merge_1.df[rowSums(is.na(merge_1.df))<length(merge_1.df),]
  #     merge_1.df <- merge_1.df[,colSums(is.na(merge_1.df))<nrow(merge_1.df)]
  #     # Ref: https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
       merge_1.df[is.na(merge_1.df)] <- ""

  ##### Export Result WithoutFilter #####
  ## Note ## Need to remove the quote
    write.table(merge_1.df,paste0(ResFld.Name,"/",Res.Name,'_WithoutFilter.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    write.table(merge_1.df,paste0(ResFld.Name,"/",Res.Name,'_WithoutFilter.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
  
##### Filter by Keywords* #####
  ## EMT
  merge_1_EMT1.df <- merge_1.df[grepl("EMT",merge_1.df[,1], ignore.case=TRUE),]
  merge_1_EMT2.df <- merge_1.df[grepl("trans",merge_1.df[,1], ignore.case=TRUE)
                                        & grepl("epithelial",merge_1.df[,1], ignore.case=TRUE),]
  merge_1_EMT.df <- smartbind(merge_1_EMT1.df,merge_1_EMT2.df)

  rm(merge_1_EMT1.df,merge_1_EMT2.df)
  
  ## DNA Repair
  merge_1_DNARepair.df <- merge_1.df[grepl("DNA",merge_1.df[,1], ignore.case=TRUE)
                               & grepl("Repair",merge_1.df[,1], ignore.case=TRUE),]
  
  ## Zinc
  merge_1_Zinc.df <- merge_1.df[grepl("Zinc",merge_1.df[,1], ignore.case=TRUE),]

  ## Methyl
  merge_1_Methyl.df <- merge_1.df[grepl("Methyl",merge_1.df[,1], ignore.case=TRUE),]

  ## Combine all index
  merge_1_AllIndex.df <- smartbind(merge_1_EMT.df,merge_1_Zinc.df,merge_1_DNARepair.df,merge_1_Methyl.df)
  merge_1_AllIndex.df <- merge_1_AllIndex.df[!duplicated(merge_1_AllIndex.df[,2]), ]
  
  ##### Export files WithFilter #####
    ## Note ## Need to remove the quote
    ## EMT
    write.table(merge_1_EMT.df,
                paste0(ResFld.Name,"/",Res.Name,'_EMT.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    ## Zinc
    write.table(merge_1_Zinc.df,
                paste0(ResFld.Name,"/",Res.Name,'_Zinc.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    ## DNA Repair
    write.table(merge_1_DNARepair.df,
                paste0(ResFld.Name,"/",Res.Name,'_DNARepair.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    ## Methyl
    write.table(merge_1_Methyl.df,
                paste0(ResFld.Name,"/",Res.Name,'_Methyl.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
  
    ## All index
    write.table(merge_1_AllIndex.df,
                paste0(ResFld.Name,"/",Res.Name,'_AllIndex.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    write.table(merge_1_AllIndex.df,
                paste0(ResFld.Name,"/",Res.Name,'_AllIndex.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
    
##### Filter multiple terms individually in different genesets  #####
    Filter.set <- c("Methyl","Zinc","AKT")
    SaveEach <- TRUE
  
    
    for (i in 1:length(Filter.set)) {
      if(i==1){
        merge_1_F.df <- merge_1.df[grepl(Filter.set[i],merge_1.df[,1], ignore.case=TRUE),]
        merge_1_AllIndex2.df <- merge_1_F.df
      }else{
        merge_1_F.df <- merge_1.df[grepl(Filter.set[i],merge_1.df[,1], ignore.case=TRUE),]
        merge_1_AllIndex2.df <- smartbind(merge_1_AllIndex2.df,merge_1_F.df)
      }
      
      # SaveEach
      if(SaveEach==TRUE){
        write.table(merge_1_F.df,
                    paste0(ResFld.Name,"/",Res.Name,'_',Filter.set[i],'.gmt'),
                    row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
      }
        
      rm(merge_1_F.df)
    }
    
      ## Remove duplicated
      merge_1_AllIndex2.df <- merge_1_AllIndex2.df[!duplicated(merge_1_AllIndex2.df[,2]), ]
      
      ##### Export files WithFilter #####
      ## All index
      write.table(merge_1_AllIndex2.df,
                  paste0(ResFld.Name,"/",Res.Name,'_AllIndex2.gmt'),
                  row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
      write.table(merge_1_AllIndex2.df,
                  paste0(ResFld.Name,"/",Res.Name,'_AllIndex2.txt'),
                  row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t')
      
##### Filter multiple terms individually or combined in different genesets  #####
      
            
#################################################################################  
  
   

# 
# ##### Other #####
# 
# ## Function setting
# 
# ## Call function
# filePath <- ""
# # Import R files in the same folder
# getFilePath <- function(fileName) {
#   # Absolute path of project folder
#   # path <- setwd("~")
#   path <- setwd(getwd()) 
#   # Combine strings without gaps
#   # <<- Assigning values to global variable
#   filePath <<- paste0(path ,"/" , fileName)  
#   # Load file
#   sourceObj <- source(filePath)
#   return(sourceObj)
# }
# getFilePath("HSsymbol2MMsymbol.R")
# 
# ## Geneset from GSEA
# H.all <- read.delim(paste0(PathName,"/h.all.v7.4.symbols.gmt"),header = F)
# 
# H.all.list <- list()
# for (i in c(1:length(H.all[,1]))) {
#   H.all.list.ori <- as.data.frame(t(H.all[i,3:length(H.all[i,])]))
#   colnames(H.all.list.ori)[[1]] <- c("Gene")
#   H.all.list.ori <- HSsymbol2MMsymbol(H.all.list.ori,"Gene")
#   
#   # Delete NA(or 0)
#   H.all.list.ori <- H.all.list.ori[H.all.list.ori$MM.symbol!=0,]
# 
#   H.all.list.ori <- unique(H.all.list.ori$MM.symbol)
#   H.all.list[[i]] <- as.character(H.all.list.ori)
#  
#   rm(H.all.list.ori)
#   names(H.all.list)[[i]] <- H.all[i,1]
# }
# 
# 
