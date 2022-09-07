##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  #library(dplyr)
  #library(plyr) #merge.df<- rbind.fill(merge.df,new_1)

  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  if(!require("gtools")) install.packages("gtools")
  library(gtools)


##### Current path and new folder setting* #####
  InputFolder = "Customized_GSEAGenesets_Pathway3D"
  OutputFolder <- paste0(InputFolder,"_",Sys.Date(),"_CF") ## Generate output folder automatically
  dir.create(OutputFolder)

  ExpFilName <- "3DPathway" #"ComB" # Combine


##### Import files & Combine df #####
  # target.dir <- list.dirs(InputFolder)[-1]
  list.files <- list.files(InputFolder,full.names = T)
  Nfiles = length(list.files)

  for(i in 1:Nfiles){
    if(i==1){
      # Deal with different number of columns
      merge.df <- read.delim2(list.files[1],
                              col.names = 1:max(count.fields(list.files[1])),
                              header = F,sep = "\t")
    }else{
    new_1 <- read.delim2(paste0(list.files[i]),
                         col.names = 1:max(count.fields(list.files[i])),
                         header = F,sep = "\t")
    merge.df <- smartbind(merge.df,new_1)
    }

  }
  rm(new_1,i)

  #### Clean up df ####
  ## Remove duplicated
    merge.df <- merge.df[!duplicated(merge.df[,2]), ]

  # ## Remove NA (Have set in the write.table)
  # # Ref: https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
  #   merge.df[is.na(merge.df)] <- ""

  ##### Export Result of Combine #####
  ## Note ## Need to remove the quote
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName ,'.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

##### Filter by Keywords* #####
  ExpFilName_KW <- "EMT" # Export file name
  Keyword.lt <- list("EMT", c("trans","epithelial"))
  # Keyword.lt <- list("EMT", c("trans"))

  ## Filter and combine
    for(i in 1:length(Keyword.lt)){
      if(length(Keyword.lt[[i]])==1){
        merge_FLT_Temp.df <- merge.df[grepl(Keyword.lt[[i]],merge.df[,1], ignore.case=TRUE),]
      }else if(length(Keyword.lt[[i]])==2){
        merge_FLT_Temp.df <- merge.df[grepl(Keyword.lt[[i]][1],merge.df[,1], ignore.case=TRUE)
                                  & grepl(Keyword.lt[[i]][2],merge.df[,1], ignore.case=TRUE),]
      }else{
        merge_FLT_Temp.df <- merge.df[grepl(Keyword.lt[[i]][1],merge.df[,1], ignore.case=TRUE)
                                      & grepl(Keyword.lt[[i]][2],merge.df[,1], ignore.case=TRUE),]
        print(paste0("In",i,": Only the first 2 elements will be used"))
      }

      if(i==1){
        merge_FLT.df <- merge_FLT_Temp.df

      }else{
        merge_FLT.df <- smartbind(merge_FLT.df,merge_FLT_Temp.df)

      }
    }
    rm(i,merge_FLT_Temp.df)

  #### Clean up df ####
  ## Remove duplicated
  merge_FLT.df <- merge_FLT.df[!duplicated(merge_FLT.df[,2]), ]


  ##### Export Result #####
  ## Note ## Need to remove the quote
  write.table(merge_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName,'_',ExpFilName_KW ,'.txt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
  write.table(merge_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName,'_',ExpFilName_KW ,'.gmt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

  ### (pending) How to add conditions to a logical vector with a loop [r]
  ## Intersect all
  ## Ref: https://stackoverflow.com/questions/8817533/loop-of-a-loop-in-r
  ## Add conditions to a logical vector with a loop [r]
  ## https://stackoverflow.com/questions/40994881/add-conditions-to-a-logical-vector-with-a-loop-r

##### Choose specific Genesets* #####
  ExpFilName_SPEC = "SPEC" # Export file name

  Int_Path.set <- c(
    "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
    "REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY",
    "HALLMARK_E2F_TARGETS",
    "REACTOME_DNA_REPLICATION",
    "REACTOME_G2_M_CHECKPOINTS",
    "KEGG_OLFACTORY_TRANSDUCTION",
    "KEGG_CALCIUM_SIGNALING_PATHWAY"
  )

  merge_SPEC.df <- merge.df[merge.df[,1] %in% Int_Path.set,]

  ##### Export Result #####
  ## Note ## Need to remove the quote
  write.table(merge_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',ExpFilName,'_',ExpFilName_SPEC ,'.txt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
  write.table(merge_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',ExpFilName,'_',ExpFilName_SPEC ,'.gmt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")


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
