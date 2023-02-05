## MSigDB Collections from GSEA
## Human
## http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
## Mouse
## https://www.gsea-msigdb.org/gsea/msigdb/mouse_geneset_resources.jsp

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  if(!require("gtools")) install.packages("gtools")
  library(gtools)

  #library(dplyr)
  #library(plyr) #merge.df<- rbind.fill(merge.df,new_1)

##### Current path and new folder setting* #####
  OutputFileName <- "ComB"
  InputFolder = "Cust_GSEA_Genesets_Test"
  OutputFolder <- paste0("Input_Genesets/", InputFolder, "_", OutputFileName)
  dir.create(OutputFolder) ## Generate output folder

##### Import files & Combine df #####
  # target.dir <- list.dirs( paste0("Input_Genesets/", InputFolder) )[-1]
  list.files <- list.files(paste0("Input_Genesets/", InputFolder),full.names = T)
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
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

##### Filter by Keywords* #####
  OutputFileName_KW <- "EMT" # Export file name
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
        print(paste0("In",i,": Only the first 2 elements will be used"))   ## 可以嘗試用迴圈的方式 ##整體改成用Apply寫
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
  write.table(merge_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_KW ,'.txt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
  write.table(merge_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_KW ,'.gmt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

  ### (pending) How to add conditions to a logical vector with a loop [r]
  ## Intersect all
  ## Ref: https://stackoverflow.com/questions/8817533/loop-of-a-loop-in-r
  ## Add conditions to a logical vector with a loop [r]
  ## https://stackoverflow.com/questions/40994881/add-conditions-to-a-logical-vector-with-a-loop-r

##### Choose specific Genesets* #####
  OutputFileName_SPEC = "SPEC" # Export file name

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
  write.table(merge_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_SPEC ,'.txt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
  write.table(merge_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_SPEC ,'.gmt'),
              row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")


#################################################################################

##### Convert GeneName #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  source("FUN_HSsymbol2MMsymbol.R")

  # # Inport Geneset
  # Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
  #                            col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
  #                            header = F,sep = "\t")


  ## 有新版更好的作法

  # Convert Human gene to mouse
  Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
  for (i in 1:nrow(Pathway.all)) {
    #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
    colnames(PathwayN)="Test"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
    Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)

  }

  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
  Pathway.all.MM[is.na(Pathway.all.MM)] <- ""
  Pathway.all.MM[Pathway.all.MM == 0] <- ""

  #### Save RData ####
  save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))



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
