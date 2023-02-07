## MSigDB Collections from GSEA
## http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  if(!require("gtools")) install.packages("gtools")
  library(gtools)

##### Condition setting* #####
  Set_Species = "Hs" #"Hs","Mm" # Homo sapiens(HS);"Mus musculus"(Mm)
  Set_LoadGeneBy = c("Default","symb") # c("Default","symb"),c("Default","entrez"),"Customize"
  Set_UpdateGeneName = TRUE # FALSE

  OutputFileName_KW <- "EMT" # Export file name of key word(KY)
  Set_Fiter = TRUE # FALSE
  Set_Fiter_KW.lt = list("EMT",c("trans","epithelial"), c("trans","epithelial","GOBP")) # "Default"

  OutputFileName_CTGY <- "C2" # Export file name of Category(CTGY)
  Set_Fiter_CTGY <- "C2"  # "Default"


## -[] Add setting record



##### Current path and new folder setting* #####
  OutputFileName <- "ComB"

  InputFolder
  if(Set_LoadGeneBy[1] == "Customize"){
    InputFolder = "Cust_GSEA_Genesets_Test"

  }else if(Set_LoadGeneBy[1] == "Default"){
    InputFolder = "Input_Genesets"

  }

  OutputFolder <- paste0("Input_Genesets/", InputFolder, "_", OutputFileName)
  dir.create(OutputFolder) ## Generate output folder

##### Import files & Combine df #####
  #### Import default RData ####
  load(paste0("Input_Genesets/Genesets_Default.RData"))

  ## clean up objects
  if(Set_Species == "Hs"){
    rm(list = str_subset(objects(), pattern = "_Mm_"))
  }else if(Set_Species == "Mm"){
    rm(list = str_subset(objects(), pattern = "_Hs_"))
  }else{
    print("Please Set Species in Hs or Mm")
  }

  if(Set_LoadGeneBy[1] == "Default" && Set_LoadGeneBy[2] =="symb"){
    rm(list = str_subset(objects(), pattern = "_entrez_"))
  }else if(Set_LoadGeneBy[1] == "Default" && Set_LoadGeneBy[2] == "entrez"){
    rm(list = str_subset(objects(), pattern = "_symb_"))
  }else if(Set_LoadGeneBy[1] == "Customize"){
    rm(list = str_subset(objects(), pattern = "_entrez_"))
    rm(list = str_subset(objects(), pattern = "_symb_"))
  }

  assign("GSEAGeneSet_MetaData.df", get(str_subset(objects(), pattern = "_XML.df")))
  rm(list = str_subset(objects(), pattern = "_XML.df"))

  if(Set_LoadGeneBy[1] == "Default"){
    assign("GSEAGeneSet.df", get(str_subset(objects(), pattern = "_gmt.df")))
    rm(list = str_subset(objects(), pattern = "_gmt.df"))

  }else if(Set_LoadGeneBy[1] == "Customize"){
    #### Import Customization ####
    # target.dir <- list.dirs( paste0("Input_Genesets/", InputFolder) )[-1]
    FilesList.set <- list.files(paste0("Input_Genesets/", InputFolder),full.names = T) %>%
                     str_subset(., pattern = "\\.gmt$")

    Nfiles = length(FilesList.set)

    for(i in 1:Nfiles){
      if(i==1){
        # Deal with different number of columns
        GSEAGeneSet.df <- read.delim2(FilesList.set[1],
                                      col.names = 1:max(count.fields(FilesList.set[1])),
                                      header = F,sep = "\t")
      }else{
        new_1 <- read.delim2(paste0(FilesList.set[i]),
                             col.names = 1:max(count.fields(FilesList.set[i])),
                             header = F,sep = "\t")
        GSEAGeneSet.df <- smartbind(GSEAGeneSet.df,new_1)
      }

    }
    rm(new_1,i)
  }




  #### Clean up df ####
  ## Remove duplicated
    GSEAGeneSet.df <- GSEAGeneSet.df[!duplicated(GSEAGeneSet.df[,2]), ]

  # ## Remove NA (Have set in the write.table) # Ref: https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
  #   GSEAGeneSet.df[is.na(GSEAGeneSet.df)] <- ""

##### Update gene name ####


##### Export Result of Combine #####
  ## Note ## Need to remove the quote
    # write.table(GSEAGeneSet.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.txt'),
    #             row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(GSEAGeneSet.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

##### Filter by Keywords* #####
  OutputFileName_KW <- "EMT" # Export file name
  Keyword.lt <- list("EMT",c("trans","epithelial"), c("trans","epithelial","GOBP"))

    for(i in 1:length(Keyword.lt)){
      if(length(Keyword.lt[[i]])==1){
        GSEAGeneSet_FLT_Temp.df <- GSEAGeneSet.df[grepl(Keyword.lt[[i]],GSEAGeneSet.df[,1], ignore.case=TRUE),]
      }else if(length(Keyword.lt[[i]])==2){
        GSEAGeneSet_FLT_Temp.df <- GSEAGeneSet.df[grepl(Keyword.lt[[i]][1],GSEAGeneSet.df[,1], ignore.case=TRUE)
                                  & grepl(Keyword.lt[[i]][2],GSEAGeneSet.df[,1], ignore.case=TRUE),]
      }else if(length(Keyword.lt[[i]])==3){
        GSEAGeneSet_FLT_Temp.df <- GSEAGeneSet.df[grepl(Keyword.lt[[i]][1],GSEAGeneSet.df[,1], ignore.case=TRUE)
                                      & grepl(Keyword.lt[[i]][2],GSEAGeneSet.df[,1], ignore.case=TRUE)
                                      & grepl(Keyword.lt[[i]][3],GSEAGeneSet.df[,1], ignore.case=TRUE),]
      }else{
        GSEAGeneSet_FLT_Temp.df <- GSEAGeneSet.df[grepl(Keyword.lt[[i]][1],GSEAGeneSet.df[,1], ignore.case=TRUE)
                                      & grepl(Keyword.lt[[i]][2],GSEAGeneSet.df[,1], ignore.case=TRUE)
                                      & grepl(Keyword.lt[[i]][3],GSEAGeneSet.df[,1], ignore.case=TRUE),]
        print(paste0("In",i,": Only the first 3 elements will be used"))   ## 可以嘗試用條件式+迴圈的方式 ##整體改成用Apply寫
      }

      if(i==1){
        GSEAGeneSet_FLT.df <- GSEAGeneSet_FLT_Temp.df

      }else{
        GSEAGeneSet_FLT.df <- smartbind(GSEAGeneSet_FLT.df,GSEAGeneSet_FLT_Temp.df)

      }
    }
    rm(i,GSEAGeneSet_FLT_Temp.df)

  #### Clean up df ####
  ## Remove duplicated
  GSEAGeneSet_FLT.df <- GSEAGeneSet_FLT.df[!duplicated(GSEAGeneSet_FLT.df[,2]), ]


  ##### Export Result #####
  ## Note ## Need to remove the quote
    # write.table(GSEAGeneSet_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_KW ,'.txt'),
    #             row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(GSEAGeneSet_FLT.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_KW ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

  ### (pending) How to add conditions to a logical vector with a loop [r]
  ## Intersect all
  ## Ref: https://stackoverflow.com/questions/8817533/loop-of-a-loop-in-r
  ## Add conditions to a logical vector with a loop [r]
  ## https://stackoverflow.com/questions/40994881/add-conditions-to-a-logical-vector-with-a-loop-r


#################################################################################################################################


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

  GSEAGeneSet_SPEC.df <- GSEAGeneSet.df[GSEAGeneSet.df[,1] %in% Int_Path.set,]

  ##### Export Result #####
  ## Note ## Need to remove the quote
    # write.table(GSEAGeneSet_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_SPEC ,'.txt'),
    #             row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(GSEAGeneSet_SPEC.df, paste0(OutputFolder,"/",InputFolder,'_',OutputFileName,'_',OutputFileName_SPEC ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")



#### Save RData ####
  save.image(paste0("Input_Genesets/", InputFolder,".RData"))

#################################################################################
  #### TO-do list ####
  ## 加入其他篩選條件
  ## 創建Geneset by 自己的實驗或線上數據(文字型,matrix型)
  ## 內文文字篩選
  ## 更新基因名稱
  ## 條件篩選可以嘗試用條件式+迴圈的方式
  ## 整體改成用Apply寫



